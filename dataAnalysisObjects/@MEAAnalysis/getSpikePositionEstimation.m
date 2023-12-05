function [data,hand]=getSpikePositionEstimation(obj,varargin)
% [data,h]=spikePositionEstimation(obj,varargin)
% Function purpose : Estimates spike position
% Last updated : 26/05/16
data=[];
obj.checkFileRecording;

%default parameters
parseObj = inputParser;
addParameter(parseObj,'smoothingDuration',0.5); %[ms]
addParameter(parseObj,'maxDistForElectrodes',150); %[uM]
addParameter(parseObj,'electrodePitch',[]); %[uM]
addParameter(parseObj,'useSameTimePoint',true); %if true takes the time point of max channel of every neuron and uses this time for all channels, if false checks for the max of each ch seperately
addParameter(parseObj,'dXY',1); %[uM]
addParameter(parseObj,'dZ',1); %[uM]
addParameter(parseObj,'dL',1); %[uM]
addParameter(parseObj,'dV0',1); %[uV]
addParameter(parseObj,'maxZ',200); %[uM]
addParameter(parseObj,'minMaxV0',[15 300]); %[uV] - for global minima
addParameter(parseObj,'minMaxL',[10 100]); %[uM] - for global minima
addParameter(parseObj,'maxPotAtCellBody',250);
addParameter(parseObj,'gridColor',[0.8 0.8 0.8]);
addParameter(parseObj,'cellDiameter',15); %[um]
addParameter(parseObj,'modelType','reducedModelOpt');%,'realModelOpt','globalMinima','reducedModelOpt';
addParameter(parseObj,'smoothSpikes',false);
addParameter(parseObj,'usePostProcessing',false);
addParameter(parseObj,'L0',0);
addParameter(parseObj,'V0',0);
addParameter(parseObj,'Xc',[]);
addParameter(parseObj,'Yc',[]);
addParameter(parseObj,'decayModel','XYZV');
addParameter(parseObj,'maxIter',5e6);
addParameter(parseObj,'useSpikeExtremum',false);
addParameter(parseObj,'LGrid',45); %decay constants to use as initial conditions in optimization
addParameter(parseObj,'V0Grid',15:5:250); %Potential at r=0 (from cell) to use as initial conditions in optimization
addParameter(parseObj,'Icell',10.^(-10:0.05:-7)); %current density in initial axon segment
addParameter(parseObj,'neuronIdentity',[]); %the identity of the neuron - 1-not classified, 2- inhibitory, 3-excitatory
addParameter(parseObj,'WF',[]); %the waveform to triangulate [nNeurons,nCh,nSpikeSamples]
addParameter(parseObj,'STWaveformFile',[],@isstr)

addParameter(parseObj,'plot3D',true);
addParameter(parseObj,'h3D',[]);
addParameter(parseObj,'plot3DGrid',true);
addParameter(parseObj,'plot3DElectrodeNumbers',true);
addParameter(parseObj,'plot3DAddScaleBar',true);
addParameter(parseObj,'plot3DTransparency',0.5);
addParameter(parseObj,'plot3DRemoveNeuronOutSideMEA',false);

addParameter(parseObj,'elecNumFontSize',5);
addParameter(parseObj,'plotTriangulationOnWaveforms',false);
addParameter(parseObj,'plot1NeuronPerPlotInTriangulation',false);
addParameter(parseObj,'plotFit',false);

addParameter(parseObj,'overwrite',false,@isnumeric);
addParameter(parseObj,'saveFileName',[],@isstr);
addParameter(parseObj,'onlyPlot',false,@isnumeric);
addParameter(parseObj,'inputParams',false,@isnumeric);

parseObj.parse(varargin{:});
if parseObj.Results.inputParams
    disp(parseObj.Results);
    return;
end
%evaluate all input parameters in workspace
for i=1:numel(parseObj.Parameters)
    eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
end

%make parameter structure
par=parseObj.Results;

%check if analysis was already done done
outputDataWhenFinished=0;

if isempty(saveFileName)
    saveFileName=obj.files.(mfilename);
end

if ~overwrite
    if exist(saveFileName,'file')
        if onlyPlot
            load(saveFileName);
            electrodePitch=par.electrodePitch;
        else
            if nargout==1
                data=load(saveFileName);
            else
                disp('Analysis results already exist for this method, use overwrite if needed');
            end
            return;
        end
    end
else
    if nargout==1
        outputDataWhenFinished=1;
    end
end

if isempty(WF)
    %populate grid sorter object if does not exist
    STfileExist=obj.checkFileRecording(obj.files.getSpikeTrigWF,'Missing spike triggered average data...run getSpikeTrigWF');
    if ~STfileExist & isempty(STWaveformFile)
        if isempty(obj.gridSorterObj)
            obj=obj.populateGridSorterObj;
        end
        gridSortFileExist=obj.checkFileRecording(obj.gridSorterObj.sortingFileNames.STWaveformFile,'Missing post processing in spike sorting, rerun sorting');
        if gridSortFileExist
            WF=load(obj.gridSorterObj.sortingFileNames.STWaveformFile,'avgHPWF'); %load data
            preMs=obj.gridSorterObj.postPreFilteredWindow;
            warning('Getting spike triggered data from grid sorter!!!! In the future, run getSpikeTrigWF from MEAAnalysis');
        else
            error('No spike triggered data in either file or gridsorter!!!');
        end
    else
        if isempty(STWaveformFile)
           STWaveformFile=obj.files.getSpikeTrigWF;
        end
        WF=load(STWaveformFile,'avgHPWF','par');
        preMs=WF.par.preFilteredWindow;
    end
else %use waveforms from an input WF
    preMs=WF.preMs;
    if ~isstruct(WF)
        WF2=WF;clear WF;
        WF.avgHPWF=WF2;
        clear WF2;
        disp('Taking pre spike duration from grid sorter object!!!');
        obj=obj.populateGridSorterObj;
        preMs=obj.gridSorterObj.postPreFilteredWindow;
    end
    disp('Using highpass waveforms from input!');
end

%get recording information
Fs=obj.currentDataObj.samplingFrequency;
ch=obj.currentDataObj.channelNumbers;
En=obj.currentDataObj.chLayoutNumbers;

%check if ch numbers vector was given and if not, take all channels in array in serial order
if isempty(ch)
    ch=unique(En(~isnan(En)));
    if numel(ch)~=nCh
        error('Number of channels in waveform does not match the channels in En');
    end
else
    if size(ch,1)==1
        ch=ch';
    end
end

%get electrode positions in physical space
if ~isempty(electrodePitch)
    [Xc,Yc]=obj.currentDataObj.getElectrodePositions(electrodePitch);
else
    [Xc,Yc]=obj.currentDataObj.getElectrodePositions;
    electrodePitch=obj.currentDataObj.electrodePitch;
end

%% Output list of default variables
hand=[];testData=[];
if isempty(WF.avgHPWF) %run on artificial data
    nNeurons=10;
    nCh=numel(unique(En(~isnan(En))));
    nSpikeSamples=preMs*2*Fs/1000;
    
    neuronIdentity=round(rand(1,nNeurons))+2;
    tX=rand(1,nNeurons)*100-50;
    tY=rand(1,nNeurons)*100-50;
    tZ=rand(1,nNeurons)*100;
    
    %tV0=V0+rand(1,nNeurons)*10;
    tV0=rand(1,nNeurons)*(minMaxV0(2)-minMaxV0(1))+minMaxV0(1);
    %tL=L0+rand(1,nNeurons)*50-25;
    tL=rand(1,nNeurons)*(minMaxL(2)-minMaxL(1))+minMaxL(1);
    tCh=ceil(rand(1,nNeurons)*nCh);
    tCh(tCh==0)=1;
    
    smoothingSamples=round(smoothingDuration*Fs/1000);
    preSpkSamples=preMs*Fs/1000;
    WF.avgHPWF=zeros([nNeurons,nCh,nSpikeSamples]);
    for i=1:nNeurons
        for j=1:nCh
            WF.avgHPWF(i,j,preSpkSamples:(preSpkSamples+3*smoothingSamples))=...
                -tV0(i).*exp(-sqrt(  (Xc(tCh(i))+tX(i)-Xc(j)).^2 + (Yc(tCh(i))+tY(i)-Yc(j)).^2+(tZ(i)-0).^2   )./tL(i));
            %WF.avgHPWF(i,j,preSpkSamples:(preSpkSamples+3*smoothingSamples))=...
            %    -tV0(i).*exp(-sqrt(  (Xc(tCh(i))+tX(i)-Xc(j)).^2 + (Yc(tCh(i))+tY(i)-Yc(j)).^2   )./tL(i));
        end
        %h=axes;activityTracePhysicalSpacePlot(h,1:nCh,squeeze(WF.avgHPWF(i,:,:)),En,'traceColor',[0.2 0.2 0.8],'scaling','none');
    end
    
    %save parameter for output
    testData.neuronIdentity=neuronIdentity;
    testData.tX=tX;
    testData.tY=tY;
    testData.tZ=tZ;
    testData.tV0=tV0;
    testData.tL=tL;
    testData.tCh=tCh;
end

% Find neuron with reliable spiking
[nNeurons,nCh,nSpikeSamples]=size(WF.avgHPWF);

% calculate pre and post samples and their position in the waveform
preSpkSamples=preMs*Fs/1000;
postSpkSamples=nSpikeSamples-preSpkSamples;


%% model with x,y,z,V
%according to the current grid arrangement positive Y is down on real MEA space and positive X is right on real MEA space

if useSpikeExtremum % find spike extremum ch and position
    [vMaxAll,pMaxAll]=max(abs(WF.avgHPWF),[],3);
    [vMax,pChMax]=max(vMaxAll,[],2);
else %find spike minimum
    [vMaxAll,pMaxAll]=min(WF.avgHPWF,[],3);
    [vMax,pChMax]=min(vMaxAll,[],2);
    vMax=abs(vMax);
    vMaxAll=abs(vMaxAll);
end
pMaxSampleInMaxCh=pMaxAll( sub2ind(size(pMaxAll), 1:nNeurons, pChMax') );

%smooth spike shapes with local linear regression
if smoothSpikes
    smoothingSamples=round(smoothingDuration*Fs/1000);
    for i=1:nNeurons
        for j=1:nCh
            WF.avgHPWF(i,j,:) = smooth(WF.avgHPWF(i,j,:),smoothingSamples,'loess');
        end
    end
end

if ~onlyPlot
    if strcmp(modelType,'globalMinima')
        %define search limits
        X=-electrodePitch:dXY:electrodePitch;
        Y=-electrodePitch:dXY:electrodePitch;
        Z=0:dZ:maxZ;
        L=minMaxL(1):dL:minMaxL(2);
        V0=minMaxV0(1):dV0:minMaxV0(2);
        
        %initialize arrays
        V0All=ones(numel(X),numel(Y),numel(Z),numel(L),numel(V0));
        XAll=V0All;
        YAll=V0All;
        ZAll=V0All;
        LAll=V0All;
        
        X5D(:,1,1,1,1)=X;
        Y5D(1,:,1,1,1)=Y;
        Z5D(1,1,:,1,1)=Z;
        L5D(1,1,1,:,1)=L;
        V05D(1,1,1,1,:)=V0;
        
        XAll=bsxfun(@times,XAll,X5D);
        YAll=bsxfun(@times,YAll,Y5D);
        ZAll=bsxfun(@times,ZAll,Z5D);
        LAll=bsxfun(@times,LAll,L5D);
        V0All=bsxfun(@times,V0All,V05D);
    end
    
    % main loop
    estX=zeros(1,nCh);
    estY=zeros(1,nCh);
    estV0=zeros(1,nCh);
    estZ=zeros(1,nCh);
    allVPeak=zeros(nNeurons,nCh);
    eValAll=zeros(1,nNeurons);
    
    hWait=waitbar(0,'Triangulating neurons...');
    for i=1:nNeurons
        waitbar(i/nNeurons,hWait);
        
        pXY=sqrt((Xc-Xc(pChMax(i))).^2+(Yc-Yc(pChMax(i))).^2)<maxDistForElectrodes;
        xTmp=Xc(pXY)-Xc(pChMax(i));
        yTmp=Yc(pXY)-Yc(pChMax(i));
        zTmp=zeros(size(yTmp));
        
        if useSameTimePoint
            V0tmp=abs(WF.avgHPWF(i,pXY,pMaxSampleInMaxCh(i)));
        else
            V0tmp=vMaxAll(i,pXY);
        end
        allVPeak(i,:)=WF.avgHPWF(i,:,pMaxSampleInMaxCh(i));
        
        if strcmp(modelType,'globalMinima')
            
            score=zeros(size(XAll));
            for k=1:numel(xTmp)
                score=score+(       V0tmp(k)-V0All.*exp(-sqrt((xTmp(k)-XAll).^2+(yTmp(k)-YAll).^2+(zTmp(k)-ZAll).^2)./LAll)         ).^2;
                %score=score+(log(V0tmp(k)./V0All)  +  sqrt((xTmp(k)-XAll).^2+(yTmp(k)-YAll).^2+(zTmp(k)-ZAll).^2)./L0).^2;
            end
            %figure;imagesc(Z,squeeze(V0),squeeze(score(find(X==0),find(X==0),:,:)));xlabel('Z');ylabel('V_0');colorbar
            
            [bestVal(i),bestPlace]=min(score(:));
            [bestXp,bestYp,bestZp,bestLp,bestV0p]=ind2sub(size(V0All),bestPlace);
            Xr(i)=X(bestXp);
            Yr(i)=Y(bestYp);
            Zr(i)=Z(bestZp);
            Lamda(i)=L(bestLp);
            V(i)=V0(bestV0p);
            
        elseif strcmp(modelType,'realModelOpt')
            
            eVal=Inf;
            for j=1:numel(Icell)
                [estTmpTmp,eValTmp] = fminsearch(@(x) potEst(x,xTmp*1e-6,yTmp*1e-6,V0tmp*1e-6),[0,0,0,Icell(j)],optimset('MaxFunEvals',1000000,'Display','off'));
                allEst(j)=eValTmp;
                if eValTmp<eVal
                    estTmp=estTmpTmp;
                    eVal=eValTmp;
                end
            end
            %plot(Icell,allEst,'.');set(gca,'XScale','log');pause;
            
            Xr(i)=estTmp(1);
            Yr(i)=estTmp(2);
            Zr(i)=estTmp(3);
            eVal(i)=eVal;
            I(i)=estTmp(4);

            
        elseif strcmp(modelType,'reducedModelOpt')
            
            eVal=Inf;
            for j=1:numel(V0Grid)
                for k=1:numel(LGrid)
                    switch decayModel
                        case 'XYZ'
                            %score=XYZ_Model(X,xTmp,yTmp,zTmp,V0tmp,V0,L,maxZ)
                            initCond=[0,0,0];
                            funcHand=@(x) XYZ_Model(x,xTmp,yTmp,zTmp,maxZ,V0tmp,V0Grid(j),LGrid(k));
                            pV=0;pL=0;
                        case 'XYZLV'
                            %score=XYZLV_Model(X,xTmp,yTmp,zTmp,V0tmp,vMax,maxPotAtCellBody,maxZ)
                            initCond=[0,0,0,LGrid(k),max(vMax(i),V0Grid(j))];
                            funcHand=@(x) XYZLV_Model(x,xTmp,yTmp,zTmp,maxZ,V0tmp,vMax(i),maxPotAtCellBody);
                            pV=5;pL=4;
                        case 'XYZL'
                            %score=XYZL_Model(X,xTmp,yTmp,zTmp,V0tmp,V0,maxZ)
                            initCond=[0,0,0,LGrid(k)];
                            funcHand=@(x) XYZL_Model(x,xTmp,yTmp,zTmp,maxZ,V0tmp,V0Grid(j));
                            pV=0;pL=4;
                        case 'XYZV'
                            %score=XYZV_Model(X,xTmp,yTmp,zTmp,V0tmp,L,vMax,maxPotAtCellBody,maxZ)
                            initCond=[0,0,0,max(vMax(i),V0Grid(j))];
                            funcHand=@(x) XYZV_Model(x,xTmp,yTmp,zTmp,maxZ,V0tmp,LGrid(k),vMax(i),maxPotAtCellBody);
                            pV=4;pL=0;
                        case 'powerDecay' %here V0 is the multiplier and L is the decay constant
                            %score=powerDecay_Model(X,xTmp,yTmp,zTmp,maxZ,V0tmp,V0,L)
                            initCond=[0,0,0];
                            funcHand=@(x) powerDecay_Model(x,xTmp,yTmp,zTmp,maxZ,V0tmp,V0Grid(j),LGrid(k));
                            pV=0;pL=0;
                        otherwise
                            error('The selected decay model is not supported! please choose from:XYZ,XYZL,XYZV,XYZLV,powerDecay');
                    end
                    
                    [estTmpTmp,eValTmp] = fminsearch(funcHand,initCond,optimset('MaxFunEvals',maxIter,'Display','off'));
                    
                    if eValTmp<eVal
                        estTmp=estTmpTmp;
                        eVal=eValTmp;
                        tmpV0=V0Grid(j);
                        tmpL=LGrid(k);
                    end
                end
            end
            
            eValAll(i)=eVal;
            Lamda(i)=tmpL;%Lamda(i)=estTmp(4);
            V(i)=tmpV0; %V(i)=estTmp(5);
            Xr(i)=estTmp(1);
            Yr(i)=estTmp(2);
            Zr(i)=estTmp(3);
            if pV>0
                V(i)=estTmp(pV);
            end
            if pL>0
                Lamda(i)=estTmp(pL);
            end

        end
        
        Xrs(i)=Xr(i)+Xc(pChMax(i));
        Yrs(i)=Yr(i)+Yc(pChMax(i));
        Zrs(i)=Zr(i);
        electrodeX=Xc;
        electrodeY=Yc;
    end
    close(hWait);
    
    save(saveFileName,'par','Xrs','Yrs','Zrs','Xr','Yr','Zr','Lamda','V','eValAll','electrodeX','electrodeY','pChMax','allVPeak');
    if outputDataWhenFinished
        data=load(saveFileName);
    end
    %notice that all variables after electrodeX are only required for future plotting
else
    load(saveFileName,'Xrs','Yrs','Zrs');
end




%{
figure;subplot(2,3,1);plot(tL-Lamda,'.');title('L');subplot(2,3,2);plot(tX-Xr,'.');title('X');subplot(2,3,3);plot(tY-Yr,'.');title('Y');subplot(2,3,4);plot(tZ-Zr,'.');title('Z');subplot(2,3,5);plot(tV0-V,'.');title('V0');
%}

if plotFit
    if nNeurons>=40
        nPlotsInRow=10;
        nPlotsInCol=6;
    else
        nPlotsInRow=8;
        nPlotsInCol=5;
    end
    hand.fFit(1)=figure('position',[40 80 1500 900]);
    c=1;
    for i=1:nNeurons
        iFig=mod(i-1,nPlotsInRow*nPlotsInCol)+1;
        subaxis(nPlotsInCol,nPlotsInRow,iFig,'S',0.005,'M',0.01);
        d=sqrt((Xrs(i)-electrodeX).^2+(Yrs(i)-electrodeY).^2+(Zrs(i)).^2);
        plot(d,vMaxAll(i,:),'.');hold on;
        if strcmp(par.decayModel,'powerDecay')
            plot(sort(d),V(i)./(sort(d)).^Lamda(i),'r');
        else
            plot(sort(d),V(i).*exp(-sort(d)/Lamda(i)),'r');
        end
        
        xlim([0 maxDistForElectrodes+100]);
        ylim([0 max(vMaxAll(i,:))+3]);
        text(0,0,[num2str(i) ' ,V0=' num2str(V(i),3) ' ,\lambda=' num2str(Lamda(i),3),' ,Z=' num2str(Zrs(i),3)],'HorizontalAlignment','left','verticalAlignment','Bottom');
        set(gca,'XTickLabel',[],'YTickLabel',[]);
        if mod(i,nPlotsInRow*nPlotsInCol)==0
            c=c+1;
            hand.fFit(c)=figure('position',[40 80 1500 900]);
        end
    end
end

if plot3D

    if isempty(neuronIdentity) %if no an empty neuron identity is given
        neuronIdentity=ones(1,nNeurons); %set all neuorns with unindentified polarity
    else
        neuronIdentity(neuronIdentity==0)=1;
    end
    
    %create electrode grid with numbers for plotting
    [max_grid_y,max_grid_x]=size(En);
    XL=[0.5 max_grid_x+0.5]*electrodePitch;
    YL=[0.5 max_grid_y+0.5]*electrodePitch;
    xM=(1.5:(max_grid_y-0.5))'*electrodePitch;
    yM=(1.5:(max_grid_x-0.5))'*electrodePitch;
    
    %plot electrode grid with numbers
    if isempty(h3D)
        hand.f3D=figure;
        hand.hPlot3D=axes;
    else
        hand.hPlot3D=h3D;
        hand.f3D=hand.hPlot3D.Parent;
    end
    hold on;
    
    if plot3DGrid
        hand.hPlot1=line([XL(1) XL(2)],[xM xM],'LineWidth',1,'Color',gridColor,'Parent',hand.hPlot3D);
        hand.hPlot2=line([yM yM],[YL(1) YL(2)],'LineWidth',1,'Color',gridColor,'Parent',hand.hPlot3D);
    end
    
    if plot3DElectrodeNumbers
        hand.hPlot3=text(Xc,Yc,num2str(ch),'fontsize',elecNumFontSize,'horizontalAlignment','center','Parent',hand.hPlot3D);
    end
    
    if plot3DRemoveNeuronOutSideMEA
        neurons2Plot=find(Xrs>=XL(1) & Xrs<=XL(2) & Yrs>=YL(1) & Yrs<=YL(2));
    else
        neurons2Plot=1:nNeurons;
    end
    
    %plot cell bodies
    [Xsp,Ysp,Zsp] = sphere(6);
    Xsp=Xsp*cellDiameter;Ysp=Ysp*cellDiameter;Zsp=Zsp*cellDiameter;
    cmap=[0.8 0.8 0.8;1 0 0;0 0 1];
    for i=neurons2Plot
        %scatter3(translation(Channels(i),1)+bestY0{i}(j)/100,translation(Channels(i),2)+bestX0{i}(j)/100,bestZ0{i}(j)/100,[],i,'filled');
        surf(Xrs(i)+Xsp,Yrs(i)+Ysp,Zrs(i)+Zsp,'FaceColor', cmap(neuronIdentity(i),:),'EdgeColor','none');
    end
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    box on;
    axis equal;
    camlight('right')
    set(gcf,'color','w');
    xlim([electrodePitch/2 max(Xc)+electrodePitch/2]);
    ylim([electrodePitch/2 max(Yc)+electrodePitch/2]);
    zlim([-cellDiameter maxZ]);
    
    view([-30 20]);
    zlabel('Z [\mum]');
    
    alpha(plot3DTransparency);
    
    if plot3DAddScaleBar
        [hScaleBar]=addScaleBar(hand.hPlot3D,'scaleFac',2,'scaleBarAxes','x','XUnitStr','[\mum]');
    end
    %set(gcf,'PaperPositionMode','auto');print estimatedPositionOfNeurons3DCtr0001 -djpeg -r300;
end

if plotTriangulationOnWaveforms
    
    %plot cell bodies
    cellDiameter=12;
    localGridSize=5;
    
    plotSpheres=false;
    if plotSpheres
        [Xsp,Ysp,Zsp] = sphere(6);
        Xsp=Xsp*cellDiameter;Ysp=Ysp*cellDiameter;Zsp=Zsp*cellDiameter;
    end
    normZ=Zrs/max(Zrs);
    cMap=jet(64);
    
    
    arrayExt=(localGridSize-1)/2;
    [nRowsTmp,nColsTmp]=size(En);
    
    EnExt=NaN(nRowsTmp+arrayExt*2,nColsTmp+arrayExt*2);
    EnExt(1+arrayExt:end-arrayExt,1+arrayExt:end-arrayExt)=En;
    
    if plot1NeuronPerPlotInTriangulation
        nPlotsInRow=1;
        nPlotsInCol=1;
    else
        nPlotsInRow=7;
        nPlotsInCol=4;
    end
    
    c=0;
    for i=1:nNeurons
        nSubAxis=mod(i-1,nPlotsInRow*nPlotsInCol)+1;
        if nSubAxis==1
            c=c+1;
            if plot1NeuronPerPlotInTriangulation
                hand.fSingleNeurons(c)=figure('position',[70 70 750 450]);
            else
                hand.fSingleNeurons(c)=figure('position',[70 70 1500 900]);
            end
        end
        
        h=subaxis(hand.fSingleNeurons(c),nPlotsInCol,nPlotsInRow,nSubAxis,'S',0.01,'M',0.002);hold on;
        
        if plot1NeuronPerPlotInTriangulation
            localGrid=En;
            [yLoc,xLoc]=find(localGrid==pChMax(i));
        else
            [x,y]=find(EnExt==pChMax(i));
            %find the surrounding channels on which feature extraction will be performed
            localGrid=EnExt(x-arrayExt:x+arrayExt,y-arrayExt:y+arrayExt);
            [yLoc,xLoc]=find(localGrid==pChMax(i));
        end
        
        [tmpCh,pInCh]=intersect(ch, localGrid(~isnan(localGrid(:))) );
        
        %plot( (Xc(tmpCh)-Xc(pChMax(i)) )/electrodePitch + arrayExt + 1 , (Yc(tmpCh)-Yc(pChMax(i)) )/electrodePitch + arrayExt + 1 ,'o','color',[0.8 0.8 0.8]);hold on;
        activityTracePhysicalSpacePlot(h,tmpCh,squeeze(WF.avgHPWF(i,pInCh,:)),localGrid,'lockXYRatio',1,'traceColor',[0.9 0.2 0.2],'DrawGrid',1,...
            'gridLineWidth',1,'drawElectrodeCenters',true,'electrodeColor',[0 0.7 0],'electrodeSize',6,'scaling','noOverlap');
        %plot( (Xrs(i)-Xc(pChMax(i)) )/electrodePitch + arrayExt + 1 , (Yrs(i)-Yc(pChMax(i)) )/electrodePitch + arrayExt + 1 ,'o','color',[0.9 0.2 0.2],'lineWidth',3);
        
        if plotSpheres
            surface((Xrs(i)+Xsp-Xc(pChMax(i)) )/electrodePitch + xLoc,...
                (Yrs(i)+Ysp-Yc(pChMax(i)) )/electrodePitch + yLoc,...
                (Zrs(i)+Zsp)/electrodePitch,'FaceColor',cMap(round(normZ(i)*63)+1,:,:),'EdgeColor','none','FaceAlpha',0.4);
            camlight left;
        else
            hold on;plot((Xrs(i)-Xc(pChMax(i)) )/electrodePitch + xLoc,(Yrs(i)-Yc(pChMax(i)) )/electrodePitch + yLoc,'o','MarkerSize',8,'Color',[0.3 0.2 0.8]);
            text(0.5,0.5,['Z=' num2str(round(Zr(i))) ', L=' num2str(round(Lamda(i))) ', V0/V=' num2str(round(V(i)/vMax(i)))]);
        end
        
    end
end


                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYZLV_Model(x,xTmp,yTmp,zTmp,V0tmp,vMax(i),maxPotAtCellBody,maxZ),[0,0,0,V0Grid(j),vMax(i)],optimset('MaxFunEvals',1000000,'Display','off'));
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYZV_Model(x,xTmp,yTmp,zTmp,vMax(i),maxZ),[0,0,0,vMax(i)],optimset('MaxFunEvals',1000000,'Display','off'));
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYZL_Model(x,xTmp,yTmp,zTmp,V0tmp,V0Grid(j),maxZ),[0,0,0,L0],optimset('MaxFunEvals',100000,'Display','off'));
                
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYZV_Model(x,xTmp,yTmp,zTmp,V0tmp,L0,vMax(i),maxPotAtCellBody,maxZ),[0,0,0,vMax(i)],optimset('MaxFunEvals',1000000,'Display','off'));
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYV_Model(x,xTmp,yTmp,V0tmp,LGrid(i),vMax(i),maxZ),[0,0,vMax(i)],optimset('MaxFunEvals',1000000,'Display','off'));
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYZ_Model(x,xTmp,yTmp,zTmp,V0tmp,V0Grid(j),L0,maxZ),[0,0,0],optimset('MaxFunEvals',1000000,'Display','off'));



function score=powerDecay_Model(X,xTmp,yTmp,zTmp,maxZ,V0tmp,V0,L)
score=mean( ( V0tmp- V0./sqrt((xTmp-X(1)).^2+(yTmp-X(2)).^2+(zTmp-X(3)).^2).^L ).^2 );
if X(3)<0 || X(3)>maxZ
    score=score*2;
end

function score=XYZ_Model(X,xTmp,yTmp,zTmp,maxZ,V0tmp,V0,L)
score=mean( ( V0tmp-V0.*exp(-sqrt((xTmp-X(1)).^2+(yTmp-X(2)).^2+(zTmp-X(3)).^2)./L) ).^2 );
if X(3)<0 || X(3)>maxZ
    score=score*2;
end

function score=XYZV_Model(X,xTmp,yTmp,zTmp,maxZ,V0tmp,L,vMax,maxPotAtCellBody)
score=mean( ( V0tmp-X(4).*exp(-sqrt((xTmp-X(1)).^2+(yTmp-X(2)).^2+(zTmp-X(3)).^2)./L) ).^2 );
if X(3)<=0 || X(3)>maxZ || X(4)<=vMax || X(4)>=maxPotAtCellBody
    score=score*2;
end

function score=XYZL_Model(X,xTmp,yTmp,zTmp,maxZ,V0tmp,V0)
score=mean( ( V0tmp-V0.*exp(-sqrt((xTmp-X(1)).^2+(yTmp-X(2)).^2+(zTmp-X(3)).^2)./X(4)) ).^2 );
if X(3)<0 || X(3)>maxZ
    score=score*2;
end

function score=XYZLV_Model(X,xTmp,yTmp,zTmp,maxZ,V0tmp,vMax,maxPotAtCellBody)
score=mean( (  V0tmp-X(5).*exp(-sqrt((xTmp-X(1)).^2+(yTmp-X(2)).^2+(zTmp-X(3)).^2)./X(4))         ).^2);
if X(3)<0 || X(3)>maxZ || X(5)<vMax || X(5)>maxPotAtCellBody
    score=score*2;
end

function score=LogXYZV_Model(X,xTmp,yTmp,zTmp,maxZ,V0tmp,L,vMax)
score=mean((log(V0tmp./X(4))  +  sqrt((xTmp-X(1)).^2+(yTmp-X(2)).^2+(zTmp-X(3)).^2)./L).^2);
if X(3)<0 || X(3)>maxZ || X(4)<vMax
    score=score*2;
end

function score=XYV_Model(X,xTmp,yTmp,V0tmp,L,vMax)
score=sum( ( V0tmp-X(3).*exp(-sqrt((xTmp-X(1)).^2+(yTmp-X(2)).^2)./L) ).^2 );
if X(3)<vMax || X(3)>maxZ
    score=score*2;
end

function score=potEst(X,Xe,Ye,Ve)
%[0,0,0,LGrid(j),mxV0]
%x=[X Y Z I]
%P0=@(Xe,Ye) 2.*Fih(Xe-X,Ye-Y,-Z) +2.*W_TS.*sum(    Fih(Xe-X,Ye-Y,-Z+2.*(1:nTerms).*h)  +  Fih(Xe-X,Ye-Y,-Z-2.*(1:nTerms).*h)   );
nTerms=20;
W_TS=-0.6667;
h=400;
sigmaT=0.3;

Fih=@(u,v,w) X(4)./(   (4*pi*sigmaT) .* sqrt(u.^2+v.^2+w.^2)   );

for n=1:nTerms
    tmp(n,:)=(W_TS.^n) .* (  Fih(Xe-X(1),Ye-X(2),-X(3)+2.*n.*h)  +  Fih(Xe-X(1),Ye-X(2),-X(3)-2.*n.*h));
end

P0=2.*Fih(Xe-X(1),Ye-X(2),-X(3)) +2.*sum(tmp);

score=sum((P0-Ve).^2);

if X(3)<0
    score=score*2;
end