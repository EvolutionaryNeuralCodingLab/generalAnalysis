function [est,hand,testData]=spikePositionEstimation(WF,ch,preMs,Fs,En,varargin)
% [est,h,testData]=spikePositionEstimation(WF,ch,preMs,Fs,En,varargin)
% Function purpose : Calculate parameters of spike shapes
%
% Function recives :    WF - average spike waveforms over all electrodes in ch [Double [neurons x ch x samples]
%                       ch - the channel numbers of the channels in avgWFWaveform [NChannels,Time] - the raw voltage samples of all channels
%                       preMs - pre spike time in WF [ms]
%                       Fs - sampling frequency of the WFs
%                       En - electrode layout
%                       neuronIdentity [1 x neurons] - neuron identity (1=none, 2=Inhibitory, 3=Excitatory)]
%
%                       varargin ('property name','property value')
%
% Function give back :  par - a structure of output parameters
%                       hand - a structure of handles from generated plots
%
% Last updated : 29/07/15

%default parameters
smoothingDuration=0.5; %[ms]
maxDistForElectrodes=190; %[uM]
electrodePitch=100; %[uM]
useSameTimePoint=true; %if true takes the time point of max channel of every neuron and uses this time for all channels, if false checks for the max of each ch seperately
dXY=1; %[uM]
dZ=1; %[uM]
dL=1;
dV0=1; %[uV]
maxZ=200; %[uM]
minMaxV0=[100 300]; %[uV]
minMaxL=[50 150]; %[uM]

gridColor=[0.8 0.8 0.8];
cellDiameter=15; %[um]
modelType='reducedModelOpt';%,'realModelOpt','globalMinima','reducedModelOpt';
smoothSpikes=false;
L0=70;
V0=160;
Xc=[];
Yc=[];
useSpikeExtremum=true;
est=[]; %if this is given as input, no calculation is made but plots are still calculated
LGrid=25:1:125; %decay constants to use as initial conditions in optimization
V0Grid=100:5:300;
Icell=10.^(-10:0.05:-7); %current density in initial axon segment
neuronIdentity=[];

plot3D=true;
elecNumFontSize=5;
plotTriangulationOnWaveforms=true;
plot1NeuronPerPlotInTriangulation=false;

%% Output list of default variables
%print out default arguments and values if no inputs are given
if nargin==0
    defaultArguments=who;
    for i=1:numel(defaultArguments)
        eval(['defaultArgumentValue=' defaultArguments{i} ';']);
        if numel(defaultArgumentValue)==1
            disp([defaultArguments{i} ' = ' num2str(defaultArgumentValue)]);
        else
            fprintf([defaultArguments{i} ' = ']);
            disp(defaultArgumentValue);
        end
    end
    return;
end

%Collects all options
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};']);
end

hand=[];testData=[];

if isempty(WF) %run on artificial data
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
    
    %Build inverse map between electrode and location
    [meshX,meshY]=meshgrid(1:size(En,1),1:size(En,2));
    Xc(En(~isnan(En)))=meshX(~isnan(En))*electrodePitch;
    Yc(En(~isnan(En)))=meshY(~isnan(En))*electrodePitch;
    
    smoothingSamples=round(smoothingDuration*Fs/1000);
    preSpkSamples=preMs*Fs/1000;
    WF=zeros([nNeurons,nCh,nSpikeSamples]);
    for i=1:nNeurons
        for j=1:nCh
            WF(i,j,preSpkSamples:(preSpkSamples+3*smoothingSamples))=...
                -tV0(i).*exp(-sqrt(  (Xc(tCh(i))+tX(i)-Xc(j)).^2 + (Yc(tCh(i))+tY(i)-Yc(j)).^2+(tZ(i)-0).^2   )./tL(i));
            %WF(i,j,preSpkSamples:(preSpkSamples+3*smoothingSamples))=...
            %    -tV0(i).*exp(-sqrt(  (Xc(tCh(i))+tX(i)-Xc(j)).^2 + (Yc(tCh(i))+tY(i)-Yc(j)).^2   )./tL(i));
        end
        %h=axes;activityTracePhysicalSpacePlot(h,1:nCh,squeeze(WF(i,:,:)),En,'traceColor',[0.2 0.2 0.8],'scaling','none');
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
[nNeurons,nCh,nSpikeSamples]=size(WF);

if isempty(neuronIdentity) %if no an empty neuron identity is given
    neuronIdentity=ones(1,nNeurons);
else
    neuronIdentity(neuronIdentity==0)=1;
end

% calculate pre and post samples and their position in the waveform
preSpkSamples=preMs*Fs/1000;
postSpkSamples=nSpikeSamples-preSpkSamples;


%% model with x,y,z,V
%according to the current grid arrangement positive Y is down on real MEA space and positive X is right on real MEA space

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

%Build inverse map between electrode and location
[meshX,meshY]=meshgrid(1:size(En,1),1:size(En,2));
if isempty(Xc)
    Xc(En(~isnan(En)))=meshX(~isnan(En))*electrodePitch;
    Xc=Xc(ch);
end
if isempty(Yc)
    Yc(En(~isnan(En)))=meshY(~isnan(En))*electrodePitch;
    Yc=Yc(ch);
end

% find spike extremum ch and position
if useSpikeExtremum
    [vMaxAll,pMaxAll]=max(abs(WF),[],3);
    [vMax,pChMax]=max(vMaxAll,[],2);
else
    [vMaxAll,pMaxAll]=min(WF,[],3);
    [vMax,pChMax]=min(vMaxAll,[],2);
    vMax=abs(vMax);
    vMaxAll=abs(vMaxAll);
end
pMaxSampleInMaxCh=pMaxAll( sub2ind(size(pMaxAll), 1:nNeurons, pChMax') );

if isempty(est) %if no previous estimation was entered
    
    %smooth spike shapes with local linear regression
    if smoothSpikes
        smoothingSamples=round(smoothingDuration*Fs/1000);
        for i=1:nNeurons
            for j=1:nCh
                WF(i,j,:) = smooth(WF(i,j,:),smoothingSamples,'loess');
            end
        end
    end
    
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
    
    hWait=waitbar(0,'Triangulating neurons...');
    for i=1:nNeurons
        waitbar(i/nNeurons,hWait);
        
        pXY=sqrt((Xc-Xc(pChMax(i))).^2+(Yc-Yc(pChMax(i))).^2)<maxDistForElectrodes;
        xTmp=Xc(pXY)-Xc(pChMax(i));
        yTmp=Yc(pXY)-Yc(pChMax(i));
        zTmp=zeros(size(yTmp));
        
        if useSameTimePoint
            V0tmp=abs(WF(i,pXY,pMaxSampleInMaxCh(i)));
        else
            V0tmp=vMaxAll(i,pXY);
        end
        
        if strcmp(modelType,'globalMinima')
            
            score=zeros(size(XAll));
            for k=1:numel(xTmp)
                score=score+(       V0tmp(k)-V0All.*exp(-sqrt((xTmp(k)-XAll).^2+(yTmp(k)-YAll).^2+(zTmp(k)-ZAll).^2)./LAll)         ).^2;
                %score=score+(log(V0tmp(k)./V0All)  +  sqrt((xTmp(k)-XAll).^2+(yTmp(k)-YAll).^2+(zTmp(k)-ZAll).^2)./L0).^2;
            end
            %figure;imagesc(Z,squeeze(V0),squeeze(score(find(X==0),find(X==0),:,:)));xlabel('Z');ylabel('V_0');colorbar
            
            [bestVal(i),bestPlace]=min(score(:));
            [bestXp,bestYp,bestZp,bestLp,bestV0p]=ind2sub(size(V0All),bestPlace);
            est.X(i)=X(bestXp);
            est.Y(i)=Y(bestYp);
            est.Z(i)=Z(bestZp);
            est.L(i)=L(bestLp);
            est.V0(i)=V0(bestV0p);
            
        elseif strcmp(modelType,'realModelOpt')

            eVal=Inf;
            for j=1:numel(Icell)
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYZLV_Model(x,xTmp,yTmp,zTmp,V0tmp,mxV0),[0,0,0,LGrid(j),mxV0],optimset('MaxFunEvals',1000000,'Display','off'));
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYZL_Model(x,xTmp,yTmp,zTmp,V0tmp,V0),[0,0,0,LGrid(j)],optimset('MaxFunEvals',100000));
                
                [estTmpTmp,eValTmp] = fminsearch(@(x) potEst(x,xTmp*1e-6,yTmp*1e-6,V0tmp*1e-6),[0,0,0,Icell(j)],optimset('MaxFunEvals',1000000,'Display','off'));
                allEst(j)=eValTmp;
                if eValTmp<eVal
                    estTmp=estTmpTmp;
                    eVal=eValTmp;
                end
            end
            %plot(Icell,allEst,'.');set(gca,'XScale','log');pause;
            
            est.X(i)=estTmp(1);
            est.Y(i)=estTmp(2);
            est.Z(i)=estTmp(3);
            est.eVal(i)=eVal;
            est.I(i)=estTmp(4);

        elseif strcmp(modelType,'reducedModelOpt')
            
            eVal=Inf;
            mxV0=max(V0tmp);
            for j=1:numel(LGrid)
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYZLV_Model(x,xTmp,yTmp,zTmp,V0tmp,mxV0),[0,0,0,LGrid(j),mxV0],optimset('MaxFunEvals',1000000,'Display','off'));
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYZV_Model(x,xTmp,yTmp,zTmp,mxV0),[0,0,0,mxV0],optimset('MaxFunEvals',1000000,'Display','off'));
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYZL_Model(x,xTmp,yTmp,zTmp,V0tmp,V0),[0,0,0,LGrid(j)],optimset('MaxFunEvals',100000));
                
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYZV_Model(x,xTmp,yTmp,zTmp,V0tmp,L0,vMax(i)),[0,0,0,vMax(i)],optimset('MaxFunEvals',1000000));
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYV_Model(x,xTmp,yTmp,V0tmp,LGrid(i),vMax(i)),[0,0,vMax(i)],optimset('MaxFunEvals',1000000));
                [estTmpTmp,eValTmp] = fminsearch(@(x) XYZ_Model(x,xTmp,yTmp,zTmp,V0tmp,V0,LGrid(j)),[0,0,0],optimset('MaxFunEvals',1000000,'Display','off'));
                
                if eValTmp<eVal
                    estTmp=estTmpTmp;
                    eVal=eValTmp;
                    tmpL=LGrid(j);
                end
            end
            
            est.X(i)=estTmp(1);
            est.Y(i)=estTmp(2);
            est.Z(i)=estTmp(3);
            est.eVal(i)=eVal;
            est.L(i)=tmpL;%est.L(i)=estTmp(4);
            est.V0(i)=V0; %est.V0(i)=estTmp(5);
            
        elseif strcmp(modelType,'reducedModelOptV0')
            
            eVal=Inf;
            mxV0=max(V0tmp);
            for j=1:numel(V0Grid)
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYZLV_Model(x,xTmp,yTmp,zTmp,V0tmp,mxV0),[0,0,0,LGrid(j),mxV0],optimset('MaxFunEvals',1000000,'Display','off'));
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYZV_Model(x,xTmp,yTmp,zTmp,mxV0),[0,0,0,mxV0],optimset('MaxFunEvals',1000000,'Display','off'));
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYZL_Model(x,xTmp,yTmp,zTmp,V0tmp,V0),[0,0,0,LGrid(j)],optimset('MaxFunEvals',100000));
                
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYZV_Model(x,xTmp,yTmp,zTmp,V0tmp,L0,vMax(i)),[0,0,0,vMax(i)],optimset('MaxFunEvals',1000000));
                %[estTmpTmp,eValTmp] = fminsearch(@(x) XYV_Model(x,xTmp,yTmp,V0tmp,LGrid(i),vMax(i)),[0,0,vMax(i)],optimset('MaxFunEvals',1000000));
                [estTmpTmp,eValTmp] = fminsearch(@(x) XYZ_Model(x,xTmp,yTmp,zTmp,V0tmp,V0Grid(j),L0),[0,0,0],optimset('MaxFunEvals',1000000,'Display','off'));
                
                if eValTmp<eVal
                    estTmp=estTmpTmp;
                    eVal=eValTmp;
                    tmpV0=V0Grid(j);
                end
            end
            
            est.X(i)=estTmp(1);
            est.Y(i)=estTmp(2);
            est.Z(i)=estTmp(3);
            est.eVal(i)=eVal;
            est.L(i)=L0;%est.L(i)=estTmp(4);
            est.V0(i)=tmpV0; %est.V0(i)=estTmp(5);
            
        end
        est.Xrs(i)=est.X(i)+Xc(pChMax(i));
        est.Yrs(i)=est.Y(i)+Yc(pChMax(i));
        est.Zrs(i)=est.Z(i);
        est.electrodeX=Xc;
        est.electrodeY=Yc;
    end
    close(hWait);
end
%{
figure;subplot(2,3,1);plot(tL-est.L,'.');title('L');subplot(2,3,2);plot(tX-est.X,'.');title('X');subplot(2,3,3);plot(tY-est.Y,'.');title('Y');subplot(2,3,4);plot(tZ-est.Z,'.');title('Z');subplot(2,3,5);plot(tV0-est.V0,'.');title('V0');
%}
if plot3D
    
    %create electrode grid with numbers for plotting
    [max_grid_y,max_grid_x]=size(En);
    XL=[0.5 max_grid_x+0.5]*electrodePitch;
    YL=[0.5 max_grid_y+0.5]*electrodePitch;
    xM=(1.5:(max_grid_y-0.5))'*electrodePitch;
    yM=(1.5:(max_grid_x-0.5))'*electrodePitch;
    
    %plot electrode grid with numbers
    hand.f3D=figure;
    hand.hPlot3D=axes;
    hPlot1=line([XL(1) XL(2)],[xM xM],'LineWidth',1,'Color',gridColor,'Parent',hand.hPlot3D);
    hPlot2=line([yM yM],[YL(1) YL(2)],'LineWidth',1,'Color',gridColor,'Parent',hand.hPlot3D);
    hPlot3=text(Xc,Yc,num2str(ch),'fontsize',elecNumFontSize,'horizontalAlignment','center','Parent',hand.hPlot3D);
    hold on;
    
    %plot cell bodies
    [Xsp,Ysp,Zsp] = sphere(6);
    Xsp=Xsp*cellDiameter;Ysp=Ysp*cellDiameter;Zsp=Zsp*cellDiameter;
    cmap=[0.2 0.4 0.6;0.8 0.2 0.2;0.2 0.2 0.8];
    for i=1:nNeurons
        %scatter3(translation(Channels(i),1)+bestY0{i}(j)/100,translation(Channels(i),2)+bestX0{i}(j)/100,bestZ0{i}(j)/100,[],i,'filled');
        surf(est.Xrs(i)+Xsp,est.Yrs(i)+Ysp,est.Zrs(i)+Zsp,'FaceColor', cmap(neuronIdentity(i),:),'EdgeColor','none');
    end
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    box on;
    axis equal;
    camlight left;
    set(gcf,'color','w');
    xlim([electrodePitch/2 max(Xc)+electrodePitch/2]);
    ylim([electrodePitch/2 max(Yc)+electrodePitch/2]);
    zlim([0 maxZ]);
    
    view([-30 20]);
    zlabel('Z [\mum]');
    %set(gcf,'PaperPositionMode','auto');print estimatedPositionOfNeurons3DCtr0001 -djpeg -r300;
end
    
if plotTriangulationOnWaveforms
    
    %plot cell bodies
    cellDiameter=12;
    plotSpheres=false;
    if plotSpheres
        [Xsp,Ysp,Zsp] = sphere(6);
        Xsp=Xsp*cellDiameter;Ysp=Ysp*cellDiameter;Zsp=Zsp*cellDiameter;
    end
    normZ=est.Zrs/max(est.Zrs);
    cMap=jet(64);
    
    localGridSize=5;
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
        activityTracePhysicalSpacePlot(h,tmpCh,squeeze(WF(i,pInCh,:)),localGrid,'lockXYRatio',1,'traceColor',[0.9 0.2 0.2],'DrawGrid',1,...
        'gridLineWidth',1,'drawElectrodeCenters',true,'electrodeColor',[0 0.7 0],'electrodeSize',6,'scaling','noOverlap');
        %plot( (est.Xrs(i)-Xc(pChMax(i)) )/electrodePitch + arrayExt + 1 , (est.Yrs(i)-Yc(pChMax(i)) )/electrodePitch + arrayExt + 1 ,'o','color',[0.9 0.2 0.2],'lineWidth',3);
        
        if plotSpheres
            surface((est.Xrs(i)+Xsp-Xc(pChMax(i)) )/electrodePitch + xLoc,...
                (est.Yrs(i)+Ysp-Yc(pChMax(i)) )/electrodePitch + yLoc,...
                (est.Zrs(i)+Zsp)/electrodePitch,'FaceColor',cMap(round(normZ(i)*63)+1,:,:),'EdgeColor','none','FaceAlpha',0.4);
            camlight left;
        else
            hold on;plot((est.Xrs(i)-Xc(pChMax(i)) )/electrodePitch + xLoc,(est.Yrs(i)-Yc(pChMax(i)) )/electrodePitch + yLoc,'o','MarkerSize',8,'Color',[0.3 0.2 0.8]);
            text(0.5,0.5,['Z=' num2str(round(est.Z(i))) ', L=' num2str(round(est.L(i))) ', V0/V=' num2str(round(est.V0(i)/vMax(i)))]);
        end

    end
end

function score=XYZV_Model(X,xTmp,yTmp,zTmp,V0tmp,L,vMax)
score=sum( ( V0tmp-X(4).*exp(-sqrt((xTmp-X(1)).^2+(yTmp-X(2)).^2+(zTmp-X(3)).^2)./L) ).^2 );
if X(3)<0 || X(4)<vMax
    score=score*2;
end

function score=XYV_Model(X,xTmp,yTmp,V0tmp,L,vMax)
score=sum( ( V0tmp-X(3).*exp(-sqrt((xTmp-X(1)).^2+(yTmp-X(2)).^2)./L) ).^2 );
if X(3)<vMax
    score=score*2;
end

function score=XYZ_Model(X,xTmp,yTmp,zTmp,V0tmp,V0,L)
score=sum( ( V0tmp-V0.*exp(-sqrt((xTmp-X(1)).^2+(yTmp-X(2)).^2+(zTmp-X(3)).^2)./L) ).^2 );
if X(3)<0
    score=score*2;
end

function score=LogXYZV_Model(X,xTmp,yTmp,zTmp,V0tmp,L,vMax)
score=sum((log(V0tmp./X(4))  +  sqrt((xTmp-X(1)).^2+(yTmp-X(2)).^2+(zTmp-X(3)).^2)./L).^2);
if X(3)<0 || X(4)<vMax
    score=score*2;
end

function score=XYZL_Model(X,xTmp,yTmp,zTmp,V0tmp,V0)
score=sum( ( V0tmp-V0.*exp(-sqrt((xTmp-X(1)).^2+(yTmp-X(2)).^2+(zTmp-X(3)).^2)./X(4)) ).^2 );
if X(3)<0
    score=score*2;
end

function score=XYZVL_Model(X,xTmp,yTmp,zTmp,V0tmp)
score=sum( (  V0tmp-X(4).*exp(-sqrt((xTmp-X(1)).^2+(yTmp-X(2)).^2+(zTmp-X(3)).^2)./X(5))         ).^2);
if X(3)<0
    score=score*2;
end

function score=XYZLV_Model(X,xTmp,yTmp,zTmp,V0tmp,vMax)
score=sum( (  V0tmp-X(5).*exp(-sqrt((xTmp-X(1)).^2+(yTmp-X(2)).^2+(zTmp-X(3)).^2)./X(4))         ).^2);
if X(3)<0 || X(5)<vMax
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