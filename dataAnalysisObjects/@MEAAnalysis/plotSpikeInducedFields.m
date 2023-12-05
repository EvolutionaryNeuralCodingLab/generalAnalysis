function [hand,out]=plotSpikeInducedFields(obj,varargin)

obj.checkFileRecording;

parseObj = inputParser;

addParameter(parseObj,'outArg',{});

addParameter(parseObj,'pNeu',[],@isnumeric); %the positions of the neurons to plot
addParameter(parseObj,'angle2LateralRight',0,@isnumeric); %The angle in degrees to rotate the axes so that lateral is on the left
addParameter(parseObj,'reverseDirction',0,@isnumeric); %whether flip direction is needed to make anterior to top when lateral is left
addParameter(parseObj,'spikeEdgeBand',0,@isnumeric); %the band arround the electrode array edge (from farthest electrode) to include in stats of neuron positions (if positive, takes also neurons outside the array)
addParameter(parseObj,'fieldEdgeBand',-50,@isnumeric); %the band arround the electrode array edge (from farthest electrode) to include in stats of field positions (if positive, takes also neurons outside the array)

addParameter(parseObj,'hand',[]);
addParameter(parseObj,'dAngle4Plot',30,@isnumeric);
addParameter(parseObj,'maxFields4Plot',375,@isnumeric);
addParameter(parseObj,'plotElectrodeNames',true,@isnumeric);
addParameter(parseObj,'plotFieldVectors',true,@isnumeric);
addParameter(parseObj,'polarPlot',true,@isnumeric);
addParameter(parseObj,'polarPlotPlotNeuronNumber',false,@isnumeric);
addParameter(parseObj,'plotNeuronNumbersAllFields',false,@isnumeric);
addParameter(parseObj,'polarPlotRemoveOuliers',false,@isnumeric);
addParameter(parseObj,'plotFieldMapAllNeurons',false,@isnumeric);
addParameter(parseObj,'polarPlotDistanceRange',[0 500],@isnumeric);
addParameter(parseObj,'plotFieldVectorsOrientation',true,@isnumeric);
addParameter(parseObj,'deleteAngleText',true,@isnumeric);
addParameter(parseObj,'EC',[0 0.2 0.8],@isnumeric);
addParameter(parseObj,'IC',[0.95 0.1 0.05],@isnumeric);

addParameter(parseObj,'normalizeColorCode',true,@isnumeric);

addParameter(parseObj,'extrapolateMaxima',true,@isnumeric);
addParameter(parseObj,'markerSizeAllFields',15,@isnumeric);

addParameter(parseObj,'fileNameSpikeInducedFields',[],@isstr);
addParameter(parseObj,'overwrite',false,@isnumeric);
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

[funName] = dbstack;funName=funName.name;funName=strsplit(funName,'.');funName=funName{end};
saveFileName=obj.files.(funName);

obj.checkFileRecording(obj.files.getSpikePositionEstimation,'Spike position estimation file missing, please run getSpikePositionEstimation');

if isempty(fileNameSpikeInducedFields)
    obj.checkFileRecording(obj.files.getSpikeInducedFields,'Spike induced fields file missing, please run getSpikeInducedFields');
    load(obj.files.getSpikeInducedFields,'fieldPar','pMaxField');
else
    load(fileNameSpikeInducedFields,'fieldPar','pMaxField');
end

%populate grid sorter object
%obj=populateGridSorterObj(obj);

%assign parameters
ch=obj.currentDataObj.channelNumbers;
En=obj.currentDataObj.chLayoutNumbers;
%Fs=obj.currentDataObj.samplingFrequency(1);
%preSpikeMs=obj.gridSorterObj.postPreRawWindow;


[Xc,Yc]=obj.currentDataObj.getElectrodePositions;
%load(obj.files.getSpikePositionEstimation,'Xrs','Yrs','Zrs');
Xrs=fieldPar.X(1,:);Yrs=fieldPar.Y(1,:);


nNeu=numel(fieldPar.mag);

%% Plotting results

rotAngRad=angle2LateralRight/180*pi;
fieldPar.anglePolar=fieldPar.angle+rotAngRad;
if reverseDirction
    fieldPar.anglePolar=-fieldPar.anglePolar;
end

%modify this to include
pExcit=find(fieldPar.classIE==3); %excitatory
pInhib=find(fieldPar.classIE==2);
    
%select sub population
if isempty(pNeu)
    pNeuron=true(1,nNeu);
else
    pNeuron=false(1,nNeu);
    pNeuron(pNeu)=true;
end

%calculate points outside electrode area
[pBound] = boundary(Xc',Yc'); %Calculate bounding points
mX=mean([Xc;Yc],2);
[teta,r]=cart2pol((Xc(pBound)-mX(1)),Yc(pBound)-mX(2));

rSpike=r+spikeEdgeBand;
[xB,yB]=pol2cart(teta,rSpike);
inSpike = inpolygon(Xrs,Yrs,xB+mX(1),yB+mX(2));
%plot(xB+mX(1),yB+mX(2));hold on;plot(Xrs,Yrs,'.');plot(Xrs(inSpike),Yrs(inSpike),'o');plot(Xc(pBound),Yc(pBound));

rSIF=r+fieldEdgeBand;
[xB,yB]=pol2cart(teta,rSIF);
inSIF = inpolygon(fieldPar.Xfield,fieldPar.Yfield,xB+mX(1),yB+mX(2));
%plot(xB+mX(1),yB+mX(2));hold on;plot(fieldPar.Xfield,fieldPar.Yfield,'.');plot(fieldPar.Xfield(inSIF),fieldPar.Yfield(inSIF),'o');plot(Xc(pBound),Yc(pBound));

fieldPar.edgeNeurons=false(1,nNeu);%initialization
fieldPar.edgeNeurons(~inSIF | ~inSpike)=true;

if ~isempty(polarPlotDistanceRange)
    pDist=fieldPar.mag>=polarPlotDistanceRange(1) & fieldPar.mag<=polarPlotDistanceRange(2);
else
    pDist=true(1,nNeu);
end
pI=intersect(find(pDist & ~fieldPar.edgeNeurons & pNeuron),pInhib);
pE=intersect(find(pDist & ~fieldPar.edgeNeurons & pNeuron),pExcit);

if polarPlot
    %prepare for plotting
    f=figure('position',[100 100 500 500]);
    P = panel(f);
    P.pack(2,2);
    P.margin=8;
    
    angleBins=(dAngle4Plot/360/2*pi):(dAngle4Plot/360*pi):(pi*2);
    if polarPlotRemoveOuliers
        maximalMag=median(fieldPar.mag([pI pE]))+6*mad(fieldPar.mag([pI pE]),1);
    else
        maximalMag=max(fieldPar.mag([pI pE]));
    end
    %inhibitory
    hand.polarAx(1,1)=P(1, 1).select();
    hRose=rose(fieldPar.anglePolar(pI),angleBins);
    hRose.Color=IC;
    XdataRose = get(hRose,'Xdata');XdataRose=reshape(XdataRose,[4,numel(XdataRose)/4]);
    YdataRose = get(hRose,'Ydata');YdataRose=reshape(YdataRose,[4,numel(YdataRose)/4]);
    hPatch=patch(XdataRose,YdataRose,IC);    
    set(gca,'color','k');
    xlimRose=xlim;
    
    %calculate mean angle
    mAngleI=angle(mean(exp(1i*fieldPar.anglePolar(pI))));
    %calculate mean magnitude by projecting all vectors on the unity vector in the direction of the average angle
    [mU,mV]=pol2cart(mAngleI,1);
    [U,V]=pol2cart(fieldPar.anglePolar(pI),fieldPar.mag(pI));
    mMagI=mean([U' V']*[mU;mV]);
    
    %[UI,VI]=pol2cart(mAngleI,mMagI);
    hold on;
    [UI,VI]=pol2cart(mAngleI,xlimRose(2));
    hand.hCompI(1)=compass(UI,VI,'r');
    hand.hCompI(1).Color=IC;
    
    hand.polarAx(1,2)=P(1, 2).select();
    polar(0,maximalMag,'-k');hold on; %set scale for polar plot
    hand.PolarI=polar(fieldPar.anglePolar(pI),fieldPar.mag(pI),'.r');
    hand.PolarI.Color=IC;
    
    xlimPolar=xlim;
    [UI,VI]=pol2cart(mAngleI,xlimPolar(2));
    hand.hCompI(2)=compass(UI,VI,'r');
    hand.hCompI(2).Color=IC;
    
    %excitatory
    hand.polarAx(2,1)=P(2, 1).select();
    hRose=rose(fieldPar.anglePolar(pE),angleBins);
    XdataRose = get(hRose,'Xdata');XdataRose=reshape(XdataRose,[4,numel(XdataRose)/4]);
    YdataRose = get(hRose,'Ydata');YdataRose=reshape(YdataRose,[4,numel(YdataRose)/4]);
    hPatch=patch(XdataRose,YdataRose,EC);
    set(gca,'color','k');
    
    %calculate mean angle
    mAngleE=angle(mean(exp(1i*fieldPar.anglePolar(pE))));
    %calculate mean magnitude by projecting all vectors on the unity vector in the direction of the average angle
    [mU,mV]=pol2cart(mAngleE,1);
    [U,V]=pol2cart(fieldPar.anglePolar(pE),fieldPar.mag(pE));
    mMagE=mean([U' V']*[mU;mV]);
    
    hold on;
    xlimRose=xlim;
    [UE,VE]=pol2cart(mAngleE,xlimRose(2));
    hand.hCompE(1)=compass(UE,VE,'b');
    hand.hCompE(1).Color=EC;
    
    hand.polarAx(2,2)=P(2, 2).select();
    polar(0,maximalMag,'-k');hold on; %set scale for polar plot
    hand.PolarE=polar(fieldPar.anglePolar(pE),fieldPar.mag(pE),'.b');
    hand.PolarE.Color=EC;

    
    [UE,VE]=pol2cart(mAngleE,xlimPolar(2));
    hand.hCompE(2)=compass(UE,VE,'b');
    hand.hCompE(2).Color=EC;
    
    %delete angle text
    if deleteAngleText
        textObj1=findall(hand.polarAx(1,1),'type','text');
        textObj2=findall(hand.polarAx(1,2),'type','text');
        textObj3=findall(hand.polarAx(2,1),'type','text');
        textObj4=findall(hand.polarAx(2,2),'type','text');
        
        %txt=get(textObj,'string');
        delete([textObj1(1:12);textObj2(1:12);textObj3(1:12);textObj4(1:12)]);
    end
    
    if polarPlotPlotNeuronNumber
        text(hand.polarAx(2,2),hand.PolarE.XData',hand.PolarE.YData',num2str(pE'),'FontSize',8);
        text(hand.polarAx(1,2),hand.PolarI.XData',hand.PolarI.YData',num2str(pI'),'FontSize',8);
    end
end

%DSI=(prefered - (prefered+pi))/(prefered + (prefered+pi))
if plotFieldVectors
    f=figure('position',[100 100 700 700]);
    hand.hVec=axes;
    hand.hVec.WarpToFill='off'; %to avoid error in arrow3 function
    
    if plotElectrodeNames
        hand.electrodeText=text(Xc,Yc,num2str(ch'),'fontsize',8,'Parent',hand.hVec,'horizontalAlignment','center');
        xlim([min(Xc)-obj.currentDataObj.electrodePitch max(Xc)+obj.currentDataObj.electrodePitch]);
        ylim([min(Yc)-obj.currentDataObj.electrodePitch max(Yc)+obj.currentDataObj.electrodePitch]);
        hold(hand.hVec,'on');
    end
    
    %hQ=quiver(Xc(neuronNames(1,:)),Yc(neuronNames(1,:)),intdX,intdY,'filled','lineWidth',2,'MaxHeadSize',0.1,'color','k','MarkerSize',2,'MarkerFaceColor','k');
    [tmpX,tmpY]=pol2cart(fieldPar.angle,50);
    
    nInhib2Display=numel(pI);
    cMapR=flipud([ones(1,60);(0:0.01:0.59);(0:0.01:0.59)]');
    normColorI = floor(min(fieldPar.mag(pI)./maximalMag,1).*(size(cMapR,1)-1))+1;
    if ~isempty(pI)
        hand.hArrowI=arrow3([Xrs(pI);Yrs(pI)]',[Xrs(pI)+tmpX(pI);Yrs(pI)+tmpY(pI)]','^r2',0.7,1);hold on;
        for i=1:nInhib2Display
            hand.hArrowI(i+1).FaceColor=cMapR(normColorI(i),:,:);
        end
    end
    nExcit2Display=numel(pE);
    cMapB=flipud([(0:0.01:0.59);(0:0.01:0.59);ones(1,60)]');
    normColorE = floor(min(fieldPar.mag(pE)./maximalMag,1).*(size(cMapR,1)-1))+1;
    if ~isempty(pE)
        hand.hArrowE=arrow3([Xrs(pE);Yrs(1,pE)]',[Xrs(pE)+tmpX(pE);Yrs(1,pE)+tmpY(pE)]','^b2',0.7,1);
        for i=1:nExcit2Display
            hand.hArrowE(i+1).FaceColor=cMapB(normColorE(i)     ,:,:);
        end
    end
    xlabel('X [\mum]','FontSize',14);
    ylabel('Y [\mum]','FontSize',14);
    
    if plotFieldVectorsOrientation
        scale=100;shift=150;
        rotMat=[cos(rotAngRad) -sin(rotAngRad);sin(rotAngRad) cos(rotAngRad)];
        x1=[0 -1;-1 0;0 1;1 0];
        x2=[0 1;1 0;0 -1;-1 0];
        xT(:,1)=[0;-1.2;0;1.2];
        xT(:,2)=[-1.2;0;1.2;0];
        x1=x1*rotMat*scale+shift;
        x2=x2*rotMat*scale+shift;
        xT=xT*rotMat*scale+shift;

        hand.hOrientation=arrow3(x1,x2,'^k2',0.7,1);
        if reverseDirction
            text(xT(:,1),xT(:,2),{'A','M','P','L'},'horizontalAlignment','center','verticalAlignment','middle')
        else
            text(xT(:,1),xT(:,2),{'P','M','A','L'},'horizontalAlignment','center','verticalAlignment','middle')
        end
        
    end
end

if plotFieldMapAllNeurons
    nNeurons=numel(fieldPar.mag);
    if normalizeColorCode
        Ilim=0;
    else
        Ilim=[min(fieldPar.val(:)) max(fieldPar.val(:))]; 
    end
    n=ceil(sqrt(min(maxFields4Plot,nNeurons)/3/5));%define images in a 3 x 5 ratio
    xPlots=n*5;
    yPlots=n*3;
    nPlotPerPage=xPlots*yPlots;
    
    f=figure;
    P = panel(f);
    P.pack(yPlots,xPlots);
    P.margin=0.001;
    
    for i=1:nNeurons
        hand.hAllFieldAxes(i)=P(ceil(i/xPlots),i-(ceil(i/xPlots)-1)*xPlots).select();
        IntensityPhysicalSpacePlot(1:120,fieldPar.val(i,:),obj.currentDataObj.chLayoutNumbers,'plotElectrodeNumbers',0,'plotGridLines',0,'plotColorBar',0,'markerSize',markerSizeAllFields,'h',hand.hAllFieldAxes(i),'Ilim',Ilim);
        
        text(Xc(pChMax(i))/electrodePitch-0.5,Yc(pChMax(i))/obj.currentDataObj.electrodePitch-0.5,'o','horizontalAlignment','center','fontsize',6);
        if plotNeuronNumbersAllFields
            text(0,0,num2str(i),'horizontalAlignment','left','verticalAlignment','bottom','fontsize',6);
        end
        line( [Xc(pChMax(i)) fieldPar.Xfield(i)]/obj.currentDataObj.electrodePitch - 0.5 , [Yc(pChMax(i)) fieldPar.Yfield(i)]/obj.currentDataObj.electrodePitch - 0.5 ,'color','k');
    end
end

out=[];
for i=1:numel(outArg)
    eval(['out.(outArg{i})=' outArg{i} ';']);
end
