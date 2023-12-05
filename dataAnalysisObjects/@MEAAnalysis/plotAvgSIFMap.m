function [hand,out]=plotAvgSIFMap(obj,varargin)

obj.checkFileRecording;

parseObj = inputParser;

addParameter(parseObj,'hand',[]);
addParameter(parseObj,'pNeu',[],@isnumeric); %the positions of the neurons to plot
addParameter(parseObj,'gridRes',10,@isnumeric); %um
addParameter(parseObj,'gridRadius',500,@isnumeric); %um
addParameter(parseObj,'equalWeight4EveryNeuron',0,@isnumeric);
addParameter(parseObj,'polarPlotDistanceRange',[0 500],@isnumeric);
addParameter(parseObj,'outArg',{});


addParameter(parseObj,'angle2LateralRight',0,@isnumeric); %The angle in degrees to rotate the axes so that lateral is on the left
addParameter(parseObj,'reverseDirction',0,@isnumeric); %whether flip direction is needed to make anterior to top when lateral is left
addParameter(parseObj,'spikeEdgeBand',0,@isnumeric); %the band arround the electrode array edge (from farthest electrode) to include in stats of neuron positions (if positive, takes also neurons outside the array)
addParameter(parseObj,'fieldEdgeBand',-50,@isnumeric); %the band arround the electrode array edge (from farthest electrode) to include in stats of field positions (if positive, takes also neurons outside the array)

addParameter(parseObj,'plotAvgField',true,@isnumeric);
addParameter(parseObj,'contoursPos',0.4:0.1:1,@isnumeric); %um
addParameter(parseObj,'EC',[0 0.2 0.8],@isnumeric);
addParameter(parseObj,'IC',[0.95 0.1 0.05],@isnumeric);

addParameter(parseObj,'fileNameSpikeInducedFields',[],@isstr);
addParameter(parseObj,'fileNameSTWaveform',[],@isstr);
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

%populate grid sorter object
%obj=populateGridSorterObj(obj);

%assign parameters
ch=obj.currentDataObj.channelNumbers;
En=obj.currentDataObj.chLayoutNumbers;
%Fs=obj.currentDataObj.samplingFrequency(1);
%preSpikeMs=obj.gridSorterObj.postPreRawWindow;

%[Xc,Yc]=obj.currentDataObj.getElectrodePositions;
%load(obj.files.getSpikePositionEstimation,'Xrs','Yrs','Zrs');

if isempty(fileNameSpikeInducedFields)
    obj.checkFileRecording(obj.files.getSpikeInducedFields,'Spike induced fields file missing, please run getSpikeInducedFields');
    load(obj.files.getSpikeInducedFields,'fieldPar','pMaxField');
else
    load(fileNameSpikeInducedFields,'fieldPar','pMaxField');
end
nNeu=numel(fieldPar.mag);

%select sub population
if isempty(pNeu)
    pNeuron=true(1,nNeu);
else
    pNeuron=false(1,nNeu);
    pNeuron(pNeu)=true;
end

%calculate points outside electrode area
[Xc,Yc]=obj.currentDataObj.getElectrodePositions;
Xrs=fieldPar.X(1,:);Yrs=fieldPar.Y(1,:);

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

pI=pNeuron & fieldPar.classIE==2 & pDist & ~fieldPar.edgeNeurons;
pE=pNeuron & fieldPar.classIE==3 & pDist & ~fieldPar.edgeNeurons;

xx=-gridRadius:gridRes:gridRadius;
[X,Y]=meshgrid(xx,xx);

[theta,r]=cart2pol(X,Y);
rotAngRad=-angle2LateralRight/180*pi; %rotate sampling point to the opposite dirction of the angle
if reverseDirction
    theta=-theta;
end
theta=theta+rotAngRad;
[Xrot,Yrot]=pol2cart(theta,r);

%rotMat=[cos(rotAngRad) -sin(rotAngRad);sin(rotAngRad) cos(rotAngRad)];
%tmpX=rotMat*[mX(:)';mY(:)'];
intAll=nan(size(X,1),size(Y,2),nNeu);
for i=1:nNeu
    intAll(:,:,i)=interp2(fieldPar.XintGrid-fieldPar.X(1,i), fieldPar.YintGrid-fieldPar.Y(1,i),squeeze(fieldPar.interpAll(:,:,i)),Xrot,Yrot);
    %imagesc(fieldPar.XintGrid(1,:),fieldPar.YintGrid(1,:),squeeze(fieldPar.interpAll(:,:,1)));
end

if equalWeight4EveryNeuron
    mn=min(min(intAll,[],1),[],2);
    mx=max(max(intAll,[],1),[],2);
    intAll=bsxfun(@rdivide,bsxfun(@minus,intAll,mn),mx-mn);
end

avgFieldI=squeeze(nanmean(intAll(:,:,pI),3));
avgFieldE=squeeze(nanmean(intAll(:,:,pE),3));

nAvgFieldI=(avgFieldI-min(avgFieldI(:)))/(max(avgFieldI(:))-min(avgFieldI(:)));
nAvgFieldE=(avgFieldE-min(avgFieldE(:)))/(max(avgFieldE(:))-min(avgFieldE(:)));

if plotAvgField
    hand.fig=figure('position',[200 200 1000 400]);
    hand.ax(1)=subplot(1,2,1);
    imagesc(xx,xx,nAvgFieldE);hold on;
    line([0 0],[-gridRadius gridRadius],'color','r');
    line([-gridRadius gridRadius],[0 0],'color','r');
    contour(xx,xx,nAvgFieldE,contoursPos,'k');
    [~,pMax]=max(nAvgFieldE(:));
    [maxY,maxX]=ind2sub(size(nAvgFieldE),pMax);
    plot(xx(maxX),xx(maxY),'.k');
    
    set(hand.ax(1),'YDir','normal');
    axis equal;
    xlabel('\mum');ylabel('\mum');
    title('Excitatory','color',EC);
    
    hand.ax(2)=subplot(1,2,2);
    imagesc(xx,xx,nAvgFieldI);hold on;
    line([0 0],[-gridRadius gridRadius],'color','r');
    line([-gridRadius gridRadius],[0 0],'color','r');
    contour(xx,xx,nAvgFieldI,contoursPos,'k');
    [~,pMax]=max(nAvgFieldI(:));
    [maxY,maxX]=ind2sub(size(nAvgFieldI),pMax);
    plot(xx(maxX),xx(maxY),'.k');
    
    set(hand.ax(2),'YDir','normal');
    axis equal;
    xlabel('\mum');ylabel('\mum');
    title('Inhibitory','color',IC);
end

out=[];
for i=1:numel(outArg)
    eval(['out.(outArg{i})=' outArg{i} ';']);
end
