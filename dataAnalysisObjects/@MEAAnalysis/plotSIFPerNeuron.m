function [hand,out]=plotSIFPerNeuron(obj,varargin)

obj.checkFileRecording;

parseObj = inputParser;

addParameter(parseObj,'outArg',{});

%input large matrices instead of calculating
addParameter(parseObj,'avgRawWF',[],@isnumeric);
addParameter(parseObj,'lowpassWFBaseline',[],@isnumeric);
addParameter(parseObj,'baselineSubstractedSIF',[],@isnumeric);
%addParameter(parseObj,'lowpassWF',[],@isnumeric);
%addParameter(parseObj,'neuronNames',[]);
addParameter(parseObj,'pNeu',[],@isnumeric); %the positions of the neurons to plot
addParameter(parseObj,'preSpikeMs',[]);
addParameter(parseObj,'saveFileName',[]);

addParameter(parseObj,'hand',[]);
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

obj.checkFileRecording(obj.files.getSpikePositionEstimation,'Spike position estimation file missing, please run getSpikePositionEstimation');

%populate grid sorter object
obj=populateGridSorterObj(obj);

%assign parameters
ch=obj.currentDataObj.channelNumbers;
En=obj.currentDataObj.chLayoutNumbers;
Fs=obj.currentDataObj.samplingFrequency(1);
preSpikeMs=obj.gridSorterObj.postPreRawWindow;

%load SIF variables
toLoad={'fieldPar','pMaxField'};
if isempty(baselineSubstractedSIF)
    toLoad=[toLoad 'baselineSubstractedSIF'];
end
%if isempty(lowpassWF)
%    toLoad=[toLoad 'lowpassWF'];
%end
if isempty(fileNameSpikeInducedFields)
    obj.checkFileRecording(obj.files.getSpikeInducedFields,'Spike induced fields file missing, please run getSpikeInducedFields');
    load(obj.files.getSpikeInducedFields,toLoad{:});
else
    load(fileNameSpikeInducedFields,toLoad{:});
end

if isempty(avgRawWF)
    if isempty(fileNameSTWaveform)
        if exist(obj.files.getSpikeTrigWF,'file')
            avgRawWF=load(obj.files.getSpikeTrigWF,'avgRawWF','neuronNames','par');
        else
            avgRawWF=obj.gridSorterObj.getSTData('avgRawWF','neuronNames');
            avgRawWF.par.preRawWindow=obj.gridSorterObj.postPreFilteredWindow;
            warning('Getting spike triggered data from grid sorter!!!! In the future, run getSpikeTrigWF from MEAAnalysis');
        end
    else
        avgRawWF=load(fileNameSTWaveform,'avgRawWF','neuronNames','P');
    end
    preSpikeMs=avgRawWF.par.preRawWindow;
    neuronNames=avgRawWF.neuronNames;
    avgRawWF=avgRawWF.avgRawWF;
end


if isempty(baselineSubstractedSIF)
    baselineSubstractedSIF=zeros(size(avgRawWF));
end

populateGridSorterObj(obj);
if ~exist('neuronNames','var')
    ic=obj.gridSorterObj.getSortedData('ic');
    neuronNames=ic.ic(1:2,:);
end

[Xc,Yc]=obj.currentDataObj.getElectrodePositions;
Xrs=fieldPar.X(1,:);Yrs=fieldPar.Y(1,:);

[nNeu,nCh,nSamples]=size(avgRawWF);
timeVec=(1:nSamples)/Fs*1000-preSpikeMs;

%% Plotting results
%temporary plot for debugging
IE=['?';'I';'E'];
pTmp=find(timeVec==0);
spikeMarker=ones(120,1)*nan(1,numel(timeVec));
spikeMarker(:,pTmp)=min(baselineSubstractedSIF(:));
spikeMarker(:,pTmp+1)=max(baselineSubstractedSIF(:));

minMaxXPos=[min(Xc) max(Xc)];
minMaxYPos=[min(Yc) max(Yc)];

%calculate max spike amplitude
spikePeakWidthMs=1;
spikePeakWidthSamples=spikePeakWidthMs*Fs/1000;
preSpikeSamples=preSpikeMs*Fs/1000;
maxSpikeAmp=max( abs(    bsxfun(@minus, avgRawWF(:,:,(preSpikeSamples-spikePeakWidthSamples/2):(preSpikeSamples+spikePeakWidthSamples/2)),...
    mean(avgRawWF(:,:,(preSpikeSamples-3*spikePeakWidthSamples/2):(preSpikeSamples-spikePeakWidthSamples)),3)  )) ,[],3);

%select neurons to plot
if isempty(pNeu)
    pNeu=1:nNeu;
end

cMap=[0.5 0.5 0.5;IC;EC];
f=figure('position',[10 50 1500 700]);
for i=pNeu
    hA(1)=subplot(2,5,[1 7]);
    [hPlot,scaleFac]=activityTracePhysicalSpacePlot(hA(1),obj.currentDataObj.channelNumbers,squeeze(avgRawWF(i,:,:)),En,'traceColor',[0.5 0.5 0.5],'gridLineWidth',0.5);hold on;
    activityTracePhysicalSpacePlot(hA(1),obj.currentDataObj.channelNumbers,squeeze(baselineSubstractedSIF(i,:,:)),En,'traceColor',[0 0 0],'scaleFac',scaleFac,'DrawElectrodeNumbers',0,'DrawGrid',0);
    %activityTracePhysicalSpacePlot(hA(1),obj.currentDataObj.channelNumbers,squeeze(lowpassWF(i,:,:)),En,'scaleFac',scaleFac,'DrawElectrodeNumbers',0,'DrawGrid',0);
    %activityTracePhysicalSpacePlot(hA(1),obj.currentDataObj.channelNumbers,spikeMarker,En,'scaleFac',scaleFac,'traceColor',[0.7 0.7 0.7],'gridLineWidth',0.5);
    tl=title(['Neu=' num2str(neuronNames(:,i)') ',idx=' num2str(i) ',Cls=' num2str(IE(fieldPar.classIE(i))) ',Mxch=' num2str(pMaxField(i)) ',PS=' num2str(fieldPar.polarityScore(i),2)]);
    tl.Color=cMap(fieldPar.classIE(i),:);
    hA(1).OuterPosition([1 3])=[-0.04 0.42];
    
    hA(2)=subplot(2,5,[3 9]);
    hA(2).Clipping='off';
    imagesc(fieldPar.XintGrid(1,:),fieldPar.YintGrid(:,1),squeeze(fieldPar.interpAll(:,:,i)));hold on;
    set(hA(2),'YDir','normal');
    hCB=colorbar('position',[ 0.7463    0.5943    0.0064    0.3314]);
    xlabel('[\mum]');
    ylabel('[\mum]');
    hA(2).OuterPosition([1 3])=[0.39 0.37];
    
    plot(Xc,Yc,'.g')
    hTmp=arrow([fieldPar.X(1,i);fieldPar.Y(1,i)]',[fieldPar.X(2,i);fieldPar.Y(2,i)],'Width',4);
    if exist('testPos','var')
        plot(testPos{i}(1,:),testPos{i}(2,:),'*r');
    end
    
    hA(3)=subplot(2,5,5);
    hA(3).OuterPosition([1 3])=[0.8 0.18];
    hCB2=IntensityPhysicalSpacePlot(ch,fieldPar.val(i,:),En,'h',hA(3),'plotElectrodeNumbers',0,'plotGridLines',0,'markerSize',50,'plotColorBar',0);hold on;
    hCB2=IntensityPhysicalSpacePlot(ch,maxSpikeAmp(i,:),En,'h',hA(3),'plotElectrodeNumbers',0,'plotGridLines',0,'markerSize',25);
    set(hCB2,'position',[0.97    0.7800    0.0051    0.1457],'YTick',[]);
    title('out=SIF , in=spk');
    
    hA(4)=subplot(2,5,10);
    [hPlot,scaleFac]=activityTracePhysicalSpacePlot(hA(4),obj.currentDataObj.channelNumbers,squeeze(avgRawWF(i,:,(preSpikeSamples-50):(preSpikeSamples+100))),En,'traceColor',[0 0 0],'gridLineWidth',0.5);
    hA(4).OuterPosition([3 4])=[0.23 0.5];
    title('Spike WF');
    
    drawnow update;
    pause;
    delete(hA);
end

out=[];
for i=1:numel(outArg)
    eval(['out.(outArg{i})=' outArg{i} ';']);
end
