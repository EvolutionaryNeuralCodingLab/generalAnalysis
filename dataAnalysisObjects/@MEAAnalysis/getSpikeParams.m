function [data,hand]=getSpikeParams(obj,varargin)
% [data,hand]=getSpikeParams(obj,varargin)
% Function purpose : Get parameters related to spike shape
%
% Last updated : 20/10/16

obj.checkFileRecording;
parseObj = inputParser;

addParameter(parseObj,'smoothingDuration',0.5,@isnumeric); %[ms]
addParameter(parseObj,'WF',[]); %%WF can be either a string with the location of the data or the average highpass waveforms

addParameter(parseObj,'plotWidth',false,@isnumeric);
addParameter(parseObj,'plotWidthSpkNumbers',false,@isnumeric); %plot spike numbers for each point
addParameter(parseObj,'plotWidthAddNoise',true,@isnumeric); %add noise a noise with a max shift of the sampling frequency to eliminate quantization effects in plot
addParameter(parseObj,'plotWidthColorBar',true,@isnumeric); %add color bar
addParameter(parseObj,'plotWidthLogColormap',true,@isnumeric); %plot colorbar in log scale

addParameter(parseObj,'hand',[]);
addParameter(parseObj,'fieldPar',[]);
addParameter(parseObj,'saveFileName',[]);
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
if isempty(saveFileName)
    saveFileName=obj.files.(funName);
end

%load data if existing or contrinue with analysis
if exist(saveFileName,'file') & ~overwrite
    if nargout==1
        data=load(saveFileName);
    else
        disp('Spike parameter analysis already exists. use overwrite to run again');
    end
    return;
end

%populate grid sorter object
obj=populateGridSorterObj(obj);
if isempty(WF)
    avgSIF=obj.gridSorterObj.getSTData({'avgHPWF','neuronNames','nSpkTotal'});
    for j=1:size(avgSIF.avgHPWF,1),WF(j,:)=avgSIF.avgHPWF(j,avgSIF.neuronNames(1,j),:);end
else %WF can be either a string with the location of the data or the average highpass waveforms
    if isnumeric(WF)
        avgSIF=obj.gridSorterObj.getSTData({'neuronNames','nSpkTotal'});
    elseif isstr(WF)
        avgSIF=load(WF,'avgHPWF','neuronNames','nSpkTotal');
        WF=[];
        for j=1:size(avgSIF.avgHPWF,1),WF(j,:)=avgSIF.avgHPWF(j,avgSIF.neuronNames(1,j),:);end
    elseif isstruct(WF)
        avgSIF=WF;WF=[];
        for j=1:size(avgSIF.avgHPWF,1),WF(j,:)=avgSIF.avgHPWF(j,avgSIF.neuronNames(1,j),:);end
    end
end
nSpkTotal=avgSIF.nSpkTotal';
neuronNames=avgSIF.neuronNames;

preSpkSamples=obj.gridSorterObj.postPreFilteredWindow;
Fs=obj.currentDataObj.samplingFrequency;

nNeu=size(avgSIF.avgHPWF,1);
%% Main code
[nNeu,nSpikeSamples]=size(WF);

% calculate pre and post samples and their position in the waveform
preSpkSamples=preSpkSamples*Fs/1000;
postSpkSamples=nSpikeSamples-preSpkSamples;
pPre=1:preSpkSamples;
pPost=(preSpkSamples+1):(preSpkSamples+postSpkSamples);

%smooth spike shapes with local linear regression
smoothingSamples=round(smoothingDuration*Fs/1000);
spikeShapeSmooth=zeros([nNeu,nSpikeSamples]);
for i=1:nNeu
    spikeShapeSmooth(i,:) = smooth(WF(i,:),smoothingSamples,'loess');
end

%extract spike maxima
extermalAmp=max(abs(spikeShapeSmooth),[],2);
maxAmp=max(spikeShapeSmooth,[],2);
minAmp=min(spikeShapeSmooth,[],2);

%extract width from crossing at half amplitude
spkWidthSamples=zeros(nNeu,1);
widthVolt=(maxAmp+minAmp)/2;
for i=1:nNeu
    pDown=find(spikeShapeSmooth(i,pPre)>widthVolt(i),1,'last');
    pUp=preSpkSamples+find(spikeShapeSmooth(i,pPost)>widthVolt(i),1,'first');
    if ~(isempty(pDown) || isempty(pUp))
        spkWidthSamples(i)=(pUp-pDown);
    else
        spkWidthSamples(i)=NaN;
    end
end
width=spkWidthSamples/Fs*1000;

%calculate width from time difference between minimal and maximal derivative
spikeShapeDer=diff(spikeShapeSmooth,[],2);
[~,derMinP]=min(spikeShapeDer(:,pPre),[],2);
[~,derMaxP]=max(spikeShapeDer(:,pPost(1:end-1)),[],2);
widthDer=( (preSpkSamples+derMaxP)-derMinP)/Fs*1000;

if plotWidth
    if plotWidthAddNoise
        noiseAdd=(rand(nNeu,2)-0.5)/Fs*1000;
    else
        noiseAdd=zeros(nNeu,2);
    end
    
    f=figure;h.widthAxes=axes;
    if plotWidthLogColormap
        h.widthScatter=scatter(widthDer+noiseAdd(:,1),width+noiseAdd(:,2),25,log(extermalAmp));
    else
        h.widthScatter=scatter(widthDer+noiseAdd(:,1),width+noiseAdd(:,2),25,extermalAmp);
    end
    xlabel('Width (derivative)');
    ylabel('Width @ half maxima');
    
    if plotWidthSpkNumbers
        text(widthDer+noiseAdd(:,1),width+noiseAdd(:,2),num2str((1:nNeu)'),'fontSize',8);
    end
    
    if plotWidthColorBar
        h.widthCb=colorbar;
        if plotWidthLogColormap
            cbTicks=get(h.widthCb,'YTick');
            set(h.widthCb,'YTickLabel',num2str(exp(cbTicks)',3));
        end
        set(gca,'position',[0.1300    0.1100    0.7164    0.8150]);
        set(h.widthCb,'position',[0.8607    0.1071    0.0179    0.8167]);
        ylabel(h.widthCb,'Extermal Amp [\muV]');
    end

end

%save all data
fprintf('\nSaving results...');
save(saveFileName,'extermalAmp','maxAmp','minAmp','width','widthDer','neuronNames','nSpkTotal');

fprintf('Done!\n');

%save(saveFileName,'pMaxField','-append');

