function [data,hand]=getSpikeInducedFields(obj,varargin)
% [data,hand]=getSpikeInducedFields(obj,varargin)
% Function purpose : Calculate distribution of post spike fields (PSF)
%
% Last updated : 20/10/16

obj.checkFileRecording;

parseObj = inputParser;

addParameter(parseObj,'fieldPar',[]);

addParameter(parseObj,'electrodePitch',[],@isnumeric);
addParameter(parseObj,'nearestNeighborsDistance',190,@isnumeric);
addParameter(parseObj,'SIFStartMs',5,@isnumeric);
addParameter(parseObj,'SIFEndMs',30,@isnumeric);
addParameter(parseObj,'preSpikePeakMs',2,@isnumeric); %this can be larger since the exact spike time is defined by the algorithm
addParameter(parseObj,'postSpike0CrossLimMs',20,@isnumeric);
addParameter(parseObj,'medianFilterLengthMs',7,@isnumeric);
addParameter(parseObj,'spikePeakWidthMs',1,@isnumeric);
addParameter(parseObj,'badChannels',[]);
addParameter(parseObj,'usePostProcessing',false); %to use the old STAs from the postProcessingMethod

addParameter(parseObj,'smartInterSmoothness',0.0001,@isnumeric); %smoothing [0 1] - higher values fit is close to data (no low pass), 0.0000005 - more low pass
addParameter(parseObj,'weightFunctionStdMs',7,@isnumeric);
addParameter(parseObj,'maxPostSpikeWidthMs',3,@isnumeric);
addParameter(parseObj,'stdThresholdCrossingSpikeInitiation',4,@isnumeric);
addParameter(parseObj,'preSpikeMinInitiationMs',1.5,@isnumeric);
addParameter(parseObj,'preSpikeMaxInitiationMs',0.5,@isnumeric);
addParameter(parseObj,'postSpikeCorrMs',10,@isnumeric); %5
addParameter(parseObj,'postSpikeReturn2BaselineMs',50,@isnumeric); %
addParameter(parseObj,'baselineSubstructionMethod','preSpikeAverage',@isstr); % 'spikeOnsetOnLowpass','preSpikeAverage','interpolatedBaseline'
addParameter(parseObj,'useInterpolatedBaseline4StrongShifts',true,@isnumeric); %if true uses interplation only in cases of extreme change in baseline during field
addParameter(parseObj,'saveLowPassBaseline',false,@isnumeric);
addParameter(parseObj,'saveBaselineSubstractedSIF',false,@isnumeric);
addParameter(parseObj,'saveBaselineSubstractedRaw',false,@isnumeric);

addParameter(parseObj,'IEclassificationMethod','templateCorr',@isstr); %'templateCorr','SIFAmpPolarity';
addParameter(parseObj,'ECorrTh',0,@isnumeric);
addParameter(parseObj,'ICorrTh',0,@isnumeric);
addParameter(parseObj,'classIE',true,@isnumeric); %[true,false,vec]if false, all assumed inhibitory, can also be a vector with excitatory (3) and inhibitory (2) classifications (or 0 for require classification)

addParameter(parseObj,'PSFMethod','maxBaselineSubstracted',@isstr);% 'max','maxBaselineSubstracted','oldNormMax','integralBaselineSubstracted',
addParameter(parseObj,'fieldPositionMethod','interpolatedMaximaSpline',@isstr);%'maxima','interpolatedMaxima','COM','interpolatedMaximaSpline'

addParameter(parseObj,'fieldEdgeBand',10,@isnumeric); %the band arround the electrode array edge (from farthest electrode) to include in stats of field positions (if positive, takes also neurons outside the array)

addParameter(parseObj,'triangulateCellPosition',true,@isnumeric);

%variables calculated by the algorithm
addParameter(parseObj,'lowpassWF',[],@isnumeric);
addParameter(parseObj,'lowpassWFBaseline',[],@isnumeric);
addParameter(parseObj,'fileNameSTWaveform',[],@isstr);

%parameters for running on external data sets
addParameter(parseObj,'cellPosition',[],@isnumeric); % [2 x nNeurons] correction to position based on spike shape [um]
addParameter(parseObj,'avgRawWF',[],@isnumeric); %[nNeurons,nCh,nSamples] to give STA as direct input instead of loading
addParameter(parseObj,'neuronNames',[]);
addParameter(parseObj,'preSpikeMs',[]);
addParameter(parseObj,'saveFileName',[]);

%debuggin plots
addParameter(parseObj,'hand',[]); %handle for plots
addParameter(parseObj,'plotIEClass',0,@isnumeric);
addParameter(parseObj,'plotMaxWFAll',0,@isnumeric);
addParameter(parseObj,'plotLowPassFilters',false,@isnumeric);
addParameter(parseObj,'plotClassPerNeuron',false,@isnumeric);

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
        disp('spike induced fields analysis already exists. use overwrite to run again');
    end
    return;
end

%check that triangulation was already performed
if triangulateCellPosition && isempty(cellPosition)
    obj.checkFileRecording(obj.files.getSpikePositionEstimation,'Single cell position estimation file missing, please run getSpikePositionEstimation');
end

%populate grid sorter object
obj=populateGridSorterObj(obj);
if isempty(avgRawWF)
    disp('loading average spike triggered waveforms...');
    if isempty(fileNameSTWaveform)
        if exist(obj.files.getSpikeTrigWF,'file')
            load(obj.files.getSpikeTrigWF,'avgRawWF','neuronNames','par');
            preSpikeMs=par.preRawWindow;
        else
            if ~obj.gridSorterObj.sortingFileNames.STWaveformExist
                disp('Post processing data for spike sorting does not exist, please run spikePostProcessing method in grid sorter');
                return;
            else
                load(obj.gridSorterObj.sortingFileNames.STWaveformFile,'avgRawWF','neuronNames');
                preSpikeMs=obj.gridSorterObj.postPreRawWindow;
                warning('Getting spike triggered data from grid sorter!!!! In the future, run getSpikeTrigWF from MEAAnalysis');
            end
        end
    else
        avgRawWF=load(fileNameSTWaveform,'avgRawWF','neuronNames','par');
        preSpikeMs=avgRawWF.par.preRawWindow;
        neuronNames=avgRawWF.neuronNames;
        avgRawWF=avgRawWF.avgRawWF;
    end
    fprintf('done\n');
end
if isempty(neuronNames) %for old cases where neuronNames did not exist in STWaveformFile - Should be removed in the future
    load(obj.gridSorterObj.sortingFileNames.postProcessingAnalysisFile,'neuronNames');
end

%assign parameters
ch=obj.currentDataObj.channelNumbers;
En=obj.currentDataObj.chLayoutNumbers;
Fs=obj.currentDataObj.samplingFrequency(1);

%% Main code - general calculations
SIFStartSamples=SIFStartMs*Fs/1000;
SIFEndSamples=SIFEndMs*Fs/1000;
preSpikeSamples=preSpikeMs*Fs/1000;
preSpikePeakSamples=preSpikePeakMs*Fs/1000;
spikePeakWidthSamples=spikePeakWidthMs*Fs/1000;
maxPostSpikeWidthSamples=maxPostSpikeWidthMs*Fs/1000;
postSpikeCorrSamples=postSpikeCorrMs*Fs/1000;
weightFunctionStdSamples=weightFunctionStdMs*Fs/1000;
preSpikeMinInitiationSamples=preSpikeMinInitiationMs*Fs/1000;
preSpikeMaxInitiationSamples=preSpikeMaxInitiationMs*Fs/1000;
        
medianFilterSamples=round(medianFilterLengthMs*Fs/1000/2)*2+1; %has to be an odd number
postSpike0CrossLimSamples=postSpike0CrossLimMs*Fs/1000;
return2BaselineSamples=postSpikeReturn2BaselineMs*Fs/1000+preSpikeSamples;

[nNeurons,nCh,nSamples]=size(avgRawWF);
timeVec=(1:nSamples)/Fs*1000-preSpikeMs;

%Build inverse map between electrode and location
if ~isempty(electrodePitch)
    [Xc,Yc]=obj.currentDataObj.getElectrodePositions(electrodePitch);
else
    [Xc,Yc]=obj.currentDataObj.getElectrodePositions;
    electrodePitch=obj.currentDataObj.electrodePitch;
end

% get the channel with max spike for extimating spike remove segment (amplitude calculated relative to the segment just before the spike)
maxSpikeAmp=max( abs(    bsxfun(@minus, avgRawWF(:,:,(preSpikeSamples-spikePeakWidthSamples/2):(preSpikeSamples+spikePeakWidthSamples/2)),...
    mean(avgRawWF(:,:,(preSpikeSamples-3*spikePeakWidthSamples/2):(preSpikeSamples-spikePeakWidthSamples)),3)  )) ,[],3);
[~,pMaxSpikeElec]=max( maxSpikeAmp ,[],2);
for i=1:nCh
    pNeighbors{i}=find(sqrt((Xc-Xc(i)).^2+(Yc-Yc(i)).^2)<=nearestNeighborsDistance);
end

if strcmp(baselineSubstructionMethod,'interpolatedBaseline') || useInterpolatedBaseline4StrongShifts || saveLowPassBaseline
    calculateInterplatedBaseline=1;
else
    calculateInterplatedBaseline=0;
end

%% pre-process the input waveforms
if isempty(lowpassWF) || (isempty(lowpassWFBaseline) && calculateInterplatedBaseline)
    fprintf('Removing spikes from average WFs: ');

    lowpassWF=zeros(size(avgRawWF));
    if calculateInterplatedBaseline
        lowpassWFBaseline=zeros(size(avgRawWF));
    end
    
    pBaseline=1:(preSpikePeakSamples-preSpikeMinInitiationSamples);
    preExtension=round(preSpikePeakSamples); %extend the detection point by a few samples

    for i=1:nNeurons
        fprintf('%d,',i);
        pNeighborsMaxSpkElec=pNeighbors{pMaxSpikeElec(i)};
        if numel(pNeighborsMaxSpkElec)==1
            warning('There was only one nearest neighbohr used, consider changing ''nearestNeighborsDistance'' parameters!!!');
        end
        %extract the initiation segment right before the spike (on nearest neigbohrs) and remove the mean of the initial part of this segment so that all segments start at 0
        spikeInitiationWF=squeeze(avgRawWF(i,pNeighborsMaxSpkElec,(preSpikeSamples-preSpikePeakSamples):preSpikeSamples));%
        spikeInitiationWF=bsxfun(@minus,spikeInitiationWF,mean(spikeInitiationWF(:,pBaseline),2) );
        %calculate the spike onset (tr) according to where std increases rapidely over different electrode
        stdProfile=std(spikeInitiationWF);
        pSpikeOnset=min([preSpikePeakSamples-preSpikeMaxInitiationSamples,find(stdProfile > mean(stdProfile(pBaseline)) + stdThresholdCrossingSpikeInitiation*std(stdProfile(pBaseline)),1,'first')-1]);
        
        pSpikeStart(i)=(preSpikeSamples-preSpikePeakSamples+pSpikeOnset-preExtension);
        pSpikeSoftEnd=(preSpikeSamples+maxPostSpikeWidthSamples);
        
        %weights for slow synaptic potential extraction
        w1=ones(1,nSamples);
        w1(pSpikeStart(i) : pSpikeSoftEnd)=0;
        w1((pSpikeSoftEnd+1):(pSpikeSoftEnd+weightFunctionStdSamples*3))=1-exp(-( (1:weightFunctionStdSamples*3)/weightFunctionStdSamples).^2);
        lowpassWF(i,:,:) = csaps(1:nSamples,squeeze(avgRawWF(i,:,:)),smartInterSmoothness,1:nSamples,w1);
        
        %weights for baseline extraction
        if calculateInterplatedBaseline
            w2=ones(1,nSamples);
            w2(pSpikeStart(i) : return2BaselineSamples)=0;
            w2((return2BaselineSamples+1):end)=1-exp(-( (1:(nSamples-return2BaselineSamples))/weightFunctionStdSamples/2).^2);
            w2(1:pSpikeStart(i))=1-exp(-( (pSpikeStart(i):-1:1)/weightFunctionStdSamples/3).^2);
            lowpassWFBaseline(i,:,:) = csaps(1:nSamples,squeeze(lowpassWF(i,:,:)),1e-6,1:nSamples,w2);
        end
        %lowpassWFBaseline(i,:,:) = lowpassWF(i,:,:);
        
        %lowpassWFBaseline(i,:,(pSpikeStart(i)-50):1200) = interp1([1:(pSpikeStart(i)-50) 1200:2000],squeeze(lowpassWF(i,:,[1:(pSpikeStart(i)-50) 1200:2000]))',(pSpikeStart(i)-50):1200)';
        
        %plotting
        if plotLowPassFilters
            h(1)=subplot(2,3,1);
            plot(timeVec,squeeze(avgRawWF(i,pMaxSpikeElec(i),:)),'r');hold on;
            if calculateInterplatedBaseline
                plot(timeVec,squeeze(lowpassWFBaseline(i,pMaxSpikeElec(i),:)),'color',[0.5 0.5 0.5]);
            end
            plot(timeVec,squeeze(lowpassWF(i,pMaxSpikeElec(i),:)),'b');
            plot(timeVec,(w1-1)*50);
            plot(timeVec,(w2-1)*50);

            xlabel('Time [ms]');axis tight;
            title(['Max electrode-neuron: ' num2str(neuronNames(:,i)')]);
            legend({'WF','base','LP','W_{SIF}','W_{base}'},'Location','southeast','box','off');
            
            spikeZoom=squeeze(avgRawWF(i,pNeighbors{pMaxSpikeElec(i)},(preSpikeSamples-preSpikePeakSamples):(preSpikeSamples+SIFEndSamples)));%
            spikeZoom=bsxfun(@minus,spikeZoom,mean(spikeZoom(:,200),2) );
            h(2)=subplot(2,3,4);
            plot(timeVec((preSpikeSamples-preSpikePeakSamples):(preSpikeSamples+SIFEndSamples)),spikeZoom');
            xlabel('Time [ms]');axis tight;
            title('Spike Zoom- raw baseline sub');
            
            h(3)=subplot(2,3,[2 6]);
            [hPlot,scaleFac]=activityTracePhysicalSpacePlot(h(3),obj.currentDataObj.channelNumbers,squeeze(avgRawWF(i,:,:)),En,'traceColor','r','DrawElectrodeNumbers',1);hold on;
            [hPlot]=activityTracePhysicalSpacePlot(h(3),obj.currentDataObj.channelNumbers,squeeze(lowpassWF(i,:,:)),En,'scaleFac',scaleFac);hold on;
            if calculateInterplatedBaseline
                [hPlot]=activityTracePhysicalSpacePlot(h(3),obj.currentDataObj.channelNumbers,squeeze(lowpassWFBaseline(i,:,:)),En,'scaleFac',scaleFac,'traceColor',[0.5 0.5 0.5]);
            end
            
            pause;
            delete(h);
        end
    end
    fprintf(' done!\n');
end

%Remove waveforms of bad channels
if ~isempty(badChannels)
    [~,pBadCh]=intersect(ch,badChannels);
    [~,pGoodCh]=setdiff(ch,badChannels);
    lowpassWF(:,pBadCh,:)=NaN;
    lowpassWFBaseline(:,pBadCh,:)=NaN;
end
%% Baseline substruction
if strcmp(baselineSubstructionMethod,'spikeOnsetOnLowpass')
    %baselineSubstractedSIF=bsxfun(@minus,lowpassWF,lowpassWF(:,:,preSpikeSamples-preSpikePeakSamples)); %substracting a one point value is ok since it is done on the low-passed waveform
    baselineSubstractedSIF=bsxfun(@minus,lowpassWF,lowpassWF(:,:,preSpikeSamples)); %substracting a one point value is ok since it is done on the low-passed waveform
    baselineSubstractedRaw=bsxfun(@minus,avgRawWF,lowpassWF(:,:,preSpikeSamples)); %substracting a one point value is ok since it is done on the low-passed waveform
elseif strcmp(baselineSubstructionMethod,'preSpikeAverage')
    preBaseline=median(lowpassWF(:,:,(preSpikeSamples-preSpikePeakSamples):(preSpikeSamples-preSpikeMinInitiationSamples)),3);
    baselineSubstractedSIF=bsxfun(@minus,lowpassWF,preBaseline);
    baselineSubstractedRaw=bsxfun(@minus,avgRawWF,preBaseline);
elseif strcmp(baselineSubstructionMethod,'interpolatedBaseline')
    baselineSubstractedSIF=lowpassWF-lowpassWFBaseline;
    baselineSubstractedRaw=avgRawWF-lowpassWFBaseline;
end

%use lowpass interpolation method of baseline estimation for situations in which there is an extreme trend in the field over the duration of the sif
maxBaselineShift=20;
maxBaselineShiftFractionElecs=0.2;
if useInterpolatedBaseline4StrongShifts
    preBaseline=median(lowpassWF(:,:,1:200),3);
    postBaseline=median(lowpassWF(:,:,(end-200):end),3);
    pInterpolatedBaseline=mean(abs(preBaseline-postBaseline)>maxBaselineShift,2)>maxBaselineShiftFractionElecs; %if at least 20% of electrodes have a transition of more than 30uV in baseline
    baselineSubstractedSIF(pInterpolatedBaseline,:,:)=lowpassWF(pInterpolatedBaseline,:,:)-lowpassWFBaseline(pInterpolatedBaseline,:,:);
    baselineSubstractedRaw(pInterpolatedBaseline,:,:)=avgRawWF(pInterpolatedBaseline,:,:)-lowpassWFBaseline(pInterpolatedBaseline,:,:);
end

%set the main samples for field estimation based on time creteria
pRelevantSamples=(preSpikeSamples+SIFStartSamples):(preSpikeSamples+SIFEndSamples);
%pRelevantSamples=(preSpikeSamples):(preSpikeSamples+SIFEndSamples);


%% Calculate polarity score according to max field for classification

%build extended grid
nNeighbors=2;
[nRowsTmp,nColsTmp]=size(En);
EnExt=NaN(nRowsTmp+nNeighbors*2,nColsTmp+nNeighbors*2);
EnExt(1+nNeighbors:end-nNeighbors,1+nNeighbors:end-nNeighbors)=En;
pPeakElectrode=ceil((nNeighbors+2)^2/2);

%remove bad channels
if ~isempty(badChannels)
    [~,pBadChEnExt]=intersect(EnExt,badChannels);
    EnExt(pBadChEnExt)=NaN;
end

%find max amp electrode

%[~,pSpikeElec]=min(avgRawWF(:,:,preSpikeSamples+1),[],2);

T=load('SIFtemptates.mat');
%resample template to fit the times of the measured data
T.resTemplate = interp1(T.t_ms,T.templates',timeVec,'spline')';
%put zeros in resampled template if these times are outside the IESIFTemplate time limits
T.resTemplate(:,timeVec>T.t_ms(end) | timeVec<T.t_ms(1))=0;

T.pSIF=find(T.t_ms>=SIFStartMs & T.t_ms< SIFEndMs);
T.pResSIF=preSpikeSamples+(SIFStartSamples:(SIFEndSamples-1));
excludeSpikePeakElectrode=1;
polarityThreshold=0.25;
for i=1:nNeurons
    %find electrodes arround max spike electrode
    [pX,pY]=find(EnExt==pMaxSpikeElec(i));
    pElecs=EnExt(pX-nNeighbors:pX+nNeighbors,pY-nNeighbors:pY+nNeighbors); %get electrodes in extended grid
    if excludeSpikePeakElectrode
        pElecs(pPeakElectrode)=NaN;
    end
    pElecs=pElecs(~isnan(pElecs)); %remove NaNs
    nElecs=numel(pElecs);
    
    %for more accurate extraction of field polarity - have to be checked if working
    %{
            tmp=squeeze(baselineSubstractedSIF(i,pElecs,(preSpikeSamples+SIFStartSamples):(preSpikeSamples+SIFEndSamples)));
            pIntersection=findfirst(tmp(:,2:end)>0 & tmp(:,1:end-1)<0, 2, 1);
            pIntersection(pIntersection==0)=SIFEndSamples-SIFStartSamples;
            sortedIntersection=sort(pIntersection);
            SIFEndSamplesNew(i)=(preSpikeSamples+SIFStartSamples)+sortedIntersection(round(0.2*nElecs));
            tmp=tmp(:,1:(SIFEndSamplesNew(i)-(preSpikeSamples+SIFStartSamples)));
    %}
    tmp=squeeze(baselineSubstractedSIF(i,pElecs,pRelevantSamples));
    
    %plot(squeeze(baselineSubstractedSIF(i,pElecs,:))','b');hold on;plot(pRelevantSamples,squeeze(baselineSubstractedSIF(i,pElecs,pRelevantSamples))','r');
    %plot(squeeze(avgRawWF(i,pElecs,:))','b');hold on;plot(squeeze(baselineSubstractedSIF(i,pElecs,:))','r');
    %activityTracePhysicalSpacePlot([],obj.currentDataObj.channelNumbers,squeeze(avgRawWF(i,:,:)),En,'traceColor','k','DrawElectrodeNumbers',1);

    [absSIFscore,pTmp]=max(  abs(tmp)   ,[] , 2);
    [~,pOrder]=sort(absSIFscore);
        
    if strcmp(IEclassificationMethod,'templateCorr')
        [CTemp,pMax]=max(corr(T.resTemplate(:,T.pResSIF)',squeeze(baselineSubstractedRaw(i,pElecs,T.pResSIF))'));
        tmpPol=T.polarity(pMax);
        polaritySign(i)=sign(mean(tmpPol)+eps); %check the sign by majority vote relative to threshold
        polarityConf(i)=mean([CTemp(tmpPol==polaritySign(i)) -CTemp(tmpPol~=polaritySign(i))]);
    elseif strcmp(IEclassificationMethod,'SIFAmpPolarity')
        selectedElecs=pOrder(round((nElecs*0.5):end));%take only the high 50% of fields
        topPolarities=tmp(sub2ind(size(tmp), selectedElecs, pTmp(selectedElecs)));
        polarityConf(i)=mean(sign(topPolarities));
        polaritySign(i)=sign(polarityConf(i)-polarityThreshold); %check the sign at maximum and assign  -1 or 1 depending of whether positive or negative
    end
    
    %SIFtemptates
    if polaritySign(i)>0
        [~,pTmp]=max(max(tmp,[],2));
        pMaxField(i)=pElecs(pTmp);
        polarityScore(i)=max(tmp(pTmp,:));
    else
        [~,pTmp]=min(min(tmp,[],2));
        pMaxField(i)=pElecs(pTmp);
        polarityScore(i)=min(tmp(pTmp,:));
    end
    fields4Classification(i,:)=squeeze(baselineSubstractedSIF(i,pMaxField(i),:));
    fields4ClassificationRaw(i,:)=squeeze(baselineSubstractedRaw(i,pMaxField(i),:));
    
    %polarity for verification does not work well
    %[polarityValidity(i)]=mean(sign( mean(   abs(lowpassWF(i,pElecs,(preSpikeSamples+SIFStartSamples):(preSpikeSamples+SIFEndSamples)))-...
    %       abs(lowpassWFBaseline(i,pElecs,(preSpikeSamples+SIFStartSamples):(preSpikeSamples+SIFEndSamples)))   ,3)  ));
    
    %{
        h(1)=subplot(1,3,1:2);
        [hPlot,scaleFac]=activityTracePhysicalSpacePlot(h(1),obj.currentDataObj.channelNumbers,squeeze(avgRawWF(i,:,:)),En,'traceColor','r','DrawElectrodeNumbers',1);hold on;
        [hPlot]=activityTracePhysicalSpacePlot(h(1),obj.currentDataObj.channelNumbers,squeeze(lowpassWF(i,:,:)),En,'scaleFac',scaleFac);hold on;
        [hPlot]=activityTracePhysicalSpacePlot(h(1),obj.currentDataObj.channelNumbers,squeeze(lowpassWFBaseline(i,:,:)),En,'scaleFac',scaleFac,'traceColor',[0.5 0.5 0.5]);
        
        h(2)=subplot(1,3,3);
        %title(['polarity= ' num2str(polarityScore(i)), ' , Validity= ' num2str(polarityValidity(i))]);
        
        pause;
        delete(h);
    %}
end

fieldPar.fields4Classification=fields4Classification;
fieldPar.fields4ClassificationRaw=fields4ClassificationRaw;

fieldPar.samples4Classification=pRelevantSamples;

%% inhibitory excitatory classification
%determine which neurons to classify
%classes:  3 = excitatory, 2 = inhibitory, 1 = unclassified

if numel(classIE)==1
    if classIE==0 %do not classify, but set all to be inhibitory
       classIE=3*ones(1,nNeurons); 
    elseif classIE==1 %classify all
       classIE=ones(1,nNeurons);
    end %nothing happens for the case of one neuron in recording that was already clasified in the input
    toClassify=(classIE==1);
else
    toClassify=false(1,nNeurons);
end

pNotClassified=[];
if any(toClassify)
    useScore=0;
    polarityThresh=[0 2];
    polarityConfThresh=0.8;
    if useScore
        %pExcit=find(polarityScore<=polarityThresh(1) & polarityConf>=polarityConfThresh);
        %pInhib=find(polarityScore>=polarityThresh(2) & polarityConf>=polarityConfThresh);
        %pNotClassified=find((polarityScore>polarityThresh(1) & polarityScore<polarityThresh(2)) | polarityConf<polarityConfThresh);
    else
        pExcit=find(polaritySign==-1);
        pInhib=find(polaritySign==1);
        pNotClassified=[];
        %pNotClassified=find(polarityScore>=-polarityThresh & polarityScore<=polarityThresh);
    end
    fieldPar.polarityScore=polarityScore;
    fieldPar.polarityConf=polarityConf;
    fieldPar.polaritySign=polaritySign;
    
    fieldPar.classIE=ones(1,nNeurons);
    fieldPar.classIE(pExcit)=3; %excitatory
    fieldPar.classIE(pInhib)=2; %inhibitory
    fieldPar.classIE(~toClassify)=classIE(~toClassify); %give the neurons that should not be classified their original classification
    
else
    pExcit=find(classIE==3);
    pInhib=find(classIE==2);
    fieldPar.classIE=classIE;
end

if plotMaxWFAll
    
    %define number of subplots
    maxFields4Plot=375;
    n=ceil(sqrt(min(maxFields4Plot,nNeurons)/3/5));%define images in a 3 x 5 ratio
    xPlots=n*5;
    yPlots=n*3;
    nPlotPerPage=xPlots*yPlots;
    cMap=lines(2);
    cMap=[cMap;0 0 0;0 0 0];
    
    f=figure('Position',[50 50 1800 900],'Visible','off');
    for i=1:nNeurons
        h=subaxis(f,yPlots,xPlots,i,'S',0.001,'M',0.001);
        plot(timeVec,squeeze(avgRawWF(i,pMaxField(i),:)));hold on;
        plot(timeVec,squeeze(lowpassWF(i,pMaxField(i),:)),'r');axis tight;
        set(h,'XTickLabel',[],'YTick',[],'XTick',0,'TickLength',h.TickLength*5);
        %text(h.XLim(2),h.YLim(2),[num2str(neuronNames(1,i)) '-' num2str(neuronNames(2,i))],'VerticalAlignment','top','HorizontalAlignment','right');
        text(h.XLim(1),h.YLim(1),'*','color',cMap(4-fieldPar.classIE(i),:),'FontSize',18)
        text(h.XLim(2),h.YLim(1)+0.2*diff(h.YLim),[num2str(i) ',' num2str(neuronNames(1,i)) '-' num2str(neuronNames(2,i))],'VerticalAlignment','top','HorizontalAlignment','right');
    end
    f.Visible='on';
end
%{
        for i=1:nNeurons
        f=figure;
        h=axes;
        [hPlot,scaleFac]=activityTracePhysicalSpacePlot(h,1:120,squeeze(avgRawWF(i,:,:)),En,'traceColor','r','DrawElectrodeNumbers',1);hold on;
        [hPlot]=activityTracePhysicalSpacePlot(h,1:120,squeeze(lowpassWF(i,:,:)),En,'scaleFac',scaleFac);
        title(['neuron ' num2str(neuronNames(:,i)') ', class = ' num2str(fieldPar.classIE(i))]);
        pause;
        delete(f);
        end
%}

%% calculate post spike fields
fprintf('\nCalculating PSDs...');
% include the fact that each neuron has a different end time for integration
% pRelevantSamples=(preSpikeSamples+SIFStartSamples):SIFEndSamplesNew(i);

switch PSFMethod
    case 'max' %pre value substracted
        fieldPar.val(pInhib,:)=max(lowpassWF(pInhib,:,pRelevantSamples),[],3);
        fieldPar.val(pExcit,:)=-min(lowpassWF(pExcit,:,pRelevantSamples),[],3);
        fieldPar.val(pNotClassified,:)=max(abs(lowpassWF(pNotClassified,:,pRelevantSamples)),[],3);
        
    case 'oldNormMax'
        %peak voltage normalized by pre spike peak
        %fieldPar.val(pInhib,:)=max(lowpassWF(pInhib,:,pRelevantSamples),[],3)-mean(lowpassWF(pInhib,:,1:(preSpikeSamples-preSpikePeakSamples)),3);
        %fieldPar.val(pExcit,:)=-(min(lowpassWF(pExcit,:,(1+preSpikeSamples):(preSpikeSamples+postSpike0CrossLimSamples)),[],3)-mean(lowpassWF(pExcit,:,1:(preSpikeSamples-preSpikePeakSamples)),3));
        
        fieldPar.val(pInhib,:)=max(lowpassWF(pInhib,:,pRelevantSamples),[],3)-mean(lowpassWF(pInhib,:,(preSpikeSamples-spikePeakWidthSamples/2):(preSpikeSamples+spikePeakWidthSamples/2)),3);
        fieldPar.val(pExcit,:)=-(min(lowpassWF(pExcit,:,(1+preSpikeSamples):(preSpikeSamples+postSpike0CrossLimSamples)),[],3)-mean(lowpassWF(pExcit,:,(preSpikeSamples-spikePeakWidthSamples/2):(preSpikeSamples+spikePeakWidthSamples/2)),3));
        fieldPar.val(pNotClassified,:)=max(abs(lowpassWF(pNotClassified,:,pRelevantSamples),[],3)-mean(lowpassWF(pNotClassified,:,(preSpikeSamples-spikePeakWidthSamples/2):(preSpikeSamples+spikePeakWidthSamples/2)),3));
     
    case 'maxBaselineSubstracted'   
         
        fieldPar.val(pInhib,:)=max(baselineSubstractedSIF(pInhib,:,pRelevantSamples),[],3);
        fieldPar.val(pExcit,:)=-min(baselineSubstractedSIF(pExcit,:,pRelevantSamples),[],3);
        fieldPar.val(pNotClassified,:)=max(abs(baselineSubstractedSIF(pNotClassified,:,pRelevantSamples)),[],3);
        
    case 'integralBaselineSubstracted' %!!!! Has to be rewritten to support separation between excitatory and inhibitory
        %mean voltage normalized by pre spike mean
        fieldPar.val(pInhib,:)=mean(baselineSubstractedSIF(pInhib,:,pRelevantSamples),3); %for inhibitory cells
        fieldPar.val(pExcit,:)=-mean(baselineSubstractedSIF(pExcit,:,pRelevantSamples),3); %for inhibitory cells
        fieldPar.val(pNotClassified,:)=mean(baselineSubstractedSIF(pNotClassified,:,pRelevantSamples),3);
          
    otherwise
        error('SIF calculation method not valid');
end

if plotClassPerNeuron
    
    IE=['?';'I';'E'];
    pTmp=find(timeVec==0);
    spikeMarker=ones(numel(obj.currentDataObj.channelNumbers),1)*nan(1,numel(timeVec));
    spikeMarker(:,pTmp)=min(lowpassWF(:));
    spikeMarker(:,pTmp+1)=max(lowpassWF(:));
    
    for i=1:nNeurons,
        h1=subplot(3,4,[1 11]);
        [hPlot,scaleFac]=activityTracePhysicalSpacePlot(h1,obj.currentDataObj.channelNumbers,squeeze(avgRawWF(i,:,:)),En,'traceColor','r','averageSubstruction',1);hold on;
        activityTracePhysicalSpacePlot(h1,obj.currentDataObj.channelNumbers,squeeze(baselineSubstractedSIF(i,:,:)),En,'scaleFac',scaleFac,'DrawElectrodeNumbers',1,'averageSubstruction',1);
        %activityTracePhysicalSpacePlot(h1,obj.currentDataObj.channelNumbers,squeeze(baselineSubstractedSIF(i,:,pRelevantSamples)),En,'scaleFac',scaleFac,'traceColor','k');
        activityTracePhysicalSpacePlot(h1,obj.currentDataObj.channelNumbers,spikeMarker,En,'scaleFac',scaleFac,'DrawElectrodeNumbers',1,'traceColor',[0.7 0.7 0.7]);
        title(['Neuron=' num2str(neuronNames(:,i)') 'index=' num2str(i) ', Max ch=' num2str(pMaxField(i)) ', C=' num2str(fieldPar.polaritySign(i))]);
        h2=subplot(3,4,8);hCB=IntensityPhysicalSpacePlot(ch,fieldPar.val(i,:),En,'h',h2,'plotElectrodeNumbers',0,'markerSize',30);
        title(IE(fieldPar.classIE(i)));
        pause;
        delete([h1 h2]);
    end
    
end

makeGaussianFit=0;
if makeGaussianFit
    gaussFit.mX=zeros(1,nNeurons);
    gaussFit.mY=zeros(1,nNeurons);
    gaussFit.sX=zeros(1,nNeurons);
    gaussFit.sY=zeros(1,nNeurons);
    gaussFit.A=zeros(1,nNeurons);
    gaussFit.Theta=zeros(1,nNeurons);
    for i=1:nNeurons
        [fitresult] = fmgaussfit(Xc,Yc,fieldPar.val(i,:)); %[amp, ang, sx, sy, xo, yo, zo]
        gaussFit.A(i)=fitresult(1);
        gaussFit.Theta(i)=fitresult(2);
        gaussFit.sX(i)=fitresult(3);
        gaussFit.sY(i)=fitresult(4);
        gaussFit.mX(i)=fitresult(5);
        gaussFit.mY(i)=fitresult(6);
    end
end

%calculate edge neurons
[~,pMax]=max(fieldPar.val,[],2);
[m,n]=size(En);
fieldPar.edgeNeurons=zeros(1,nNeurons);
for i=1:nNeurons
    [pX,pY]=find(En==neuronNames(1,i));
    if isempty(pX)
        disp('Warning!!! neuron name not found!');
    end
    if pX==1 || pX==n || pY==1 || pY==m
        fieldPar.edgeNeurons(i)=1;
    else
        surroundingSquare=En(pY-1:pY+1,pX-1:pX+1);
        if any(any(isnan(surroundingSquare)))
            fieldPar.edgeNeurons(i)=2;
        end
    end
end

%calculate a band arround electrode area to detect points outside array
%calculate points outside electrode area
[pBound] = boundary(Xc',Yc'); %Calculate bounding points
mX=mean([Xc;Yc],2);
[teta,r]=cart2pol((Xc(pBound)-mX(1)),Yc(pBound)-mX(2));

rSIF=r+fieldEdgeBand;
[xB,yB]=pol2cart(teta,rSIF);
confiningPolgon=[xB+mX(1);yB+mX(2)];

fprintf('\nCalculating field peak...');
switch fieldPositionMethod
    case 'interpolatedMaximaSpline' %fits a 2D polynomial on a local grid of 9 points surrounding center
        dXinterp=10; %um
        [m,n]=size(En);
        Z=nan([m,n]);
        %Z=zeros([m,n]);
        fieldCoord=zeros(2,nNeurons);
        Zidx=sub2ind([m,n],Xc(ch)/electrodePitch,Yc(ch)/electrodePitch);
        [Xmesh,Ymesh]=meshgrid(unique(Yc),unique(Xc));
        intX=min(Xc):dXinterp:max(Xc);
        intY=min(Yc):dXinterp:max(Yc);
        [fieldPar.XintGrid,fieldPar.YintGrid]=meshgrid(intY,intX);
        fieldPar.interpAll=zeros(size(fieldPar.XintGrid,1),size(fieldPar.YintGrid,2),nNeurons);
        notInInArray = ~inpolygon(fieldPar.XintGrid,fieldPar.YintGrid,confiningPolgon(1,:),confiningPolgon(2,:));

        warning('off', 'MATLAB:interp2:NaNstrip'); %suppress warning related to the existance of nan values
        for i=1:nNeurons
            F = scatteredInterpolant(Xc', Yc',fieldPar.val(i,:)');
            fieldPar.interpAll(:,:,i) = interp2(Xmesh,Ymesh,F(Xmesh,Ymesh),fieldPar.XintGrid,fieldPar.YintGrid,'spline');
            tmp=squeeze(fieldPar.interpAll(:,:,i)');
            tmp(notInInArray)=0;
            [~,pMax] = max(tmp(:));
            [fieldCoord(1,i),fieldCoord(2,i)] = ind2sub(size(tmp),pMax);
            %p = polyFit2D(Z,XGrid,YGrid,2,2);f = polyVal2D(p,XGrid,YGrid,2,2);imagesc(f)
        end
        warning('on', 'MATLAB:interp2:NaNstrip');
        fieldPar.Xfield=intX(fieldCoord(1,:));
        fieldPar.Yfield=intY(fieldCoord(2,:));
        
    case 'interpolatedMaxima' %fits a 2D polynomial on a local grid of 9 points surrounding center
        [m,n]=size(En);
        Z=nan([m,n]);
        %Z=zeros([m,n]);
        fieldCoord=zeros(2,nNeurons);
        Zidx=sub2ind([m,n],Xc(ch)/electrodePitch,Yc(ch)/electrodePitch);
        [Xmesh,Ymesh]=meshgrid(unique(Yc),unique(Xc));
        for i=1:nNeurons
            F = scatteredInterpolant(Xc', Yc',fieldPar.val(i,:)');
            Z = interp2(Xmesh, Ymesh,F(Xmesh,Ymesh),Xmesh,Ymesh,'spline');
            [fieldCoord(:,i)] = peakfit2d(Z');
            %p = polyFit2D(Z,XGrid,YGrid,2,2);f = polyVal2D(p,XGrid,YGrid,2,2);imagesc(f)
        end
        fieldPar.Xfield=fieldCoord(1,:)*electrodePitch;
        fieldPar.Yfield=fieldCoord(2,:)*electrodePitch;
        
    case 'medianCOM' %biased by array edges
        %pTmp=fieldPar.val>median(fieldPar.val,2)*ones(1,nCh);
        medSubstractedField=fieldPar.val-(median(fieldPar.val,2)*ones(1,nCh));
        fieldPar.Xfield=(sum(bsxfun(@times,medSubstractedField,Xc),2)./sum(medSubstractedField,2))';
        fieldPar.Yfield=(sum(bsxfun(@times,medSubstractedField,Yc),2)./sum(medSubstractedField,2))';
        
    case 'maxima'
        [PSF,pChPSF]=max(fieldPar.val,[],2);%location of field integral maxima
        fieldPar.Xfield=Xc(ch(pChPSF));
        fieldPar.Yfield=Yc(ch(pChPSF));
        
    case 'fitGaussian'
        [m,n]=size(En);
        Z=nan([m,n]);
        %Z=zeros([m,n]);
        fieldCoord=zeros(2,nNeurons);
        [YGrid,XGrid]=meshgrid(1:size(Z,1),1:size(Z,2));
        for i=1:nNeurons
            Z(sub2ind([m,n],Xc(ch)/electrodePitch,Yc(ch)/electrodePitch))=fieldPar.val(i,:);
            [fitresult] = fmgaussfit(XGrid,YGrid,Z);
            fieldCoord(:,i) = fitresult([5 6]);
        end
        fieldPar.Xfield=fieldCoord(1,:)*electrodePitch;
        fieldPar.Yfield=fieldCoord(2,:)*electrodePitch;
        
    case 'sumOfRegMax'
        [m,n]=size(En);
        %Z=min(fieldPar.val(:))*ones([m+2,n+2]);
        %Z0=min(fieldPar.val(:))*ones([m,n]);
        Z=zeros([m+2,n+2]);
        Z0=zeros([m,n]);
        fieldCoord=zeros(2,nNeurons);
        [YGrid,XGrid]=meshgrid(1:size(Z,1),1:size(Z,2));
        for i=1:nNeurons
                %Z0(sub2ind([m,n],Xc(ch)/electrodePitch,Yc(ch)/electrodePitch))=fieldPar.val(i,:);
                Z0(sub2ind([m,n],Xc(ch)/electrodePitch,Yc(ch)/electrodePitch))=fieldPar.val(i,:)-min(fieldPar.val(i,:));
            Z(2:end-1,2:end-1)=Z0;
            [ind] = find(imregionalmax(Z,8));
            pTmp=find(Z(ind)>(fieldPar.val(i,pMaxField(i))/2));
            nPeaks=numel(pTmp);
            ys=zeros(nPeaks,1);xs=zeros(nPeaks,1);
            for j=1:nPeaks
                K=Z((XGrid(ind(pTmp(j)))-1):(XGrid(ind(pTmp(j)))+1),(YGrid(ind(pTmp(j)))-1):(YGrid(ind(pTmp(j)))+1));
                % approximate polynomial parameter
                a = (K(2,1)+K(1,1)-2*K(1,2)+K(1,3)-2*K(3,2)-2*K(2,2)+K(2,3)+K(3,1)+K(3,3));
                b = (K(3,3)+K(1,1)-K(1,3)-K(3,1));
                c = (-K(1,1)+K(1,3)-K(2,1)+K(2,3)-K(3,1)+K(3,3));
                %d = (2*K(2,1)-K(1,1)+2*K(1,2)-K(1,3)+2*K(3,2)+5*K(2,2)+2*K(2,3)-K(3,1)-K(3,3));
                e = (-2*K(2,1)+K(1,1)+K(1,2)+K(1,3)+K(3,2)-2*K(2,2)-2*K(2,3)+K(3,1)+K(3,3));
                f = (-K(1,1)-K(1,2)-K(1,3)+K(3,1)+K(3,2)+K(3,3));
                
                % (ys,xs) is subpixel shift of peak location relative to point (2,2)
                xs(j) = (6*b*c-8*a*f)/(16*e*a-9*b^2);
                ys(j) = (6*b*f-8*e*c)/(16*e*a-9*b^2);
            end
            fieldCoord(:,i)=[mean(XGrid(ind(pTmp))-1+xs);mean(YGrid(ind(pTmp))-1+ys)];
            testPos{i}=[XGrid(ind(pTmp))-1+xs YGrid(ind(pTmp))-1+ys]'*electrodePitch;
            %testPos{i}=[XGrid(ind(pTmp))-1 YGrid(ind(pTmp))-1]'*electrodePitch;
        end
        fieldPar.Xfield=fieldCoord(1,:)*electrodePitch;
        fieldPar.Yfield=fieldCoord(2,:)*electrodePitch;
        
        %s = regionprops(L, 'Centroid');
end

%incorporate cell position
if triangulateCellPosition && isempty(cellPosition)
    cellPosEst=load(obj.files.getSpikePositionEstimation,'');
    cellPosition(1,:)=cellPosEst.Xrs;
    cellPosition(2,:)=cellPosEst.Yrs;
    %{
        figure;
        for i=1:nNeurons
        h=axes;
        [hPlot,scaleFac]=activityTracePhysicalSpacePlot(h,1:120,squeeze(avgSpkWF(i,:,:)),En,'traceColor','b','DrawElectrodeNumbers',1);hold on;title(neuronNames(:,i));
        pause;
        delete(h);
        end
    %}
elseif ~triangulateCellPosition && isempty(cellPosition) %use max spike electrode as a cell position estimator
    cellPosition(1,:)=Xc(neuronNames(1,:));
    cellPosition(2,:)=Yc(neuronNames(1,:));
end
%create projection vectors
fieldPar.X=[cellPosition(1,:);fieldPar.Xfield];
fieldPar.Y=[cellPosition(2,:);fieldPar.Yfield];

pTmp=isnan(fieldPar.Xfield);
fieldPar.X(:,pTmp)=NaN;
fieldPar.Y(:,pTmp)=NaN;

if size(cellPosition,2)==1 && size(cellPosition,1)~=2
    error('cellPosition was not entered in the correct format');
end

fieldPar.mag=sqrt((fieldPar.X(2,:)-fieldPar.X(1,:)).^2 + (fieldPar.Y(2,:)-fieldPar.Y(1,:)).^2);
fieldPar.angle=atan2(fieldPar.Y(2,:)-fieldPar.Y(1,:),fieldPar.X(2,:)-fieldPar.X(1,:));

%save all data
fprintf('\nSaving results...');
save(saveFileName,'pMaxField','lowpassWF','fieldPar','timeVec','par','-v7.3');
if saveLowPassBaseline
    save(saveFileName,'lowpassWFBaseline','-append');
end
if saveBaselineSubstractedSIF
    save(saveFileName,'baselineSubstractedSIF','-append');
end
if saveBaselineSubstractedRaw
    save(saveFileName,'baselineSubstractedRaw','-append');
end
fprintf('Done!\n');

%save(saveFileName,'pMaxField','-append');

