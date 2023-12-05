function obj=spikeDetection(obj)
%Detect spikes in raw data and save waveforms
%Todo:
%Add detection of multiple spikes

%determine quantization
obj.detectionInt2uV=obj.detectionMaxSpikeAmp/2^(obj.detectionNQuantizationBits-1);

obj=obj.getHighpassFilter;

%start spike detection
fprintf('\nRunning spike detection on %s...',obj.dataRecordingObj.dataFileNames{1});

upSamplingFactor=obj.upSamplingFrequencySpike/obj.dataRecordingObj.samplingFrequency(1);
if upSamplingFactor~=round(upSamplingFactor) %check that upsampling factor is an integer
    upSamplingFactor=round(upSamplingFactor);
    obj.upSamplingFrequencySpike=upSamplingFactor*obj.dataRecordingObj.samplingFrequency(1);
    disp(['upSampling factor was not an integer and was rounded to: ' num2str(upSamplingFactor)]);
end

gaussianityWindow=obj.detectionGaussianityWindow/1000*obj.dataRecordingObj.samplingFrequency(1);
testSamples=gaussianityWindow*1000;

preSpikeSamples=obj.detectionPreSpikeWindow/1000*obj.dataRecordingObj.samplingFrequency(1); %must be > spikePeakInterval
postSpikeSamples=obj.detectionPostSpikeWindow/1000*obj.dataRecordingObj.samplingFrequency(1); %must be > spikePeakInterval
spikeTimeShiftIntervalSamples=obj.detectionSpikeTimeShiftInterval/1000*obj.dataRecordingObj.samplingFrequency(1);
postSpikeSamplesInitial=postSpikeSamples+spikeTimeShiftIntervalSamples;

timeVec=-preSpikeSamples:postSpikeSamplesInitial;
intrpTimeVec=timeVec(1):(1/upSamplingFactor):timeVec(end);
pZeroTimeVec=find(intrpTimeVec>=0,1,'first');

preSpikeSamplesIntrp=round(preSpikeSamples*upSamplingFactor); %must be > spikePeakInterval
postSpikeSamplesIntrp=round(postSpikeSamples*upSamplingFactor); %must be > spikePeakInterval
spikeTimeShiftIntervalIntrp=spikeTimeShiftIntervalSamples*upSamplingFactor; %must be > spikePeakInterval

minimumDetectionIntervalSamplesIntrp=obj.detectionMinimumDetectionInterval/1000*obj.dataRecordingObj.samplingFrequency(1)*upSamplingFactor;
peakDetectionSmoothingSamples=round(obj.detectionPeakDetectionSmoothingWindow/1000*obj.dataRecordingObj.samplingFrequency(1)*upSamplingFactor);
peakSmoothingKernel=fspecial('gaussian', [3*peakDetectionSmoothingSamples 1] ,peakDetectionSmoothingSamples);

%determine the chunck size
if obj.detectionMaxChunkSize>obj.dataRecordingObj.recordingDuration_ms
    startTimes=0;
    endTimes=obj.dataRecordingObj.recordingDuration_ms;
else
    startTimes=0:obj.detectionMaxChunkSize:obj.dataRecordingObj.recordingDuration_ms;
    endTimes=[startTimes(2:end)+obj.detectionChunkOverlap obj.dataRecordingObj.recordingDuration_ms];
end
nChunks=numel(startTimes);

obj.nCh=numel(obj.chPar.s2r);

matFileObj=cell(1,obj.nCh);
if obj.runWithoutSaving2File %if saving data is not required
    obj.runWithoutSaving2File=true;
    spikeShapesAll=cell(obj.nCh,nChunks);
    spikeTimesAll=cell(obj.nCh,nChunks);
else
    obj.runWithoutSaving2File=false;
    for i=find(obj.sortingFileNames.spikeDetectionExist==0)
        matFileObj{i} = matfile(obj.sortingFileNames.spikeDetectionFile{i},'Writable',true);
        matFileObj{i}.spikeShapes=zeros(preSpikeSamplesIntrp+postSpikeSamplesIntrp,0,obj.chPar.nValidChExt(i),'int16');
    end
end

%initiate arrays
Th=zeros(obj.nCh,nChunks);nCumSpikes=zeros(1,obj.nCh);
fprintf('\nExtracting spikes from chunks (total %d): ',nChunks);
for j=1:nChunks
    fprintf('%d ',j);
    %get data
    MAll=squeeze(obj.filterObj.getFilteredData(obj.dataRecordingObj.getData(obj.chPar.s2r(1:obj.nCh),startTimes(j),endTimes(j)-startTimes(j))))';
    if obj.detectionRemoveAllElectrodeMedian % if there are noise signals that are shared by all electrodes, like perfusion of other global artifacts
        MAll=bsxfun(@minus,MAll,median(MAll,2));
    end
    nSamples=size(MAll,1);
    for i=find(obj.sortingFileNames.spikeDetectionExist==0) %go over all channels that require rewriting
        %get local data
        Mlong=MAll(:,obj.chPar.surChExtVec{i});
        
        %estimate channel noise
        tmpData=buffer(Mlong(1:min(testSamples,nSamples),obj.chPar.pCenterCh(i)),gaussianityWindow,gaussianityWindow/2);
        noiseSamples=tmpData(:,kurtosis(tmpData,0)<obj.detectionKurtosisNoiseThreshold);
        noiseStd=std(noiseSamples(:));
        noiseMean=mean(noiseSamples(:)); 
        Th(i,j)=noiseMean-obj.detectionSpikeDetectionThresholdStd*noiseStd;
        
        %find thershold crossings and extract spike windows
        thresholdCrossings=find(Mlong(1:end-1,obj.chPar.pCenterCh(i))>Th(i,j) & Mlong(2:end,obj.chPar.pCenterCh(i))<Th(i,j));
        thresholdCrossings=thresholdCrossings(thresholdCrossings>preSpikeSamples & thresholdCrossings<nSamples-postSpikeSamplesInitial);
        %{
        plot(Mlong(:,obj.chPar.pCenterCh(i)));hold on;line([0 size(Mlong,1)],[Th(i,j) Th(i,j)],'color','r');plot(thresholdCrossings+1,Mlong(thresholdCrossings+1,obj.chPar.pCenterCh(i)),'og');
        %}
        %extract upsample and allign spikes
        startSamplesInM=[];startSamplesInIdx=[];
        if ~isempty(thresholdCrossings)
            %upsample
            startSamplesInM(1,1,:)=(0:nSamples:(nSamples*(obj.chPar.nValidChExt(i)-1)));
            idx=bsxfun(@plus,bsxfun(@plus,thresholdCrossings',(-preSpikeSamples:postSpikeSamplesInitial)'), startSamplesInM );
            M=Mlong(idx);
            
            %upsample data
            M = interp1(timeVec, M, intrpTimeVec, 'spline');
            nSamplesShort=size(M,1);
            
            %allign spike windows to spike extrema
            Msmooth = convn(M((pZeroTimeVec+1):(pZeroTimeVec+spikeTimeShiftIntervalIntrp),:,obj.chPar.pCenterCh(i)), peakSmoothingKernel, 'same');
            [spikeAmp,shift]=min(Msmooth);
            
            spikeTimesTmp=startTimes(j)+(thresholdCrossings'+shift/upSamplingFactor)/obj.dataRecordingObj.samplingFrequency(1)*1000; %[ms]
            
            nSamplesPerCh=numel(spikeTimesTmp)*nSamplesShort;
            startSamplesInIdx(1,1,:)=(0:nSamplesPerCh:nSamplesPerCh*obj.chPar.nValidChExt(i)-1);
            
            idx=bsxfun(@plus , bsxfun(@plus,pZeroTimeVec+shift+(0:nSamplesShort:(nSamplesPerCh-1)),(-preSpikeSamplesIntrp:(postSpikeSamplesIntrp-1))') , startSamplesInIdx);
            M=M(idx);
            %{
            figure;plotShifted(reshape(permute(M,[1 3 2]),[size(M,1)*size(M,3) size(M,2)]),'verticalShift',30);line([(obj.chPar.pCenterCh(i)-1)*size(M,1) obj.chPar.pCenterCh(i)*size(M,1)],[0 0],'color','g','lineWidth',3);
            ii=2;h=axes;activityTracePhysicalSpacePlot(h,obj.chPar.surChExtVec{i},squeeze(M(:,ii,:))',obj.chPar.rEn);
            %}
            if obj.detectionRemoveSpikesNotExtremalOnLocalGrid
                %check for a minimum (negative spike peak) over all channels to detect the channel with the strongest amplitude for each spike
                [maxV,maxP]=min(   min(    M((preSpikeSamplesIntrp-minimumDetectionIntervalSamplesIntrp):(preSpikeSamplesIntrp+minimumDetectionIntervalSamplesIntrp),:,:)   ,[],1)    ,[],3);
                p=(maxP==obj.chPar.pCenterCh(i));
                M=M(:,p,:);
                spikeTimesTmp=spikeTimesTmp(p'); %the p' is important to create a 1-0 empty matrix (not 0-1) which cell2mat can handle
            end
            if obj.detectionRemoveSpikesSubMinimumDelays
                p=find(diff(spikeTimesTmp)>obj.detectionMinimumDelayBetweenSpikes)+1;
                M=M(:,p,:);
                spikeTimesTmp=spikeTimesTmp(p'); %the p' is important to create a 1-0 empty matrix (not 0-1) which cell2mat can handle
            end
            spikeTimesAll{i,j}=spikeTimesTmp;
            if ~obj.runWithoutSaving2File
                tmpSpikeCount=numel(spikeTimesTmp);
                if numel(spikeTimesTmp)>0
                    nCumSpikes(i)=nCumSpikes(i)+numel(spikeTimesTmp);
                    matFileObj{i}.spikeShapes(:,(nCumSpikes(i)-tmpSpikeCount+1):nCumSpikes(i),:)=int16(M./obj.detectionInt2uV);
                end
            else
                spikeShapesAll{i,j}=int16(M./obj.detectionInt2uV);
            end
        else
            spikeTimesAll{i,j}=[];
            if obj.runWithoutSaving2File
                spikeShapesAll{i,j}=[];
            end
        end
    end
end
clear MAll Mlong M;

if ~obj.runWithoutSaving2File %write files to disk
    for i=find(obj.sortingFileNames.spikeDetectionExist==0)
        matFileObj{i}.Th=Th(i,:);
        
        matFileObj{i}.spikeTimes=cell2mat(spikeTimesAll(i,:));
        matFileObj{i}.preSpikeSamplesIntrp=preSpikeSamplesIntrp;
        matFileObj{i}.postSpikeSamplesIntrp=postSpikeSamplesIntrp;
        matFileObj{i}.upSamplingFrequencySpike=obj.upSamplingFrequencySpike;
        matFileObj{i}.minimumDetectionIntervalSamplesIntrp=minimumDetectionIntervalSamplesIntrp;
        matFileObj{i}.detectionInt2uV=obj.detectionInt2uV;
    end
else %keep files in memory
    obj.spikeDetectionData.spikeShapes=cell(1,obj.nCh);
    obj.spikeDetectionData.spikeTimes=cell(1,obj.nCh);
    for i=find(obj.sortingFileNames.spikeDetectionExist>0)
        obj.spikeDetectionData.spikeShapes{i}=cell2mat(spikeShapesAll(i,:));
        obj.spikeDetectionData.spikeTimes{i}=cell2mat(spikeTimesAll(i,:));
    end
    obj.spikeDetectionData.Th=Th;
    
    obj.spikeDetectionData.preSpikeSamplesIntrp=preSpikeSamplesIntrp;
    obj.spikeDetectionData.postSpikeSamplesIntrp=postSpikeSamplesIntrp;
    obj.spikeDetectionData.minimumDetectionIntervalSamplesIntrp=minimumDetectionIntervalSamplesIntrp;
    obj.spikeDetectionData.upSamplingFrequencySpike=obj.upSamplingFrequencySpike;
    obj.spikeDetectionData.detectionInt2uV=obj.detectionInt2uV;
end

obj=obj.findSortingFiles; %update sorted files
