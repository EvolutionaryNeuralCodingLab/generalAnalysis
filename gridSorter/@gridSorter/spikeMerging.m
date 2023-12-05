function [obj,avgWF,stdWF,ch,mergedNeurons]=spikeMergingClusters(obj)
%Merge waveforms belonging to the same neurons which were either not clustered correctly or were associated with different channels
avgWF=cell(1,obj.nCh);
stdWF=cell(1,obj.nCh);
ch=cell(1,obj.nCh);
isNoise=cell(1,obj.nCh);
if isempty(obj.sortingDir)
    obj.runWithoutSaving2File=true;
else
    obj.runWithoutSaving2File=false;
end

%steepness=1.2;
%nonLin = @(x) x+x./sqrt(1+abs(x).^steepness);
%x=-100:100;plot(x,nonLin(x));

load(obj.sortingFileNames.spikeDetectionFile{1},'postSpikeSamplesIntrp','preSpikeSamplesIntrp');
nSamples=postSpikeSamplesIntrp+preSpikeSamplesIntrp;
preSpikeSamples4NoiseDetection=obj.mergingPreSpike4NoiseDetection*preSpikeSamplesIntrp;
postSpikeSamples4NoiseDetection=obj.mergingPostSpike4NoiseDetection*postSpikeSamplesIntrp;
pSamples4NoiseDetection=round(nSamples/2-preSpikeSamples4NoiseDetection):round(nSamples/2+postSpikeSamples4NoiseDetection);

maxSpikeShiftSamples=obj.mergingMaxSpikeShift*obj.upSamplingFrequencySpike/1000;

load(obj.sortingFileNames.avgWaveformFile,'avgClusteredWaveforms','stdClusteredWaveforms','nAvgSpk');

%classify tempates as noise if the number of points exceeding N std is smaller than a threshold
validSamplesInTemplate=cell(1,obj.nCh);
isNoise=cell(1,obj.nCh);
for n=1:obj.nCh
    for c=1:size(avgClusteredWaveforms{n},2)
        tmpSpikes1=squeeze(avgClusteredWaveforms{n}(:,c,:));
        tmpSpikesStd1=squeeze(stdClusteredWaveforms{n}(:,c,:));
        pValid=~(tmpSpikes1>-obj.mergingNStdNoiseDetection*tmpSpikesStd1/sqrt(nAvgSpk{n}(c)-1) & tmpSpikes1<obj.mergingNStdNoiseDetection*tmpSpikesStd1/sqrt(nAvgSpk{n}(c)-1));
        pValidReduces=pValid(pSamples4NoiseDetection,:);
        validSamplesInTemplate{n}{c}=sum(pValidReduces(:))/numel(pValidReduces);
        isNoise{n}{c}=validSamplesInTemplate{n}{c}<obj.mergingNoiseThreshold;
        %{
                    z
                     if validSamplesInTemplate{n}(c)<obj.mergingNoiseThreshold
                        tVec=(1:numel(tmpSpikes1))/60;
                        plot(tVec(~pValid),tmpSpikes1(~pValid),'o','color',[0.6 0.6 0.9]);hold on;
                        plot(tVec,tmpSpikes1(:),'lineWidth',1);
            
                        ylabel('Voltage [\muV]');
                        xlabel('Time [ms]');
                        axis tight;
                        title([num2str([n c]) '- noise = ' num2str(validSamplesInTemplate{n}{c}<0.4) ', score = ' num2str(validSamplesInTemplate{n}{c})]);
                        pause;hold off;
                    end
        %}
    end
end
%hist(cell2mat(cellfun(@(x) x(:),{validSamplesInTemplate{:}},'UniformOutput',0)'),100);

lags=-maxSpikeShiftSamples:maxSpikeShiftSamples;
merge=cell(obj.nCh,obj.nCh);
fracValSamples=cell(obj.nCh,obj.nCh);

%calculate the standard deviation of the mean template (from the standard deviation of the spikes when creating the template divided by sqrt(n+1))
%and take only data points exceeding a value of N*std. For these points (M1-M2)^2/(var1+var2) is calculated and compared to the threshold.
fprintf('Calculating channels to merge');
mergingList=[];
for n1=1:obj.nCh
    for n2=obj.chPar.surChExtVec{n1}(obj.chPar.pSurChOverlap{n1})
        for c1=1:size(avgClusteredWaveforms{n1},2)
            for c2=1:size(avgClusteredWaveforms{n2},2)
                
                tmpSpikes1=squeeze(avgClusteredWaveforms{n1}(:,c1,obj.chPar.pSharedCh1{n1}{n2}));
                tmpSpikes2=squeeze(avgClusteredWaveforms{n2}(:,c2,obj.chPar.pSharedCh2{n1}{n2}));
                
                tmpSpikesStd1=squeeze(stdClusteredWaveforms{n1}(:,c1,obj.chPar.pSharedCh1{n1}{n2}));
                tmpSpikesStd2=squeeze(stdClusteredWaveforms{n2}(:,c2,obj.chPar.pSharedCh2{n1}{n2}));
                
                if obj.mergingAllignWaveShapes
                    % Compute cross-correlation
                    X=nanmean(xcorrmat(tmpSpikes1,tmpSpikes2,maxSpikeShiftSamples),2);
                    [cc,d]=max(X);
                    d=lags(d);
                    if d>0
                        tmpSpikes1=tmpSpikes1(d+1:end,:);
                        tmpSpikes2=tmpSpikes2(1:end-d,:);
                        tmpSpikesStd1=tmpSpikesStd1(d+1:end,:);
                        tmpSpikesStd2=tmpSpikesStd2(1:end-d,:);
                    else
                        tmpSpikes1=tmpSpikes1(1:end+d,:);
                        tmpSpikes2=tmpSpikes2(-d+1:end,:);
                        tmpSpikesStd1=tmpSpikesStd1(1:end+d,:);
                        tmpSpikesStd2=tmpSpikesStd2(-d+1:end,:);
                    end
                end
                
                %remove from spike shapes all noise points estimated as points below obj.mergingNStdSpikeDetection standard deviation
                pValid=~(   tmpSpikes1>-obj.mergingNStdSpikeDetection*tmpSpikesStd1/sqrt(nAvgSpk{n1}(c1)-1) & tmpSpikes1<obj.mergingNStdSpikeDetection*tmpSpikesStd1/sqrt(nAvgSpk{n1}(c1)-1) & ...
                    tmpSpikes2>-obj.mergingNStdSpikeDetection*tmpSpikesStd2/sqrt(nAvgSpk{n2}(c2)-1) & tmpSpikes2<obj.mergingNStdSpikeDetection*tmpSpikesStd2/sqrt(nAvgSpk{n2}(c2)-1)   );
                
                %tmpSpikes1(tmpSpikes1>-obj.mergingNStdSpikeDetection*tmpSpikesStd1/sqrt(nAvgSpk{n1}(c1)-1) & tmpSpikes1<obj.mergingNStdSpikeDetection*tmpSpikesStd1/sqrt(nAvgSpk{n1}(c1)-1))=NaN;
                %tmpSpikes2(tmpSpikes2>-obj.mergingNStdSpikeDetection*tmpSpikesStd2/sqrt(nAvgSpk{n2}(c2)-1) & tmpSpikes2<obj.mergingNStdSpikeDetection*tmpSpikesStd2/sqrt(nAvgSpk{n2}(c2)-1))=NaN;
                
                %tmpScore=nanmedian(((tmpSpikes1(:)-tmpSpikes2(:)).^2)./(tmpSpikesStd1(:).^2+tmpSpikesStd2(:).^2));
                %merge{n1,n2}(c1,c2)=nanmean(((tmpSpikes1(:)-tmpSpikes2(:)).^2))./sqrt(nanvar(tmpSpikes1(:))+nanvar(tmpSpikes2(:)));
                merge{n1,n2}(c1,c2)=nanmean((tmpSpikes1(pValid)-tmpSpikes2(pValid)).^2)./(nanvar(tmpSpikes1(pValid))+nanvar(tmpSpikes2(pValid)));
                fracValSamples{n1,n2}(c1,c2)=sum(pValid(:))/numel(pValid);
                
                if merge{n1,n2}(c1,c2)<=obj.mergingThreshold
                    mergingList=[mergingList; [n1 c1 n2 c2 d]];
                end
                if obj.mergingTestInitialTemplateMerging
                    %plotting tests
                    f=figure;
                    tVec=(1:numel(tmpSpikes1))/60;
                    plot(tVec(~pValid),tmpSpikes1(~pValid),'o','color',[0.6 0.6 0.9]);hold on;
                    plot(tVec(~pValid),tmpSpikes2(~pValid),'o','color',[0.9 0.6 0.6]);
                    plot(tVec,tmpSpikes1(:),'lineWidth',1);hold on;
                    plot(tVec,tmpSpikes2(:),'r','lineWidth',1);
                    
                    ylabel('Voltage [\muV]');
                    xlabel('Time [ms]');
                    axis tight;
                    title(['merge=' num2str(merge{n1,n2}(c1,c2)),', T=' num2str(obj.mergingThreshold)]);
                    %pause;hold off;
                end
            end
        end
    end
end

if obj.mergingPlotStatistics
    f=figure;
    %hist(cell2mat(cellfun(@(x) x(:),{fracValSamples{:}},'UniformOutput',0)'),100);
    hist(cell2mat(cellfun(@(x) x(:),{merge{:}},'UniformOutput',0)'),100);
    line([obj.mergingThreshold obj.mergingThreshold],ylim,'color','r');
    ylabel('# pairs');
    xlabel('merging score');
    
    print([obj.sortingDir filesep 'mergingStats'],'-djpeg','-r300');
    close(f);
end

%change the format of the average waveform to support different number of channels for different neurons on the same electrode
for i=1:1:obj.nCh
    for j=1:size(avgClusteredWaveforms{i},2)
        avgWF{i}{j}=squeeze(avgClusteredWaveforms{i}(:,j,:));
        stdWF{i}{j}=squeeze(avgClusteredWaveforms{i}(:,j,:));
        ch{i}{j}=obj.chPar.surChExtVec{i};
    end
end
%f=figure;h=axes;activityTracePhysicalSpacePlot(h,ch{i}{j},avgWF{i}{j}',obj.chPar.rEn);


%mergingList=sortrows(mergingList); %rearrange mergingList according ascending order of 1st then 2nd,3rd,4th rows
if ~isempty(mergingList)
    mergingList = unique(mergingList,'rows'); %finds the unique rows to merge
    uniqueNeurons=unique([mergingList(:,[1 2]);mergingList(:,[3 4])],'rows');
    
    %associate each neuron with a number and recalculate list
    tmpList=[mergingList(:,[1 2]);mergingList(:,[3 4])];
    numericMergeList=zeros(size(tmpList,1),1);
    
    for i=1:size(uniqueNeurons,1)
        p=find(tmpList(:,1)==uniqueNeurons(i,1) & tmpList(:,2)==uniqueNeurons(i,2));
        numericMergeList(p)=i;
    end
    numericMergeList=reshape(numericMergeList,[numel(numericMergeList)/2 2]);
    %numericDelayList = sparse([numericMergeList(:,1);numericMergeList(:,2)],[numericMergeList(:,2);numericMergeList(:,1)],[mergingList(:,5);-mergingList(:,5)]);
    
    %group connected neurons (the numeric representation corresponds to the order in uniqueNeurons
    groups=groupPairs(numericMergeList(:,1),numericMergeList(:,2));
    
    %A good alternative to recalculating templates would be to take the neuron with most spikes, this solution will not require any reloading of data,
    %however, it should be checked how reliable is this relative to the recalculation solution
    
    %Recalculate templates and merging them to one template by averaging and taking into account the potential lags
    neurons2Remove=[]; %initialization
    if obj.mergingRecalculateTemplates
        nGroups=numel(groups);
        fprintf('Recalculating templates for merged groups (/%d)',nGroups);
        
        neuron=cell(1,nGroups); %initialization
        for i=1:nGroups
            fprintf('%d ',i);
            newSpikeShapes=[]; %initialization
            minimalChMap=1:obj.nCh;
            maximalChMap=[];
            nNeuron=numel(groups{i});
            neuron{i}=zeros(nNeuron,2);  %initialization
            channel=zeros(nNeuron,1);  %initialization
            %tmpDelays=full(numericDelayList(groups{i},groups{i})); %extract the relative delays as previously calculated using the cross correlation function
            %first get the neuron numbers from groups and calculat the common channels (minimalChMap) that are recorded in all these groups
            for j=1:nNeuron
                neuron{i}(j,:)=uniqueNeurons(groups{i}(j),:);
                channel(j)=neuron{i}(j,1);
                noiseSpike(j)=isNoise{neuron{i}(j,1)}{neuron{i}(j,2)};
                [commonCh,pComN1,pComN2]=intersect(minimalChMap,obj.chPar.surChExtVec{channel(j)});
                minimalChMap=minimalChMap(pComN1);
                maximalChMap=[maximalChMap obj.chPar.surChExtVec{channel(j)}];
            end
            maxChMap=unique(maximalChMap);
            
            if isempty(minimalChMap) %some noise signal will be similar across the whole array and will come out as one group -in this case the whole array is used for comparison
                fprintf('Group %d (%d neurons) not merged due to lack of common channes',i,nNeuron);
                minimalChMap=maximalChMap;
            else
                avgWaveformsSimilarNeurons=nan(nSamples,nNeuron,numel(minimalChMap));
                for j=1:nNeuron
                    [commonCh,pComN1,pComN2]=intersect(minimalChMap,obj.chPar.surChExtVec{channel(j)});
                    avgWaveformsSimilarNeurons(:,j,pComN1)=avgClusteredWaveforms{channel(j)}(:,neuron{i}(j,2),pComN2);
                end
                [allignedWaveforms,tmpDelays]=allignSpikeShapes(permute(avgWaveformsSimilarNeurons,[2 3 1]));
                
                %load the waveform shift them according to cross correlation delays and add them to one big spike waveform matrix
                for j=1:nNeuron
                    [commonCh,pComN1,pComN2]=intersect(minimalChMap,obj.chPar.surChExtVec{channel(j)});
                    %try to rewrite this part to access only specific indices in spike shapes by defining a matfile object (maybe this can increase speed)
                    load(obj.sortingFileNames.spikeDetectionFile{channel(j)},'spikeShapes','preSpikeSamplesIntrp','minimumDetectionIntervalSamplesIntrp','detectionInt2uV');
                    load(obj.sortingFileNames.clusteringFile{channel(j)},'idx');
                    
                    pRelevantSpikes=find(idx==neuron{i}(j,2));
                    
                    %tmpSpikeShapes=detectionInt2uV.*double(spikeShapes(: , pRelevantSpikes , pComN2)); %there is a conversion 2 double here - if no shifting is needed
                    
                    tmpSpikeShapes=nan(nSamples,numel(pRelevantSpikes),numel(minimalChMap)); %there is a conversion 2 double here
                    if tmpDelays(j) >= 0
                        tmpSpikeShapes(1 : (end-tmpDelays(j)) , : , pComN1)=detectionInt2uV.*double(spikeShapes( (tmpDelays(j)+1) : end , pRelevantSpikes , pComN2)); %there is a conversion 2 double here
                    else
                        tmpSpikeShapes(-tmpDelays(j)+1:end,:,pComN1)=detectionInt2uV.*double(spikeShapes(1:end+tmpDelays(j),pRelevantSpikes,pComN2)); %there is a conversion 2 double here
                    end
                    %}
                    newSpikeShapes=cat(2,newSpikeShapes,tmpSpikeShapes);
                end
                
                %newSpikeShapes=permute(allignSpikeShapes(permute(newSpikeShapes,[2 3 1])),[3 1 2]);
                %newSpikeShapes(newSpikeShapes==0)=NaN; %to not include in the averages the padding due to spike shifting
                
                avgSpikeWaveforms=nanmedian(newSpikeShapes,2);
                stdSpikeWaveforms=1.4826*nanmedian(abs(newSpikeShapes- bsxfun(@times,avgSpikeWaveforms,ones(1,size(newSpikeShapes,2),1)) ),2);
                avgSpikeWaveforms(isnan(avgSpikeWaveforms))=0;   %average waveform should not contain NaNs
                stdSpikeWaveforms(isnan(stdSpikeWaveforms))=0;   %average waveform should not contain NaNs
                
                [~,pMin]=min(min(avgSpikeWaveforms((preSpikeSamplesIntrp-minimumDetectionIntervalSamplesIntrp):(preSpikeSamplesIntrp+minimumDetectionIntervalSamplesIntrp),:,:),[],1),[],3);
                maxChannel=commonCh(pMin); %real channel
                pMergedNeuron=find(channel==maxChannel,1,'first'); %select ch with largest amp. if there are several neurons on the same channel with the same waveform just takes the first
                if isempty(pMergedNeuron) %for the special case where the peak waveform of the joined neurons sits on a different channel than any of the original neurons
                    pMergedNeuron=1;
                    fprintf('Maximum amplitude channel for group %d, channel %d, neuron %d - detected on a different channel than any of the original neurons before merging.\n',i,neuron{i}(:,1),neuron{i}(:,2));
                end
                
                %M=newSpikeShapes;plotShifted(reshape(permute(M,[1 3 2]),[size(M,1)*size(M,3) size(M,2)]),'verticalShift',30);line([(pMin-1)*size(M,1) pMin*size(M,1)],[0 0],'color','g','lineWidth',3);
                %f=figure;h=axes;activityTracePhysicalSpacePlot(h,commonCh,squeeze(avgSpikeWaveforms(:,1,:))',obj.chPar.rEn);
                %pause;hold off;
                
                [commonCh,pComN1,pComN2]=intersect(minimalChMap,obj.chPar.surChExtVec{neuron{i}(pMergedNeuron,1)});
                avgWF{neuron{i}(pMergedNeuron,1)}{neuron{i}(pMergedNeuron,2)}=avgSpikeWaveforms;
                stdWF{neuron{i}(pMergedNeuron,1)}{neuron{i}(pMergedNeuron,2)}=stdSpikeWaveforms;
                ch{neuron{i}(pMergedNeuron,1)}{neuron{i}(pMergedNeuron,2)}=commonCh;
                isNoise{neuron{i}(pMergedNeuron,1)}{neuron{i}(pMergedNeuron,2)}=any(noiseSpike);
                
                %collect all neurons to remove from waveforms
                for j=[1:(pMergedNeuron-1) (pMergedNeuron+1):numel(groups{i})]
                    neurons2Remove=[neurons2Remove ; neuron{i}(j,:)];
                end
            end
        end
        %remove all merged waveforms
        if ~isempty(neurons2Remove)
            for i=unique(neurons2Remove(:,1))'
                p=neurons2Remove(find(neurons2Remove(:,1)==i),2);
                for j=1:numel(p)
                    avgWF{i}{p(j)}=[];
                    stdWF{i}{p(j)}=[];
                    ch{i}{p(j)}=[];
                    isNoise{i}{p(j)}=[];
                end
                if all(cellfun(@(x) isempty(x),avgWF{i})) %if after removal the channel i has no neurons this channel has to be replaced by an empty value
                    avgWF{i}=[];
                    stdWF{i}=[];
                    ch{i}=[];
                    isNoise{i}=[];
                end
            end
        end
        for i=find(cellfun(@(x) ~isempty(x),avgWF))
            pEmptyNeurons=cellfun(@(x) isempty(x),avgWF{i});
            avgWF{i}(pEmptyNeurons)=[];
            stdWF{i}(pEmptyNeurons)=[];
            ch{i}(pEmptyNeurons)=[];
            isNoise{i}(pEmptyNeurons)=[];
        end
    else
        % take the average of the largest group as the average
    end
    mergedNeurons=cellfun(@(x) [x(:,1) x(:,2)],neuron,'UniformOutput',0);
else
    mergedNeurons=[];
end
if ~obj.runWithoutSaving2File
    save(obj.sortingFileNames.mergedAvgWaveformFile,'avgWF','stdWF','ch','mergedNeurons','isNoise');
end

obj=obj.findSortingFiles; %update sorted files
