function mergeUnits()
% Notice!!!! In this process waveforms are not merged. To do so, we need to go over all spikes are recalculate

%merge all waveforms
clustersInCh=cellfun(@(x) size(x,2),avgWaveform);
pClustersInActiveCh=find(clustersInCh>0);
clustersInActiveCh=clustersInCh(pClustersInActiveCh);
activeChannels=selectedChannels(pClustersInActiveCh);

%calculate indices in merged waveform matrix
endIdxClustersInActiveCh=cumsum(clustersInActiveCh);
startIdxClustersInActiveCh=[1 endIdxClustersInActiveCh(1:end-1)+1];
%build merged waveform matrix
nWaveforms=sum(clustersInActiveCh);
avgWaveformAll=zeros(size(avgWaveform{pClustersInActiveCh(1)},1),nWaveforms,size(avgWaveform{pClustersInActiveCh(1)},3));
for i=1:numel(clustersInActiveCh)
    avgWaveformAll(:,startIdxClustersInActiveCh(i):endIdxClustersInActiveCh(i),:)=avgWaveform{pClustersInActiveCh(i)};
end
%examine distances between all waveforms in merged matrix
waveformSimilarity=nan(nWaveforms);
for i=1:nWaveforms
    for j=(i+1):nWaveforms
        waveformSimilarity(i,j)=mean(mean(bsxfun(@minus,avgWaveformAll(:,i,:),avgWaveformAll(:,j,:)).^2,1),3);
    end
end
[pMerge1 pMerge2]=find(waveformSimilarity<4);

%merge all neurons below threshold
for i=1:numel(pMerge1)
    if pMerge1(i)~=pMerge2(i)
        %find the
        mergeCh1(i)=find(pMerge1(i)>=startIdxClustersInActiveCh & pMerge1(i)<=endIdxClustersInActiveCh);
        mergeCh2(i)=find(pMerge2(i)>=startIdxClustersInActiveCh & pMerge2(i)<=endIdxClustersInActiveCh);
        
        %real electrode name for merged neuron
        mergeCh1Name(i)=activeChannels(mergeCh1(i));
        mergeCh2Name(i)=activeChannels(mergeCh2(i));
        %real neuron number in electrode for merged neuron
        mergeNeu1(i)=pMerge1(i)-startIdxClustersInActiveCh(mergeCh1(i))+1;
        mergeNeu2(i)=pMerge2(i)-startIdxClustersInActiveCh(mergeCh2(i))+1;
        %merge neurons
        [t,ic]=MergeNeurons(t,ic,[mergeCh1Name(i);mergeNeu1(i)],[mergeCh2Name(i);mergeNeu2(i)]);
        %update the merging list since some of the neurons were already merged previously
        
        deleteChPlace=pMerge2(i);
        pMerge2(pMerge2==deleteChPlace)=pMerge1(i);
        pMerge1(pMerge1==deleteChPlace)=pMerge1(i);
    end
end