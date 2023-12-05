function [clusterMerged,Merge]=SpikeTempDiffMerging(spikeShapes,Clusters,Templates,crit)
%merge clusters in a group that have to be merge based on the residuals between spikes and templates.
%synthax : [clusterMerged,Merge]=SpikeTempDiffMerging(spikeShapes,Clusters,Templates,crit)
%input:
%           - spikeShapes : all the spikeShapes Nspikes-sizeSpike-NbChannels
%           - Clusters : vector containing the index of spikes in spikeShapes
%           - Templates : The templates corresponding to the clusters
%           nChannels-sizeTemplate-nTemplates
%           - crit : the value of the threshold to merge clusters(default value : 1.5)
%output :
%           - clusterMerged:vector of symbolic merging
%           - Merge : binary symmetric merging decision matrix
%           nclust-nclust.
if nargin<4
    crit=1.5;
end
numClust=size(Templates,3);
clusterMerged=1:numClust; % a group is assigned to every cluster
%build ajacent matrix of clusters
Merge=zeros(numClust,numClust);
for i=1:numClust
    for j=1:numClust
        if j~=i
            clSel=find(Clusters==i);
            cTmp1=Templates(:,:,i)';
            cTmp2=Templates(:,:,j)';
            diffSpikeTemp1=zeros(1,length(clSel));diffSpikeTemp2=zeros(1,length(clSel));
            for k=1:length(clSel)
                cSpike=reshape(spikeShapes(clSel(k),:,:),size(spikeShapes,2),size(spikeShapes,3));
                diffSpikeTemp1(k)=sum((cTmp1(:)-cSpike(:)).^2);
                diffSpikeTemp2(k)=sum((cTmp2(:)-cSpike(:)).^2);
            end
            Merge(i,j)=median(diffSpikeTemp2)/median(diffSpikeTemp1);
        end
    end
end

%Symetrize matrix
for k=1:numClust
    for m=k+1:numClust
        if Merge(k,m)<crit &&Merge(m,k)<crit
            Merge(k,m)=1;Merge(m,k)=1;
            clusterMerged(clusterMerged==k)=clusterMerged(m);
        else
            Merge(k,m)=0;Merge(m,k)=0;
        end
    end
end