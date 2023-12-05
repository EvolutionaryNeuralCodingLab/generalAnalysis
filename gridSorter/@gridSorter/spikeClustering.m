function [obj,idx,initIdx,nClusters,avgSpikeWaveforms,stdSpikeWaveforms]=spikeClustering(obj)
% Run clustering on previously extracted features
% [obj,idx,initIdx,nClusters,avgSpikeWaveforms,stdSpikeWaveforms]=spikeClustering(obj)
avgClusteredWaveforms=cell(1,obj.nCh);
stdClusteredWaveforms=cell(1,obj.nCh);

if isempty(obj.sortingDir)
    obj.runWithoutSaving2File=true;
    idxAll=cell(1,obj.nCh);
    initIdxAll=cell(1,obj.nCh);
    nClustersAll=cell(1,obj.nCh);
else
    obj.runWithoutSaving2File=false;
end

fprintf('\nClustering on channel (total %d): ',obj.nCh);
for i=find(obj.sortingFileNames.clusteringExist==0 | obj.overwriteClustering) %go over all channels in the recording
    fprintf('%d ',i);
    MaxClustersTmp=obj.clusteringMaxClusters;
    
    if ~exist(obj.sortingFileNames.featureExtractionFile{i},'file')
        warning(['No feature extraction file was found for Channel ' num2str(i) '. Clustering not performed!']);
        continue;
    else
        load(obj.sortingFileNames.featureExtractionFile{i});
        load(obj.sortingFileNames.spikeDetectionFile{i},'spikeTimes');
    end
    
    [nSpikes,nFeatures]=size(spikeFeatures);
    
    if nSpikes >= obj.clusteringMinSpikesTotal && nSpikes >= (spikeTimes(end)-spikeTimes(1))/1000*obj.clusteringMinimumChannelRate
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%  Clustering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        switch obj.clusteringMethod
            case 'kMeans'
                % set options to k-Means
                opts = statset('obj.clusteringMaxIter',obj.clusteringMaxIter);
                try
                    [initIdx] = kmeans(spikeFeatures,obj.clusteringMaxClusters,'options',opts,...
                        'emptyaction','singleton','distance','city','onlinephase','on','obj.clusteringNReplicates',obj.clusteringNReplicates,'start',obj.clusteringInitialClusterCentersMethod);
                catch %if the number of samples is too low, kmeans gives an error -> try kmeans with a lower number of clusters
                    MaxClustersTmp=round(obj.clusteringMaxClusters/2);
                    [initIdx] = kmeans(spikeFeatures,MaxClustersTmp,'options',opts,...
                        'emptyaction','singleton','distance','city','onlinephase','on','obj.clusteringNReplicates',obj.clusteringNReplicates,'start',obj.clusteringInitialClusterCentersMethod);
                end
            case 'meanShift'
                initIdx=zeros(nSpikes,1);
                out=MSAMSClustering(spikeFeatures');
                for j=1:numel(out)
                    initIdx(out{j})=j;
                end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%  Merging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        switch obj.clusteringMergingMethod
            case 'MSEdistance'
                %calculate templates
                load(obj.sortingFileNames.spikeDetectionFile{i},'spikeShapes');
                avgSpikeWaveforms=zeros(nSpikeSamples,max(1,MaxClustersTmp),nSurroundingChannels);
                for j=1:MaxClustersTmp
                    avgSpikeWaveforms(:,j,:)=median(spikes4Clustering(:,initIdx==j,:),2);
                end
                
                [gc,Merge]=obj.SpikeTempDiffMerging(permute(spikes4Clustering,[2 1 3]),initIdx,permute(avgSpikeWaveforms,[3 1 2]));
                
                if obj.saveFigures
                    f1=figure('Position',[50 50 1400 900]);
                    set(f1,'PaperPositionMode','auto');
                    if obj.fastPrinting
                        imwrite(frame2im(getframe(f1)),[obj.sortingDir filesep 'Ch_' num2str(obj.chPar.s2r(i)) 'projectionTest.jpeg'],'Quality',90);
                    else
                        print([obj.sortingDir filesep 'Ch_' num2str(obj.chPar.s2r(i)) 'projectionTest'],'-djpeg','-r300');
                    end
                    close(f1);
                end
                
            case 'projectionMeanStd'
                [gc,f1]=obj.projectionMerge(spikeFeatures,initIdx,'obj.clusteringMinNSpikesCluster',obj.clusteringMinNSpikesCluster,'obj.clusteringSTDMergeFac',obj.clusteringSTDMergeFac,'obj.clusteringMergeThreshold',obj.clusteringMergeThreshold,'obj.clusteringPlotProjection',obj.clusteringPlotProjection);
                
                if obj.saveFigures && ~isempty(f1);
                    set(f1,'PaperPositionMode','auto');
                    if obj.fastPrinting
                        imwrite(frame2im(getframe(f1)),[obj.sortingDir filesep 'Ch_' num2str(obj.chPar.s2r(i)) 'projectionTest.jpeg'],'Quality',90);
                    else
                        print([obj.sortingDir filesep 'Ch_' num2str(obj.chPar.s2r(i)) 'projectionTest'],'-djpeg','-r300');
                    end
                    close(f1);
                end
                
                if obj.clusteringRunSecondMerging
                    uniqueClusters=unique(gc);
                    nClusters=numel(uniqueClusters);
                    idx=zeros(nSpikes4Clustering,1);
                    for k=1:nClusters
                        p=find(gc==uniqueClusters(k));
                        for j=1:numel(p)
                            idx(initIdx==p(j))=k;
                        end
                    end
                    MaxClustersTmp=nClusters;
                    initIdx=idx;
                    
                    [gc,f1]=obj.projectionMerge(spikeFeatures,initIdx,'obj.clusteringMinNSpikesCluster',obj.clusteringMinNSpikesCluster,'obj.clusteringSTDMergeFac',obj.clusteringSTDMergeFac,'obj.clusteringMergeThreshold',obj.clusteringMergeThreshold,'obj.clusteringPlotProjection',obj.clusteringPlotProjection);
                    
                    if obj.saveFigures
                        set(f1,'PaperPositionMode','auto');
                        if obj.fastPrinting
                            imwrite(frame2im(getframe(f1)),[obj.sortingDir '\Ch_' num2str(obj.chPar.s2r(i)) 'projectionTest.jpeg'],'Quality',90);
                        else
                            print([obj.sortingDir '\Ch_' num2str(obj.chPar.s2r(i)) 'projectionTest'],'-djpeg','-r300');
                        end
                        close(f1);
                    end
                end
        end
        
        %reclassify clusters
        uniqueClusters=unique(gc);
        nClusters=numel(uniqueClusters);
        idx=zeros(nSpikes,1);
        for k=1:nClusters
            p=find(gc==uniqueClusters(k));
            for j=1:numel(p)
                idx(initIdx==p(j))=k;
            end
        end
        
        %calculate spikeShape statistics
        load(obj.sortingFileNames.spikeDetectionFile{i},'spikeShapes','detectionInt2uV');
        spikeShapes=double(spikeShapes) .* detectionInt2uV;
        
        [nSpikeSamples,nSpikes,nSurroundingChannels]=size(spikeShapes);
        if nSpikes>obj.featuresMaxSpikesToCluster
            spikeShapes=spikeShapes(:,1:obj.featuresMaxSpikesToCluster,:);
        end
        
        avgSpikeWaveforms=zeros(nSpikeSamples,max(1,nClusters),nSurroundingChannels);
        stdSpikeWaveforms=zeros(nSpikeSamples,max(1,nClusters),nSurroundingChannels);
        for j=1:nClusters
            pCluster=idx==j;
            avgSpikeWaveforms(:,j,:)=median(spikeShapes(:,pCluster,:),2);
            stdSpikeWaveforms(:,j,:)=1.4826*median(abs(spikeShapes(:,pCluster,:)- bsxfun(@times,avgSpikeWaveforms(:,j,:),ones(1,sum(pCluster),1)) ),2);
            nSpk(j)=numel(pCluster);
        end
    else
        fprintf('X '); %to note than no neurons were detected on this electrode
        
        idx=ones(nSpikes,1);
        initIdx=ones(nSpikes,1);
        avgSpikeWaveforms=[];
        stdSpikeWaveforms=[];
        nClusters=0;
        nSpk=0;
    end
    
    avgClusteredWaveforms{i}=avgSpikeWaveforms;
    stdClusteredWaveforms{i}=stdSpikeWaveforms;
    nAvgSpk{i}=nSpk;
    
    if ~obj.runWithoutSaving2File
        save(obj.sortingFileNames.clusteringFile{i},'idx','initIdx','nClusters','avgSpikeWaveforms','stdSpikeWaveforms');
    else
        idxAll{i}=idx;
        initIdxAll{i}=initIdx;
        nClustersAll{i}=nClusters;
    end
    
    if obj.clusteringPlotClassification && nClusters>0
        
        cmap=lines;
        
        f2=figure('Position',[100 100 1200 800],'color','w');
        PCAfeaturesSpikeShapePlot(spikeFeatures,spikeShapes,obj.upSamplingFrequencySpike,initIdx,idx,obj.chPar.En,obj.chPar.s2r(obj.chPar.surChExtVec{i}),'hFigure',f2,'cmap',cmap);
        
        f3=figure('Position',[100 100 1200 800],'color','w');
        featureSubSpacePlot(spikeFeatures,idx,'hFigure',f3,'cmap',cmap);
        
        if obj.saveFigures
            figure(f2);
            set(f2,'PaperPositionMode','auto');
            if obj.fastPrinting
                imwrite(frame2im(getframe(f2)),[obj.sortingDir '\Ch_' num2str(obj.chPar.s2r(i)) 'classification.jpeg'],'Quality',90);
            else
                print([obj.sortingDir '\Ch_' num2str(obj.chPar.s2r(i)) 'classification'],'-djpeg','-r300');
            end
            
            figure(f3);
            set(f3,'PaperPositionMode','auto');
            if obj.fastPrinting
                imwrite(frame2im(getframe(f3)),[obj.sortingDir '\Ch_' num2str(obj.chPar.s2r(i)) 'featureSpace.jpeg'],'Quality',90);
            else
                print([obj.sortingDir '\Ch_' num2str(obj.chPar.s2r(i)) 'featureSpace'],'-djpeg','-r300');
            end
            if ishandle(f2)
                close(f2);
            end
            if ishandle(f3)
                close(f3);
            end
        end
    end %plot initial classification classification
end %go over all channels

if ~obj.runWithoutSaving2File
    save(obj.sortingFileNames.avgWaveformFile,'avgClusteredWaveforms','stdClusteredWaveforms','nAvgSpk');
end

obj=obj.findSortingFiles; %update sorted files
