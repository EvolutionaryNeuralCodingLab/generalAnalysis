function obj=findSortingFiles(obj)
%locate spike sorting related files from the same recording to determine which step were already done
%check conditions for recalculating the different stages of spike sorting. Notice that for the 3 first procedures only, calcululation of a subset of uncalculated channels is possiblee

if obj.overwriteAll
    obj.overwriteClustering=1;
    obj.overwriteFeatureExtraction=1;
    obj.overwriteFitting=1;
    obj.overwriteMerging=1;
    obj.overwritePostProcessingAnalysis=1;
    obj.overwriteSpikeExtraction=1;
    obj.overwriteQualityAssessment=1;
end

obj.sortingFileNames=[];

for i=1:obj.nCh
    obj.sortingFileNames.spikeDetectionFile{i}=[obj.sortingDir filesep 'ch_' num2str(obj.chPar.s2r(i)) '_spikeDetection.mat'];
    if ~obj.overwriteSpikeExtraction
        obj.sortingFileNames.spikeDetectionExist(i)=exist(obj.sortingFileNames.spikeDetectionFile{i},'file');
    else
        obj.sortingFileNames.spikeDetectionExist(i)=0;
    end
    
    obj.sortingFileNames.featureExtractionFile{i}=[obj.sortingDir filesep 'ch_' num2str(obj.chPar.s2r(i)) '_featureExtraction.mat'];
    if ~obj.overwriteFeatureExtraction
        obj.sortingFileNames.featureExtractionExist(i)=exist(obj.sortingFileNames.featureExtractionFile{i},'file');
    else
        obj.sortingFileNames.featureExtractionExist(i)=0;
    end
    
    obj.sortingFileNames.clusteringFile{i}=[obj.sortingDir filesep 'ch_' num2str(obj.chPar.s2r(i)) '_clustering.mat'];
    if ~obj.overwriteClustering
        obj.sortingFileNames.clusteringExist(i)=exist(obj.sortingFileNames.clusteringFile{i},'file');
    else
        obj.sortingFileNames.clusteringExist(i)=0;
    end
end

obj.sortingFileNames.avgWaveformFile=[obj.sortingDir filesep 'avgClusteredWaveforms.mat'];

obj.sortingFileNames.mergedAvgWaveformFile=[obj.sortingDir filesep 'AllMergedWaveforms.mat'];
if ~obj.overwriteMerging
    obj.sortingFileNames.mergedAvgWaveformExist=exist(obj.sortingFileNames.mergedAvgWaveformFile,'file');
else
    obj.sortingFileNames.mergedAvgWaveformExist=0;
end

obj.sortingFileNames.fittingFile=[obj.sortingDir filesep 'spikeSorting.mat'];
if ~obj.overwriteFitting
    obj.sortingFileNames.fittingExist=exist(obj.sortingFileNames.fittingFile,'file');
else
    obj.sortingFileNames.fittingExist=0;
end

obj.sortingFileNames.postProcessingAnalysisFile=[obj.sortingDir filesep 'postProcessingAnalysis.mat'];
if ~obj.overwritePostProcessingAnalysis
    obj.sortingFileNames.postProcessingAnalysisExist=exist(obj.sortingFileNames.postProcessingAnalysisFile,'file');
    if obj.sortingFileNames.postProcessingAnalysisExist %check if part of the units were already processed in a previous run
        tmp=load(obj.sortingFileNames.postProcessingAnalysisFile,'postProcessingAnalysisExist');
        if ~isempty(tmp)
            if isfield(tmp,'postProcessingAnalysisExist')
                obj.sortingFileNames.postProcessingAnalysisExist=tmp.postProcessingAnalysisExist;
            end
        else
            warning('Field "postProcessingAnalysisExist" does not exist in post processing file, maybe this sorting was done with an earilier version');
        end
    end
else
    obj.sortingFileNames.postProcessingAnalysisExist=0;
end

obj.sortingFileNames.STWaveformFile=[obj.sortingDir filesep 'STWaveform.mat'];
if ~obj.overwriteSTWaveform
    obj.sortingFileNames.STWaveformExist=exist(obj.sortingFileNames.STWaveformFile,'file');
else
    obj.sortingFileNames.STWaveformExist=0;
end

obj.sortingFileNames.assessQualityFile=[obj.sortingDir filesep 'errorEstimates.mat'];
if ~obj.overwriteQualityAssessment
    obj.sortingFileNames.assessQualityExist=exist(obj.sortingFileNames.fittingFile,'file');
else
    obj.sortingFileNames.assessQualityExist=0;
end
