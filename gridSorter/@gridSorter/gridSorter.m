classdef gridSorter
    properties (Constant)
        gridSorterVersion=1.03; %the version of grid sorter
    end
    properties (SetAccess=public)
        overwriteAll=0; %if true overwrites all sorting process
        overwriteSpikeExtraction=0; %if true recalculate spike extraction
        overwriteFeatureExtraction=0; %if true recalculate  feature extraction
        overwriteClustering=0; %if true recalculate clustering
        overwriteMerging=0; %if true recalculate cluster merging
        overwriteFitting=0; %if true recalculate spike fitting
        overwritePostProcessingAnalysis=0; %if true recalculate post processign analysis
        overwriteSTWaveform=0; %if true recalculate spike triggered analysis
        overwriteQualityAssessment=0; %if true recalculate quality assessment
                
        sortingDir %the root directory of the sorted data
        runWithoutSaving2File=false; %run without saving files to hard drive (applicable only for short recordings in which all data can be kept in memory)
        fastPrinting=0; %uses fast printing of the figures (the resulting figures are low quality and require the screen to be turned on)
        saveFigures=1; %save the figures plotted during sorting to file
        
        selectedChannelSubset=[]; %the numbers of a subset of channels to be sorted out of the whole recording
        dataRecordingObj %the data recording object (must be a subclass of dataRecording, e.g. MCRackDataNeuroshare, Intan)
        sendProgressEmail = true; %whether to send an email at the end or during errors
        sendProgressEmailTo = 'shein.mark@gmail.com';
        
        upSamplingFrequencySpike=60000; %Desired upsampling frequency for spike extrema and shape calculation. Must be an integer multiple of the sampling frequency (if not it is adjusted automatically)
        
        localGridSize=3; % The length of the local grid on which sorting is calculated, e.g. if==3, sorts on a 3x3 grid surrounding the central channel (must be an odd number). If empty uses all channels
        localGridExt=1; %the overhead used by the algorithm to later merge neurons that were detected on more than one channel (according to peak amplitude).
        
        detectionRemoveAllElectrodeMedian=false; %spike detection - remove artifacts that appear on all electrodes using a median filter substruction
        detectionMaxSpikeAmp=1000; %Spike detection - spike with higher max amplitude are considered as noise and are not analyzed
        detectionNQuantizationBits=16; %Spike detection - the quantization depth for reducing the raw data resolution for reducing disk space (should be the same as the quantization of the recording hardware)
        detectionGaussianityWindow=20;%Spike detection [ms] - window for estimating noise segments according to kurtosis
        detectionPreSpikeWindow=2;%Spike detection [ms] - window before spike peak for spike shape extraction
        detectionPostSpikeWindow=3;%Spike detection [ms] - window after spike peak for spike shape extraction
        detectionSpikeTimeShiftInterval=1;%Spike detection [ms] - initial extension of the spike window used to compensate for aligning spikes to the peak
        detectionPeakDetectionSmoothingWindow=0.3; %Spike detection [ms] - the smoothing window for detection of the exact spike occurance
        detectionKurtosisNoiseThreshold=3;%Spike detection - the threshold on the kurtosis value that differentiates noise samples from spike samples
        detectionSpikeDetectionThresholdStd=5;%Spike detection - number of standard deviations above the noise level for spike detection
        detectionMaxChunkSize=2*60*1000; %Spike detection [ms] - the size of the consecutive chunks of raw data extracted (should fit memory capacity)
        detectionChunkOverlap=1; %Spike detection - overlap between analyzed chunks to not miss a spike that fell between data chunks
        detectionMinimumDetectionInterval=0.5;%Spike detection [ms] - the maximal interval on which the global spike minima is calculated for all surrounding local extended channels (spike with minima differnt than the center channel are not included)
        detectionRemoveSpikesNotExtremalOnLocalGrid=true; %Spike detection - removes all spike that have a stronger minimum on another channel surrounding the detection channel on the grid
        detectionRemoveSpikesSubMinimumDelays=true; %Spike detection - remove spikes with delays shorter than detectionMinimumDelayBetweenSpikes on the same channel
        detectionMinimumDelayBetweenSpikes=1; %Spike detection - minimum delay allowed between spike peaks on the same channel
        
        simpleButter=1; %use simple butterworth or a more complex
        
        highPassCutoff=200;
        lowPassCutoff=2500;
        filterOrder=2;
        filterDesign='butter'; %Filter - type of filter used usually 'butter' or 'elliptic'

        filterHighPassPassCutoff=250; %Filter - highpass cutoff - upper bound
        filterHighPassStopCutoff=150; %Filter - highpass cutoff - lower bound
        filterLowPassPassCutoff=2000; %Filter - lowpass cutoff - upper bound
        filterLowPassStopCutoff=2500; %Filter - highpass cutoff - lower bound
        filterAttenuationInHighpass=20; %Filter - attenuation of power [dB] for highpass between lower and upper
        filterAttenuationInLowpass=10; %Filter - attenuation of power [dB] for lowpass between lower and upper
        filterRippleInPassband=0.5; %Filter - allowed ripples [dB] in pass band
        
        featuresMaxSpikesToCluster=10000; %Feature extraction - maximum number of spikes for feature extraction (not to template matching)
        featuresNWaveletCoeff=30; %Feature extraction - number of wavelet coefficient to extract
        featuresReduceDimensionsWithPCA=true; %Feature extraction - if to reduce dimensions of extracted features with PCA
        featuresDimensionReductionPCA=6; %Feature extraction - final dimension for clustering after wavelet decomposition
        featuresSelectedWavelet='haar'; %Feature extraction - the mother wavelet used for wavelet decomposition
        featuresWTdecompositionLevel=4; %Feature extraction - the level of wavelet decomposition
        featuresFeatureExtractionMethod='wavelet'; %Feature extraction - the method for extracting features 'wavelet'/'PCA'
        featuresConcatenateElectrodes=true; %Feature extraction - concatenate all surrounding electrodes to 1 big trace before extracting features
        
        clusteringMethod='meanShift';%Clustering - the clustering method used: 'meanShift'/'kMeans'
        clusteringMergingMethod='projectionMeanStd';%Clustering - 'projectionMeanStd';%'MSEdistance'
        clusteringMaxIter=1000; %Clustering - max itteration for clustering algorithm
        clusteringNReplicates=50; %Clustering - number of clustering replicates in k-means clustering algorithm
        clusteringInitialClusterCentersMethod='cluster'; %Clustering - method for initial conditions in k-means algorithm
        clusteringMaxPointsInSilhoutte=1000; %Clustering - the maximal number of data points per cluster used for silloute quality estimation
        clusteringMergeThreshold=0.18; %Clustering - threshold for merging clusters with 'projectionMeanStd'  method, the higher the threshold the less cluster will merge
        clusteringSTDMergeFac=2; %Clustering - the number of standard deviations for checking the crossing point on both sides (the higher the more clusters separate into groups)
        clusteringRunSecondMerging=0; %Clustering - run merging again after initial merging
        clusteringMaxClusters=12; %Clustering - the maximum number of clusters for a specific channels
        clusteringMinimumChannelRate=0.01; %Clustering [Hz] - the minimal spike rate for spikes in a given clusters
        clusteringMinNSpikesCluster=10; %Clustering - the minimal number of spikes in a cluster
        clusteringMinSpikesTotal=20; %Clustering - the minimal total number of spikes per channel (if there was a lower spike count in the channel, channel is disgarded
        clusteringPlotProjection=1; %Clustering - if to plot the projection data for 'projectionMeanStd' merging method
        clusteringPlotClassification=1;%Clustering - if to plot the final clustering result for the channel
        
        mergingThreshold=0.1; %Merging - the threshold for merging similar waveforms
        mergingRecalculateTemplates=true;%Merging - recalculate templates and not use the ones provided by the clustering algorithm
        mergingNStdSpikeDetection=4;%Merging - standard deviation threshold for remove samples from the spike waveform for template comparison
        mergingAllignSpikeBeforeAveraging=true;%Merging - allign spike before comparing templates
        mergingAllignWaveShapes=true;%Merging -
        mergingTestInitialTemplateMerging=false;%Merging -
        mergingPlotStatistics=true;%Merging - if to plot merging statistics
        mergingNStdNoiseDetection=6;%Merging -
        mergingPreSpike4NoiseDetection=0.2;%Merging -
        mergingPostSpike4NoiseDetection=0.2;%Merging -
        mergingNoiseThreshold=0.2;%Merging -
        mergingMaxSpikeShift=1; %%Merging [ms] -
        
        fittingTemplateMethod=1;
        fittingMaxLag=0.5;
        fittingLagIntervalSamples=1;
        fittingSpikeNoiseEdgeInterval=0.5;
        fittingMaxMinkowskiDist=4; %threshold for the distance of spike from template (distant templates are rejected
        
        postMaxSpikes2Present=600; %Post processing - max number of spike to extract per neuron
        postPlotAllAvgSpikeTemplates=true; %Post processing - plot small templates for all neurons
        postExtractFilteredWaveformsFromSpikeTimes=true;%Post processing - extract filtered short spike shapes with template
        postExtractRawLongWaveformsFromSpikeTimes=true;%Post processing - extract raw long LFP shapes with template
        postFilteredSNRStartEnd=[-0.2 0.2];%Post processing [ms] - [pre,post] spike interval for calculating spikeSNR
        postRawSNRStartEnd=[5 20];%Post processing [ms] - [pre,post] LFP interval for calculating :LFP SNR
        postPlotFilteredWaveforms=true;%Post processing - if to plot the filtered spike waveforms
        postPlotRawLongWaveforms=true;%Post processing - if to plot the raw post spike field waveforms
        postPlotSpikeReliability=true;%Post processing - if to plot the spike reliability plot
        postPreFilteredWindow=2;%Post processing - pre spike window for plotting
        postTotalFilteredWindow=5;%Post processing - post spike window for plotting
        postPreRawWindow=20;%Post processing - pre raw window for plotting
        postTotalRawWindow=100;%Post processing - post raw window for plotting
        postRunPostProcessing=true;%if to run post processing
        postAlternativeFileName;%if to save post processing results in a different file
        
        assessNumDiscrimDim=10;
        assessMinSpikeNum=10;
        assessMaxSpikeNumber=1000;
                
        sortingFileNames %a structure with all the files names that grid sorter calculates
    end
    
    properties (SetAccess=protected)
        nCh %number of channels to be sorted
        chPar %a structure with channel arrangement parameters that are internally used by gridSorter
        arrayExt %The extension of the array
        
        detectionInt2uV %Detection - translation between uint16 to microVolts (data is save in uint16 format to save space)
        filterObj %Filter -the highpass filter object used to filter the raw data
    end
    
    properties (SetObservable, AbortSet = true, SetAccess=public)
    end
    
    properties (Hidden, SetAccess=protected)
    end
    
    methods (Hidden)
        %class constractor
        function obj=gridSorter(dataRecordingObj,varargin)
            %addlistener(obj,'visualFieldBackgroundLuminance','PostSet',@obj.initializeBackground); %add a listener to visualFieldBackgroundLuminance, after its changed its size is updated in the changedDataEvent method
            
            %Collects all options - if properties are given as a 'propertyName',propertyValue series
            for i=1:2:length(varargin)
                eval(['obj.' varargin{i} '=' 'varargin{i+1};'])
            end
            
            if nargin==0
                disp('Data recording object not enter. In later versions a GUI will be imlemented for recording selection');
                return;
            else
                if isobject(dataRecordingObj) %if input is a data recording object
                    obj.dataRecordingObj=dataRecordingObj;
                    if iscell(obj.dataRecordingObj.recordingName)
                        [~, name, ~] = fileparts(obj.dataRecordingObj.recordingName{1});
                    else
                        [~, name, ~] = fileparts(obj.dataRecordingObj.recordingName);
                    end
                    if ~iscell(obj.dataRecordingObj.recordingDir)
                        obj.sortingDir=[obj.dataRecordingObj.recordingDir filesep name '_spikeSort'];
                    else
                        obj.sortingDir=[obj.dataRecordingObj.recordingDir{1} '_spikeSort'];
                    end
                    initiateDataRecObj=0;
                elseif exist(dataRecordingObj)==2 %if input is the name of the metadata file
                    metaData=load(dataRecordingObj);
                    sortingDir=metaData.props((strcmp(metaData.props(:,1),'sortingDir')),2);
                    obj=obj.loadMetaData(sortingDir{1});
                    initiateDataRecObj=1;
                elseif exist(dataRecordingObj)==7
                    obj=obj.loadMetaData(dataRecordingObj);
                    initiateDataRecObj=1;
                else
                    disp('First input must be a data recording object of a gridSorter metaData file');
                    return;
                end
            end
            %verify that the recording dir in the recording object is the same as for grid sorter in case directory was replaced
            %if not it changes the dataObj folder assuming that the spikeSorting folder is always a subfolder of the data folder
            if ~strcmp(obj.dataRecordingObj.recordingDir,obj.sortingDir(1:numel(obj.dataRecordingObj.recordingDir)))
                fprintf('\nPlease notice that the data folder of the recording object is different than the one of gridSorter\nReplacing data folder in dataRecording object!!!\n'); 
                foldersDataRecObj = regexp(obj.dataRecordingObj.recordingDir,filesep,'split');foldersDataRecObj(strcmp(foldersDataRecObj,''))=[];
                foldersGSObj = regexp(obj.sortingDir,filesep,'split');foldersGSObj(strcmp(foldersGSObj,''))=[];
                nFolder=numel(foldersGSObj);
                while nFolder>0
                    k = find(strcmp(foldersGSObj(nFolder),foldersDataRecObj));
                    if ~isempty(k)
                        p1=regexp(obj.dataRecordingObj.recordingDir,foldersGSObj(nFolder));
                        p2=regexp(obj.sortingDir,foldersGSObj(nFolder));
                        obj.dataRecordingObj.recordingDir=[obj.sortingDir(1:p2{1}-1) obj.dataRecordingObj.recordingDir(p1{1}:end)];
                        nFolder=0;
                    else
                        nFolder=nFolder-1;
                    end
                end
                %reinitiate data recording obj
            end
            
            if initiateDataRecObj
                fullRecFiles=cellfun(@(x) [obj.dataRecordingObj.recordingDir x],obj.dataRecordingObj.dataFileNames,'UniformOutput',0);
                disp('Reinitiating data recording object');
                eval(['obj.dataRecordingObj=' class(obj.dataRecordingObj) '(fullRecFiles);']);
            end
            
            obj=obj.calculateChParameters;
            if isempty(obj.nCh)
                error('No channels in recording or no layout provided');
            end
            
            obj=obj.findSortingFiles;
            
            %make directory in recording folder with spike sorting data
            if ~exist(obj.sortingDir,'dir')
                mkdir(obj.sortingDir);
                disp(['Creating spike sorting folder: ' obj.sortingDir]);
            end
            
        end
        
    end
    
    methods
        
        function deleteSortingFiles(obj)
            d=dir([obj.sortingDir filesep 'ch_*spikeDetection.mat']);
            if numel(d)>0
                d=dir(obj.sortingDir);
                tmpFiles=cellfun(@(x) [obj.sortingDir filesep x],{d.name},'UniformOutput',0);
                delete(tmpFiles{:});
            else
                fprintf('Please recheck that sorting data exists in this folder: %s\n',obj.sortingDir);
            end
        end
        
        function data=getSortedUnit(obj,neuronNames)
            tmp=load(obj.sortingFileNames.fittingFile);
            nNeu=size(neuronNames,2);
            data=cell(1,nNeu);
            for i=1:nNeu
                pNeu=find(tmp.ic(1,:)==neuronNames(1,i) & tmp.ic(2,:)==neuronNames(2,i));
                data{i}=tmp.t(tmp.ic(3,pNeu):tmp.ic(4,pNeu));
            end
        end
        
        function data=getSortedData(obj,vars)
            if nargin==1
                data=load(obj.sortingFileNames.fittingFile);
            else
                if iscell(vars)
                    data=load(obj.sortingFileNames.fittingFile,vars{:});
                else
                    data=load(obj.sortingFileNames.fittingFile,vars);
                end
            end
        end
        
        function data=getPostProcessedData(obj,vars)
            if nargin==1
                data=load(obj.sortingFileNames.postProcessingAnalysisFile);
            else
                if iscell(vars)
                    data=load(obj.sortingFileNames.postProcessingAnalysisFile,vars{:});
                else
                    data=load(obj.sortingFileNames.postProcessingAnalysisFile,vars);
                end
            end
        end
        
        function data=getSTData(obj,vars)
            if nargin==1
                data=load(obj.sortingFileNames.STWaveformFile);
            else
                if iscell(vars)
                    data=load(obj.sortingFileNames.STWaveformFile,vars{:});
                else
                    data=load(obj.sortingFileNames.STWaveformFile,vars);
                end
            end
        end
        
        [obj,t,ic,avgWaveform]=runSorting(obj)
        
        obj=spikeDetection(obj)
        
        gridObj=detectedSpikes2TIC(gridObj,varargin);
        
        [obj,spikeFeaturesAll]=spikeFeatureExtraction(obj)
        
        [obj,idx,initIdx,nClusters,avgSpikeWaveforms,stdSpikeWaveforms]=spikeClustering(obj)
        
        [obj,avgWF,stdWF,ch,mergedNeurons]=spikeMerging(obj)
        
        [obj,t,ic]=spikeFitting(obj)
        
        obj=spikePostProcessing(obj,tStartEnd)
        
        obj=calculateChParameters(obj)
        
        obj=findSortingFiles(obj)
        
        obj=setSortingFilesMainDir(obj)
        
        obj=getSpikeTrigWF(obj,startEnd);
        
        plotSpikeTrigWF(obj,startEnd,printFolder);
        
        [publicProps]=getProperties(obj)
        
        obj=manuallySelectValidUnits(obj)
        
        [obj, er]=assessQuality(obj)
        
        [obj,metaDataExistFlag] = loadMetaData(obj,dataRecordingFolder)
        
        saveMetaData(obj)
    end %methods
    
    methods (Static)
        
        [clusterMerged,Merge]=SpikeTempDiffMerging(spikeShapes,Clusters,Templates,crit)
        
        [gc,f]=projectionMerge(spikeFeatures,initIdx,varargin)
        
    end %methods (static)
    
end
