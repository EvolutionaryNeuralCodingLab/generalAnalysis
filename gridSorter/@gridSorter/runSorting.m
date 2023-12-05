function [obj,t,ic,avgWaveform]=runSorting(obj)
%run the grid sorter (check existance of sorting stages and runs the ones that were not yet done or that needs to be overwritten)

try
    obj=obj.findSortingFiles; %update sorted files
    
    %initiate variables
    t=[];ic=[];avgWaveform=[];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  spike Detection  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    if all(obj.sortingFileNames.spikeDetectionExist)  %check for the existence of spike shapes
        disp('Sorting will be preformed on previously detected waveforms');
    else
        obj=obj.spikeDetection;
    end
    toc;
    obj.saveMetaData;  % Save meta data if something was changed in this stage
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  feature extraction  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    if all(obj.sortingFileNames.featureExtractionExist)  %check for the existence of spike shapes
        disp('Sorting will be preformed on previously extracted features');
    else
        obj=obj.spikeFeatureExtraction;
    end
    toc;
    obj.saveMetaData;  % Save meta data if something was changed in this stage
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  clustering  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    if all(obj.sortingFileNames.clusteringExist)  %check for the existence of spike shapes
        disp('Sorting will be preformed on previously extracted clusters');
    else
        obj=obj.spikeClustering;
    end
    toc;
    obj.saveMetaData;  % Save meta data if something was changed in this stage
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Merging duplicate neurons  %%%%%%%%%%%%%%%%%%%%%
    tic;
    if obj.sortingFileNames.mergedAvgWaveformExist  %check for the existence of spike shapes
        disp('Sorting will be preformed on previously merged clusters');
    else
        obj=obj.spikeMerging;
    end
    toc;
    obj.saveMetaData;  % Save meta data if something was changed in this stage
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Fitting duplicate neurons  %%%%%%%%%%%%%%%%%%%%%
    tic;
    if obj.sortingFileNames.fittingExist  %check for the existence of spike shapes
        disp('No fitting performed!!!');
    else
        obj=obj.spikeFitting;
    end
    toc;
    obj.saveMetaData;  % Save meta data if something was changed in this stage
    
    %{
    %%%%%%%%%%%%%%%  Post-processing - general final plots and accesses sorting quality  %%%%%%%%%%%%%%%%%%%%%
    tic;
    if all(obj.sortingFileNames.STWaveformExist) || ~obj.postRunPostProcessing  %check for the existence of spike shapes
        disp('No post analysis performed!!!');
    else
        obj=obj.getSpikeTrigWF;
        obj.plotSpikeTrigWF;
    end
    toc;
    obj.saveMetaData;  % Save meta data if something was changed in this stage
    %}
    
    sendMailViaGmail('shein.mark@gmail.com',['grid sorter finished!!! ' getenv('COMPUTERNAME') 'running ' obj.dataRecordingObj.recordingDir],'');
    
catch errorMsg
    if obj.sendProgressEmail
        sendMailViaGmail(obj.sendProgressEmailTo,['grid sorter error on computer ' getenv('COMPUTERNAME') 'running ' obj.dataRecordingObj.recordingDir],errorMsg.getReport);
        rethrow(errorMsg);
    end
end
%%%%%%%%%%%  Further assessment of Cluster Quality  %%%%%%%%%%%%%%%%%%%%
%{
tic
if obj.sortingFileNames.assessQualityExist  %check for the existence of all error estimates
   load(obj.sortingFileNames.assessQualityFile,'er');
   if numel(er)==size(ic,2)
       disp('No fitting performed!!!');
   else
       [obj, er]=assessQuality(obj);
   end
else
[obj, er]=assessQuality(obj);
end
toc
%}
