classdef MEAAnalysis < handle
    
    properties
        par
        relevantFieldsXls
        excelRecordingDataFileName
        nTotalRecordings
        gridSorterObj
        gridSorterFolder
        parPool4Batch = false;
        
        currentDataObj
        currentRecName
        currentPRec
        currentMEAFiles
        currentMEADir
        currentAnalysisFolder
        currentPlotFolder
        currentExpFolder
        
        files
    end
    
    properties (Constant)
        xlsSheet=1;
        titleLine=1;
        startCol=1;
        figResJPG=300;
        defaultXlsFile='\\storage.laur.corp.brain.mpg.de\Data_3\Shein-IdelsonMark\DCMEA.xlsx';
        MCRackDir='MCRackData';
    end
    
    methods
        
   [data]=getSpikePositionEstimation(obj,varargin);
   [data]=getSpikeInducedFields(obj,varargin);
   [data,hand]=getCrossCorr(obj,varargin);
   [data,hand]=getSpikeParams(obj,varargin);
   
   [hand,out]=plotSpikeInducedFields(obj,varargin);
   [hand,out]=plotAvgSIFMap(obj,varargin);
   [hand,out]=plotSIFPerNeuron(obj,varargin);
       
        %% MEAAnalysis - class constructor
        function [obj]=MEAAnalysis(xlsFile)
            if nargin==0
                obj=obj.getExcelData;
            elseif nargin==1
                obj=obj.getExcelData(xlsFile);
            end
        end
        %% getFileNames
        function [obj,fileName]=getFileNames(obj,methodName)
            %get the names of mat files associated with every method (or a specific method)
            %[obj,fileName]=getFileNames(obj,methodName)
            %   methodName - the name of the method
            %   fileName - the mat file name associated with the method and a specific recording
            if nargin==1
                methodNames=methods(obj);
                handleMethods=methods('handle');
                uniqueMethods=setdiff(methodNames,handleMethods);
                for i=1:numel(uniqueMethods)
                    obj.files.(uniqueMethods{i})=[obj.currentAnalysisFolder filesep uniqueMethods{i} '.mat'];
                end
            else
                fileName=[obj.currentAnalysisFolder filesep methodName '.mat'];
                obj.files.(methodName)=[obj.currentAnalysisFolder filesep methodName '.mat'];
            end
        end
                
        %% populate gridSorter object
        function [obj,flag]=populateGridSorterObj(obj,gridSorterObj)
            if ~isempty(obj.gridSorterObj)
                disp('gridSorter object already exists');
                flag=1;
                return;
            end
            if nargin==1
                obj.gridSorterObj=gridSorter(obj.currentDataObj);
            elseif nargin==2
                obj.gridSorterObj=gridSorterObj;
            end
            [obj.gridSorterObj,flag]=obj.gridSorterObj.loadMetaData;
        end
        
        %% batchProcessData

        function [outArgAll]=batchProcessData(obj,method,fieldName,fieldValues,varargin)
        % run analysis on multiple recordings
        % method is a string with the name of the requested method. fieldName and field value are used to select recordings, the latter is a cell array
        % varargin can be used to forward input arguments to the method.            
            
            %send e-mails when running long batch analysis
            %{
                        catch errorMsg
                if obj.sendEMailMessages
                    sendMailViaGmail(obj.email4Messages,['An error occured while running runSpikeSorting on ' getenv('COMPUTERNAME') ' session ' num2str(i) '/' num2str(nExp)],errorMsg.getReport);
                end
                rethrow(errorMsg);
            end
            if obj.sendEMailMessages
                sendMailViaGmail(obj.email4Messages,['runSpikeSorting completed on ' getenv('COMPUTERNAME') ,' session success: ' num2str(nExp)]);
            end
            %}
            nOut=nargout;
            nRec=numel(fieldValues);
            
            pMultiParam=cellfun(@(x) iscell(x),varargin(2:2:end));
            %check input validity
            for i=find(pMultiParam)
                if numel(varargin{i*2})~=nRec
                    disp(['Size of cell array args for arg: ' varargin{i*2-1} ' does not match the number of recordings']);
                    return;
                end
            end
            
            %change arguments with single value to cell arrays to fit the multi value arguments
            if any(pMultiParam)
                for i=find(~pMultiParam)
                    tmpCell=cell(1,nRec);
                    tmpCell=cellfun(@(x) varargin{i*2},tmpCell,'UniformOutput',0);
                    varargin{i*2}=tmpCell;
                end
            end
            
            fprintf(['Performing batch analysis on method ' method '\nAnalyzing recording number:']);
            if obj.parPool4Batch & nRec>1
                parfor i=1:nRec
                    tmpObj=obj.setCurrentRecording([fieldName '=' fieldValues{i}]);
                    fprintf('Analyzing recording %s...\n',fieldValues{i});
                    if any(pMultiParam)
                        tmpVaragrin=cellfun(@(x) x{i},varargin(2:2:end),'UniformOutput',0);
                        newVarargin={};
                        newVarargin(1:2:numel(varargin))=varargin(1:2:end);
                        newVarargin(2:2:numel(varargin))=tmpVaragrin;
                        if nOut>0
                            [outArgs]=tmpObj.(method)(newVarargin{:});
                            outArgAll{i}=outArgs;
                        else
                            tmpObj.(method)(newVarargin{:});
                        end
                    else
                        if nOut>0
                            [outArgs]=tmpObj.(method)(varargin{:});
                            outArgAll{i}=outArgs;
                        else
                            tmpObj.(method)(varargin{:});
                        end
                    end
                    %return all non object variables
                end
            else
                for i=1:nRec
                    tmpObj=obj.setCurrentRecording([fieldName '=' fieldValues{i}]);
                    fprintf('Analyzing recording %s...\n',fieldValues{i});
                    if any(pMultiParam)
                        tmpVaragrin=cellfun(@(x) x{i},varargin(2:2:end),'UniformOutput',0);
                        newVarargin={};
                        newVarargin(1:2:numel(varargin))=varargin(1:2:end);
                        newVarargin(2:2:numel(varargin))=tmpVaragrin;
                        if nOut>0
                            [outArgs]=tmpObj.(method)(newVarargin{:});
                            outArgAll{i}=outArgs;
                        else
                            tmpObj.(method)(newVarargin{:});
                        end
                    else
                        if nOut>0
                            [outArgs]=tmpObj.(method)(varargin{:});
                            outArgAll{i}=outArgs;
                        else
                            tmpObj.(method)(varargin{:});
                        end
                    end
                    %return all non object variables
                end
            end
        end
        
        %% getVStimParams - still under development
        function obj=getVStimParams(obj,VSFile)
            error('During development!!!!!');
            if nargin==1
                VSFileLocation=dir([obj.currentMEAFiles{1}(1:end-7) '*.vs']);
                if numel(VSFileLocation)==1
                    VSFile = importdata(VSFiles.name);
                else %try to use time stamps to evaluate the identitfy of the VSfile
                    startDate=obj.currentDataObj.startDate;
                    endDate=obj.currentDataObj.startDate;
                    d=dir(obj.currentDataObj);
                end
            end
            disp(['Extracting information from: ' VSFile]);
            
        end
        
        %% getSpikeSorting
        function obj=getSpikeSorting(obj,varargin)
            %default parameters
            parseObj = inputParser;
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
            
            obj=obj.populateGridSorterObj;
            %obj.gridSorterObj.clusteringMinSpikesTotal = 2;
            %obj.gridSorterObj.clusteringMinNSpikesCluster = 1;
            
            if strcmp(obj.par.MEA_Layout,'40_Hexa')
                obj.gridSorterObj.localGridSize = 9;
                obj.gridSorterObj.localGridExt = 3;
                
            else
                obj.gridSorterObj.localGridSize = 3;
                obj.gridSorterObj.localGridExt = 1;
            end
            
            if strcmp(obj.par.MEA_Layout,'100_9x4_FlexMEA120')
                obj.gridSorterObj.filterLowPassPassCutoff=1000;
                obj.gridSorterObj.filterLowPassStopCutoff=1500;
                obj.gridSorterObj.detectionSpikeDetectionThresholdStd=7;
            end
            
            if overwrite
                obj.gridSorterObj.overwriteAll=1;
                obj.gridSorterObj.deleteSortingFiles;
            end
            
            obj.gridSorterFolder=obj.gridSorterObj.sortingDir;
            
            % Start sorting
            obj.gridSorterObj=obj.gridSorterObj.runSorting;
            
        end %runSpikeSorting
        
        function getCSDMovie(obj)
            obj.checkFileRecording;
            
            parseObj = inputParser;
            
            electrodePitch=100;
            frameRate=30;
            videoQuality=90;
            positionRealSpace=[];
            dataType='CSD';%'pot'/'CSD'
            electrodeMarker='.';
            saveData=false; %to only calculate CSD without making movie
            dataFolder=[cd filesep 'CSDProfiles'];
            makeMovie=true;
            saveFullMatrices=false;
            
            addParameter(parseObj,'avgHPWF',[],@isnumeric);
            addParameter(parseObj,'stdHPWF',[],@isnumeric);
            addParameter(parseObj,'preSpikeMs',0.5,@isnumeric);
            addParameter(parseObj,'postSpikeMs',2.5,@isnumeric);
            addParameter(parseObj,'minDist',150,@isnumeric);
            
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end
            %make parameter structure
            par=parseObj.Results;
            
            %evaluate all input parameters in workspace
            for i=1:numel(parseObj.Parameters)
                eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
            end
            
            %obj.checkFileRecording(obj.files.getSpikePositionEstimation,'Single cell position estimation file missing, please run getSpikePositionEstimation');
            
            [funName] = dbstack;funName=funName.name;funName=strsplit(funName,'.');funName=funName{end};
            saveFileName=obj.files.(funName);
            
            if exist(saveFileName,'file') & ~overwrite
                if nargout==1
                    data=load(saveFileName);
                else
                    disp('spike SNR file already exists');
                end
                return;
            end
            
            %populate grid sorter object
            obj=populateGridSorterObj(obj);
            if ~all(obj.gridSorterObj.sortingFileNames.postProcessingAnalysisExist)
                disp('Post processing data for spike sorting does not exist, please run spikePostProcessing method in grid sorter');
                return;
            else
                if isempty(avgLongWF)
                    disp('loading average spike triggered waveforms...');
                    load(obj.gridSorterObj.sortingFileNames.postProcessingAnalysisFile,'avgLongWF');
                    fprintf('done\n');
                end
                if isempty(neuronNames)
                    load(obj.gridSorterObj.sortingFileNames.postProcessingAnalysisFile,'neuronNames');
                end
            end
            
            
            
            
            %%
            [nNeu,nCh,nSamples]=size(WF);
            timeSamples=(1:nSamples)/Fs*1000-preMs;
            spikeSample=round(preMs/1000*Fs);
            
            %calculate electrode positions
            elecPos=NaN(nCh,3);
            En2=En;
            for i=1:nCh
                [n,m]=find(En2==ch(i));
                if ~isempty(n)
                    elecPos(i,:)=[m n ch(i)];
                else
                    elecPos(i,:)=ch(i);
                end
            end
            elecPos(:,1:2)=elecPos(:,1:2)*electrodePitch;
            
            if saveData
                mkdir(dataFolder);
                disp(['Data will be saved in ' dataFolder]);
            end
            
            samplingPosX=min(elecPos(:,1)):10:max(elecPos(:,1));
            samplingPosY=min(elecPos(:,2)):10:max(elecPos(:,2));
            [XuM,YuM]=meshgrid(samplingPosX,samplingPosY);
            CSD=zeros(size(WF));
            
            for i=1:nNeu
                mM=squeeze(WF(i,:,:));
                k = kcsd2d(elecPos(:,1:2), mM(elecPos(:,3),:), 'manage_data', 0, 'X', XuM, 'Y', YuM);
                
                if strcmp(dataType,'CSD')
                    dynamics=k.CSD_est;
                elseif strcmp(dataType,'pot')
                    dynamics=k.pots_est;
                else
                    error('The parameter dataType was not chosen correctly');
                end
                
                for j=1:nCh
                    [pTmpX,pTmpY]=find(XuM==elecPos(j,1) & YuM==elecPos(j,2));
                    CSD(i,j,:)=dynamics(pTmpX,pTmpY,:);
                end
                
                %[hPlot,scaleFac]=activityTracePhysicalSpacePlot([],1:120,squeeze(CSDelec(i,:,:)),En);
                %[hPlot,scaleFac]=activityTracePhysicalSpacePlot([],1:120,squeeze(WF(i,:,:)),En);
                if saveFullMatrices
                    save([dataFolder filesep 'CSD_Neuron_' num2str(neuronNames(1,i)) '-' num2str(neuronNames(2,i))],'dynamics','dataType','preMs','Fs','XuM','YuM','elecPos');
                end
                
                if makeMovie
                    mn=min(min(min(dynamics(:,:,(spikeSample+40):end),[],1),[],2),[],3);
                    mx=max(max(max(dynamics(:,:,(spikeSample+40):end),[],1),[],2),[],3);
                    l=max(abs([mn,mx]));
                    cLim=[-l l];
                    
                    
                    writerObj = VideoWriter([dataFolder filesep dataType '_neu' num2str(neuronNames(1,i)) '_' num2str(neuronNames(2,i)) '.mp4'],'MPEG-4');
                    writerObj.FrameRate=frameRate;
                    writerObj.Quality=videoQuality;
                    open(writerObj);
                    
                    F=figure('position',[50 50 550 500],'color','w');h=axes;
                    imagesc(XuM(1,:),YuM(:,1),squeeze(dynamics(:,:,i)),cLim);set(gca,'YDir','normal');hold on;
                    plot(elecPos(:,1),elecPos(:,2),'.');
                    %text(elecPos(:,1),elecPos(:,2),num2str(elecPos(:,3)));
                    cb=colorbar;
                    set(cb,'position',[0.9167    0.7600    0.0117    0.1650],'Ticks',round([-l 0 l]));
                    cb.Label.Position=[4.2 0 0];
                    ylab=ylabel(cb,dataType);
                    ylab.Position=[4.2 0 0];
                    
                    xlabel('\mum');
                    ylabel('\mum');
                    
                    axis equal tight;
                    set(h,'nextplot','replacechildren');
                    set(F,'Renderer','zbuffer');
                    
                    if isempty(positionRealSpace)
                        [~,pPeak]=min(min(min(dynamics(:,:,(spikeSample-10):(spikeSample+10)),[],1),[],2),[],3);
                        CSDSpikePeak=squeeze(mean(dynamics(:,:,(spikeSample-10+pPeak-5):(spikeSample-10+pPeak+5)),3));
                        [ySpk,xSpk] = find(CSDSpikePeak == min(min(CSDSpikePeak)));
                        cellBodyPos=[XuM(1,xSpk) YuM(ySpk,1)];
                    else
                        cellBodyPos=positionRealSpace(:,i);
                    end
                    
                    for j=1:nSamples
                        tmpImg=squeeze(dynamics(:,:,j));
                        
                        h(1)=imagesc(XuM(1,:),YuM(:,1),tmpImg,cLim);hold on;
                        h(2)=plot(elecPos(:,1),elecPos(:,2),electrodeMarker,'color',[0.8 0.8 0.8]);
                        h(3)=line([XuM(1,1) XuM(1,end)],[cellBodyPos(2) cellBodyPos(2)],'color','k');
                        h(4)=line([cellBodyPos(1) cellBodyPos(1)],[YuM(1,1) YuM(end,1)],'color','k');
                        
                        title([num2str(timeSamples(j),'% 10.1f') ' ms']);
                        
                        frame = getframe(F);
                        writeVideo(writerObj,frame);
                        delete(h);
                    end
                    close(writerObj);
                    close(F);
                end
            end
            
        end
        
        %% getHighSNRNeurons
        function [data]=getSpikeSNR(obj,varargin)
            
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'avgHPWF',[],@isnumeric);
            addParameter(parseObj,'stdHPWF',[],@isnumeric);
            addParameter(parseObj,'preSpikeMs',0.5,@isnumeric);
            addParameter(parseObj,'postSpikeMs',2.5,@isnumeric);
            addParameter(parseObj,'minDist',150,@isnumeric);
            
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end
            %make parameter structure
            par=parseObj.Results;
            
            %evaluate all input parameters in workspace
            for i=1:numel(parseObj.Parameters)
                eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
            end
            
            %obj.checkFileRecording(obj.files.getSpikePositionEstimation,'Single cell position estimation file missing, please run getSpikePositionEstimation');
            
            [funName] = dbstack;funName=funName.name;funName=strsplit(funName,'.');funName=funName{end};
            saveFileName=obj.files.(funName);
            
            if exist(saveFileName,'file') & ~overwrite
                if nargout==1
                    data=load(saveFileName);
                else
                    disp('spike SNR file already exists');
                end
                return;
            end
            
            %populate grid sorter object
            obj=populateGridSorterObj(obj);
            if ~all(obj.gridSorterObj.sortingFileNames.STWaveformExist)
                disp('Post processing data for spike sorting does not exist, please run spikePostProcessing method in grid sorter');
            else
                if isempty(avgHPWF) || isempty(stdHPWF)
                    disp('loading average+std spike triggered waveforms...');
                    load(obj.gridSorterObj.sortingFileNames.STWaveformFile,'avgHPWF','stdHPWF','neuronNames');
                    fprintf('done\n');
                end
            end
            
            [Xc,Yc]=obj.currentDataObj.getElectrodePositions;
                
            [funName] = dbstack;funName=funName.name;funName=strsplit(funName,'.');funName=funName{end};
            saveFileName=obj.files.(funName);
            
            if exist(saveFileName,'file') & ~overwrite
                if nargout==1
                    data=load(saveFileName);
                else
                    disp('spike SNR file already exists');
                end
                return;
            end

            Fs_ms=obj.currentDataObj.samplingFrequency(1)/1000;
            tSpk=(-obj.gridSorterObj.detectionPreSpikeWindow+1/Fs_ms):(1/Fs_ms):obj.gridSorterObj.detectionPostSpikeWindow;
            pSpk=tSpk>-preSpikeMs & tSpk<postSpikeMs;
            
            for i=1:size(neuronNames,2)
                spkSNRAll(i,:)=mean(abs(avgHPWF(i,:,pSpk))./stdHPWF(i,:,pSpk),3);
                [~,pMax(i)]=max(abs(spkSNRAll(i,:)));
                d=sqrt((Xc-Xc(pMax(i))).^2+(Yc-Yc(pMax(i))).^2);
                spkSNR(i)=mean(spkSNRAll(i,d<=minDist));
            end

            %save files
            save(saveFileName,'spkSNR','par');
        end
        
        %% getDigitalTriggers
        function [data]=getDigitalTriggers(obj,varargin)
            
            obj.checkFileRecording;
            
            parseObj = inputParser;
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
            
            [funName] = dbstack;funName=funName.name;funName=strsplit(funName,'.');funName=funName{end};
            saveFileName=obj.files.(funName);
            
            if exist(saveFileName,'file') & ~overwrite
                if nargout==1
                    data=load(saveFileName);
                else
                    disp('Trigger file already exists');
                end
                return;
            end
            disp(['Getting triggers for ' obj.currentRecName]);
            [tTrig]=obj.currentDataObj.getTrigger;
            
            %save files
            save(saveFileName,'tTrig');
        end
        
        function data=getDiodeSync(obj,varargin)
            data=[];
            obj.checkFileRecording(obj.files.getDigitalTriggers,'Trigger file missing, please run getDigitalTrigger');
            T=load(obj.files.getDigitalTriggers);
            
            parseObj = inputParser;
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            addParameter(parseObj,'sessionStartDigiTriggerNumber',1,@isnumeric);
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
            saveFileName=obj.files.(funName);
            
            if exist(saveFileName,'file') & ~overwrite
                if nargout==1
                    data=load(saveFileName);
                else
                    disp('Trigger file already exists');
                end
                return;
            end
            
            disp('Syncing diode signal...');
            [frameShifts,upCross,downCross,digiTriggers,transitionNotFound]=frameTimeFromDiode(obj.currentDataObj,'sessionStartDigiTriggerNumber',sessionStartDigiTriggerNumber,'T',T.tTrig);
            save(saveFileName,'par','frameShifts','upCross','downCross','digiTriggers','transitionNotFound','-v7.3');
            
        end
        
        %% getExcelData
        function [obj]=getExcelData(obj,excelRecordingDataFileName,additionalExcelFieldNames)
            if nargin==1
                obj.excelRecordingDataFileName=obj.defaultXlsFile;
            else
                obj.excelRecordingDataFileName=excelRecordingDataFileName;
            end
            
            %General rules for formatting excel files:
            %1) All title names will be extracted as fields, except if they have a '#' as a first letter
            %2) If an line does not contain data, it has to have 1 in the exclude colomn
            %3) If a field has a 
            
            
            %get data from excel spread sheet
            [excelNum,excelTxt,excelRaw] = xlsread(obj.excelRecordingDataFileName,obj.xlsSheet,'','basic'); %basic mode used for computers without installed excel
            xlsFieldNames=excelTxt(obj.titleLine,:);
            maxRow=size(excelTxt,1);
                
            %get all fields except the ones starting with #
            obj.relevantFieldsXls=excelRaw(obj.titleLine,:);
            pRelevantFields=cellfun(@(x) x(1)~='#',obj.relevantFieldsXls);
            obj.relevantFieldsXls=obj.relevantFieldsXls(pRelevantFields);
            
            %remove rows that have the exclude field true
            pExclude=find(strcmp(xlsFieldNames,'Exclude'));
            nonExludedRows=(1:(maxRow-numel(obj.titleLine)));
            if ~isempty(pExclude)
                p2Remove=cellfun(@(x) x=='1',excelRaw(obj.titleLine+nonExludedRows,pExclude+obj.startCol-1));
                nonExludedRows(p2Remove)=[];
            end
            
            for i=1:numel(obj.relevantFieldsXls)
                pFieldInExcel=find(strcmp(xlsFieldNames,obj.relevantFieldsXls{i}));
                if ~isempty(pFieldInExcel)
                    obj.par.(obj.relevantFieldsXls{i})=excelRaw(obj.titleLine+nonExludedRows,pFieldInExcel+obj.startCol-1);
                else
                    obj.par.(obj.relevantFieldsXls{i})=[];
                end
            end
            obj.nTotalRecordings=numel(obj.par.(obj.relevantFieldsXls{1}));
            
            if isunix
                for i=1:numel(obj.par.folder)
                    obj.par.folder{i}=convertPath2LinuxMPIBR(obj.par.folder{i});
                end
            end
            disp(['Experiment data retrieved from: ' num2str(obj.excelRecordingDataFileName)]);
        end
        
        %% getRecordingNames
        function getRecordingNames(obj)
            disp(obj.par.recNames);
        end
        
        %% setCurrentRecording
        function [obj]=setCurrentRecording(obj,recName) %if recNumber is negative loads all recording related files but not recording object
            %Function: select a subset of lines from the excel table for setting the current recording in the object
            %recName is a string containing the conditions (e.g. 'recNames=Ctr2,Animal=Lizard43').
            %Name subset can also be used. e.g 'recNames=Ctr' will give all recordings with recNames field starting with Ctr
            if nargin==2
                selectedFields=regexp(recName,',','split');
                for i=1:numel(selectedFields)
                    selectedValues=regexp(selectedFields{i},'=','split');
                    pRec(i,:)=cellfun(@(x) strcmp(num2str(x),selectedValues{2}),obj.par.(selectedValues{1}));
                end
                pRec=find(all(pRec==1,1));
                nRec=numel(pRec);
                
                if nRec>0
                    %collect all files in case multiple files were inserted (comma separated) or group format and couple folder with file name
                    allFullFiles=[];
                    for i=1:nRec,
                        if ~isempty(regexp(obj.par.MEAfiles{pRec(i)},'*')) %for cases in which the files are give as 'ctrl0*.mcd'
                            multipleFiles=dir([obj.par.folder{pRec(i)} filesep obj.par.MEAfiles{pRec(i)}]);
                            multipleFiles={multipleFiles.name};
                        else
                            multipleFiles=regexp(obj.par.MEAfiles{pRec(i)},',','split');
                        end
                        allFullFiles=[allFullFiles cellfun(@(x) [obj.par.folder{pRec(i)} filesep x],multipleFiles,'UniformOutput',0)];
                    end
                    obj.currentMEAFiles=allFullFiles';
                    
                    %check which data acquisition system was used
                    if strcmp(allFullFiles{1}(end-3:end),'.mcd') %MCRack recording
                        if ispc
                            obj.currentDataObj=MCRackRecording(obj.currentMEAFiles');
                        elseif isunix
                            obj.currentDataObj=MCRackRecordingNeuroshare(obj.currentMEAFiles');
                        end
                    elseif strcmp(allFullFiles{1}(end-3:end),'.rhd') %Intan recording
                        obj.currentDataObj=Intan(obj.currentMEAFiles);
                    elseif isdir(allFullFiles{1}) %OE recording
                        obj.currentDataObj=OERecording(obj.currentMEAFiles);
                    else
                        error(['dataRecording class could not be determined from the file extension: ' num2str(allFullFiles{1}(end-3:end))]);
                    end
                    obj.gridSorterObj=[]; %clear any existing grid sorter object from the past 
                    
                    %create data object
                    obj.currentRecName=recName;
                    obj.currentPRec=pRec;
                    
                    %define related folders and construct correspondin file names
                    obj.currentAnalysisFolder=[obj.par.folder{pRec(1)} filesep 'analysis' filesep recName];
                    obj.currentPlotFolder=[obj.par.folder{pRec(1)} filesep 'plots' filesep recName];
                    obj.currentExpFolder=obj.par.folder{pRec(1)};
                    obj=obj.getFileNames;
                    
                    [stat,mess,messid]=mkdir(obj.currentAnalysisFolder); %creates analysis directory if not existing
                    [stat,mess,messid]=mkdir(obj.currentPlotFolder); %creates analysis directory if not existing
                    fprintf('Current exp. set to: %s-%s @ %s\n',obj.par.MEAfiles{pRec(1)},obj.par.MEAfiles{pRec(end)},obj.par.folder{pRec(1)});
                elseif numel(pRec)==0
                    disp('Selected recording/s were not found in recording list');
                    return;
                end
                
             elseif nargin==1
                 disp('Not enough inputs, specify a recording name or number');
                 return;
             end
        end
    
        %% checkFile - check the existance of a data file and a recording object
        function isDone=checkFileRecording(obj,fileName,message)
            if isempty(obj.currentDataObj)
                error('No data recording object selected!!!!');
            end
            isDone=0;
            if nargin>1
                if ~exist(fileName,'file')
                    if nargin==2
                        if nargout==0
                            error('Relevant analysis file missing, please first run relevant function');
                        else
                            isDone=0;
                        end
                    else
                        error(message);
                    end
                else
                    isDone=1;
                end
            end
        end

    end
    
end