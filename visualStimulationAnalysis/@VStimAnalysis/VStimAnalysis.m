classdef (Abstract) VStimAnalysis < handle

    properties
        %diodeUpCross % times [ms] of diode up crossing
        %diodeDownCross % times [ms] of diode down crossing
        startSessionTrigger = [1 2] % Trigger number for start and end of visual stimulation session
        sessionOrderInRecording = [] % the sequential position of the relevant session in case a few sessions exist in the same recording
        sessionStartTime % the start time [ms] of the stimulation session
        sessionEndTime % the end time [ms] of the stimulation session
        visualStimFolder %the folder with visual stimulation files
        visualStimAnalysisFolder %the folder with results of visual stimulation analysis (under visualStimFolder)
        visualStimulationFile %the name of the visual stimulation mat file
        stimName %the name of the visual stimulation - extracted by removing analysis from the class name
        VST %all visual stimulation properties and values
    end

    properties (SetObservable, AbortSet = true, SetAccess=public)
        dataObj %data recording object
    end

    methods (Hidden)

    end

    methods

        %class constructor - subclass name should match the visual stimulation
        function obj = initialize(obj,dataObj)
            obj.stimName=class(obj);obj.stimName=obj.stimName(1:end-8); %
            addlistener(obj, 'dataObj', 'PostSet',@(src,evnt)obj.setVisualStimFolder);
            if nargin==0
                fprintf('No recording class entered! Some functions may not work.\n Please manually populate the dataRecordingObj property');
            else
                obj.dataObj=dataObj;
            end
            obj=getStimParams(obj);
        end

        %extract visual stimulation parameters from the file saved by VStim classes when running visualStimGUI
        function obj = getStimParams(obj)
           
            VSFile=[obj.visualStimFolder filesep obj.visualStimulationFile];
            disp(['Extracting information from: ' VSFile]);
            VS=load(VSFile);

            %create structure
            if isfield(VS,'VSMetaData')
                for i=1:numel(VS.VSMetaData.allPropName)
                    obj.VST.(VS.VSMetaData.allPropName{i})=VS.VSMetaData.allPropVal{i};
                end
            else % for compatibility with old versions
                for i=1:size(VS.props,1)
                    obj.VST.(VS.props{i,1})=VS.props{i,2};
                end
            end
            if nargout>0
                VST=obj.VST;
            end

        end

        function analysisFile = getAnalysisFileName(obj)
            db=dbstack;currentMethod=strsplit(db(end).name,'.');
            analysisFile=[obj.visualStimAnalysisFolder,filesep,filesep,currentMethod{2},'.mat'];
        end

        function obj = syncDiodeTriggers(obj,params)
            arguments
                obj
                params
            end
            diode=obj.getDiodeTriggers;

            expectedFlips=numel(obj.VST.flip);
            measuredFlips=numel(diode.diodeDownCross)+numel(diode.diodeUpCross);
            fprintf('%d flips expected, %d found. Linking existing flip times with stimuli\n',expectedFlips,measuredFlips);

            obj.diodeUpCross
            obj.diodeDownCross

            obj.VST.nTotTrials
        end

        %set the visual stimulation folder as soon as a data recording object is populated
        function obj=setVisualStimFolder(obj) %the events need to be here otherwise the code does not work
            %extract the visual stimulation parameters from parameter file
            nParentFolders2Check=2;
            folderFound=false;

            %find visual stimulation folder
            tmpDir=dir([obj.dataObj.recordingDir filesep 'visualStimulation*']);
            if isempty(tmpDir) %check if not in the current data folder
                %go one folder back and look for visualStimulation folder
                fileSepTransitions=regexp(obj.dataObj.recordingDir,filesep); %look for file separation transitions
                if fileSepTransitions(end)==numel(obj.dataObj.recordingDir) %if last transition appears in the end of the folder remove this transition
                    fileSepTransitions(end)=[];
                end
                for i=1:nParentFolders2Check %repeat folder search nParentFolders2Check folders up
                    tmpCurrentFolder=obj.dataObj.recordingDir(1:fileSepTransitions(end));
                    %check parent folder for visual stimulation folder
                    tmpDir=dir([tmpCurrentFolder filesep 'visualStimulation*']);
                    if ~isempty(tmpDir)
                        VSFileLocation=[tmpCurrentFolder filesep tmpDir.name];
                        folderFound=true;
                    end
                    fileSepTransitions(end)=[];
                end
                if ~folderFound
                    % Get list of .mat files in one folder down (old
                    % location of .mat stim files)
                    [OldStimDir, ~, ~] = fileparts(obj.dataObj.recordingDir); % remove filename  
                    matFiles=dir([OldStimDir filesep '*.mat']);
                    %matFiles = dir('*.mat');

                    % Check if any .mat files exist
                    if ~isempty(matFiles)
                        % Create the new folder if it doesn't already exist
                        newFolder = 'visualStimulation';
                        if ~exist(newFolder, 'dir')
                            mkdir(newFolder);
                        end

                        % Move each .mat file into the new folder
                        for k = 1:length(matFiles)
                            oldPath = fullfile(OldStimDir, matFiles(k).name);
                            newPath = fullfile(OldStimDir, newFolder, matFiles(k).name);
                            movefile(oldPath, newPath);
                        end

                        tmpDir=dir([tmpCurrentFolder filesep 'visualStimulation*']);
                        VSFileLocation=[tmpCurrentFolder filesep tmpDir.name];         
                    else
                        error('Visual stimulation folder was not found!!! Notice the the name of the folder should be visualStimulation');
                    end
                end
            else
                VSFileLocation=[obj.dataObj.recordingDir filesep tmpDir.name];
            end

            if isfolder(VSFileLocation)
                fprintf('Visual stimulation folder found:\n%s\n',VSFileLocation);
                obj.visualStimFolder = VSFileLocation;
            else
                fprintf('Visual stimulation folder not found!!!\nPlease place visual stimulations in a folder that starts with visualStimulation and is located in the same folder as the data folder and run again');
                return;
            end

            obj=setVisualStimulationFile(obj);

        end

        function obj=setVisualStimulationFile(obj,visualStimulationfile)
            %find visual stimulation file according to recording file names and the name of the visual stimulation analysis class
            if nargin==1
                VSFiles=dir([obj.visualStimFolder filesep '*.mat']);
                dateTime=datetime({VSFiles.date},'InputFormat','dd-MMM-yyyy HH:mm:ss');
                [~,pDate]=sort(dateTime);
                VSFiles={VSFiles.name}; %do not switch with line above
                recordingsFound=0;
                for i=1:numel(VSFiles)
                    if contains(VSFiles{i},obj.stimName,'IgnoreCase',true)
                        recordingsFound=recordingsFound+1;
                        pSession=i;
                    end
                    try
                        vStimIdentifiers=split(VSFiles{i},["_","."]);
                        tmpDateTime(i)=datetime(join(vStimIdentifiers(2:8), "-"),'InputFormat','yyyy-MM-dd-HH-mm-ss-SSS');
                    catch
                        fprintf('!!!!Important!!!!!\nUnable to extract date and time from a visual stimulation file name!!!!\nPlease correct the format or remove the file and run again!!!!\n')
                    end
                end

                if recordingsFound~=1
                    fprintf('%d matchings visual stimulation files found.\n Please check the names of visual stimulation files or run setVisualStimulationFile(file) with the filename as input.\n');
                    return;
                else
                    obj.visualStimulationFile=VSFiles{pSession};
                    [~,order]=sort(tmpDateTime);
                    obj.sessionOrderInRecording=order(pSession);
                end
            else
                obj.visualStimulationFile=visualStimulationfile;
            end

            %populate properties and create folders for analysis if needed
            [~,fileWithoutExtension]=fileparts(obj.visualStimulationFile);
            obj.visualStimAnalysisFolder=[obj.visualStimFolder filesep fileWithoutExtension '_Analysis'];
            if ~isfolder(obj.visualStimAnalysisFolder)
                mkdir(obj.visualStimAnalysisFolder);
                fprintf('Visual stimulation Analysis folder created:\n%s\n',obj.visualStimAnalysisFolder);
            end
        end

        function obj=getSessionTime(obj,startEndChannel)
            if nargin == 2
                obj.startSessionTrigger=startEndChannel;
            end
            T=obj.dataObj.getTrigger;
            if isscalar(T{obj.startSessionTrigger(1)}) && isscalar(T{obj.startSessionTrigger(2)})
                obj.sessionStartTime=T{obj.startSessionTrigger(1)};
                obj.sessionEndTime=T{obj.startSessionTrigger(2)};
            else
                fprintf('There are %d session in the recording. Analysis indicated session # %d. If not please modify the sessionOrderInRecording property\n',numel(T{obj.startSessionTrigger(1)}),obj.sessionOrderInRecording);
                obj.sessionStartTime=T{obj.startSessionTrigger(1)}(obj.sessionOrderInRecording);
                obj.sessionEndTime=T{obj.startSessionTrigger(2)}(obj.sessionOrderInRecording);
            end

            sessionStartTime=obj.sessionStartTime;
            sessionEndTime=obj.sessionEndTime;
            startSessionTrigger=obj.startSessionTrigger;

            %save results in the right file
            db=dbstack;currentMethod=strsplit(db(end).name,'.');savedFileName=[obj.visualStimAnalysisFolder,filesep,currentMethod{2},'.mat'];
            save(savedFileName,'sessionStartTime','sessionEndTime','startSessionTrigger');
        end

        %Extract the frame flips from the diode signal
        function results=getDiodeTriggers(obj,params)
            arguments %(Input)
                obj
                %extractionMethod (1,1) string {mustBeMember(extractionMethod,{'diodeThreshold','digitalTriggerDiode'})} = 'diodeThreshold';
                params.extractionMethod string = 'diodeThreshold'
                params.overwrite logical = false
                params.analogDataCh = 1
            end

            if isfile(obj.getAnalysisFileName) && ~params.overwrite
                disp('Analysis already exists, use overwrite option to recalculate');
                if nargout==1
                    results=load(obj.getAnalysisFileName);
                end
                return;
            end

            fprintf('Extracting diode signal from analog channel #%d\n',params.analogDataCh);
            switch params.extractionMethod
                case "diodeThreshold"

                    if ~any(isempty([obj.sessionStartTime,obj.sessionEndTime]))
                        [A,t_ms]=obj.dataObj.getAnalogData(params.analogDataCh,obj.sessionStartTime,obj.sessionEndTime-obj.sessionStartTime); %extract diode data for entire recording

                        Th=mean(A(1:100:end));
                        diodeUpCross=t_ms(A(1:end-1)<Th & A(2:end)>=Th)+obj.sessionStartTime;
                        diodeDownCross=t_ms(A(1:end-1)>Th & A(2:end)<=Th)+obj.sessionStartTime;
                    else
                        disp('Missing start and end times!!! Please run getSessionTime before extracting triggers');
                    end

                case "digitalTriggerDiode"

            end
            %save results in the right file
            db=dbstack;currentMethod=strsplit(db(end).name,'.');savedFileName=[obj.visualStimAnalysisFolder,filesep,currentMethod{2},'.mat'];
            save(savedFileName,'params','diodeUpCross','diodeDownCross','Th');
        end

        function f=plotDiodeTriggers(obj)
            if isfile([obj.visualStimAnalysisFolder filesep 'getDiodeTriggers.mat'])
                D=load([obj.visualStimAnalysisFolder filesep 'getDiodeTriggers.mat']);
            else
                fprintf('Missing analysis: Running getDiodeTriggers is required!');return;
            end
            if isfile([obj.visualStimAnalysisFolder filesep 'getSessionTime.mat'])
                S=load([obj.visualStimAnalysisFolder filesep 'getSessionTime.mat']);
            else
                fprintf('Missing analysis: Running getSessionTime is required!');return;
            end

            if ~any(isempty([D.diodeUpCross,D.diodeDownCross,obj.sessionStartTime]))
                f=figure('Position',[100,100,1200,300]);
                h1=plot(D.diodeUpCross,ones(1,numel(D.diodeUpCross)),'^k');hold on;
                h2=plot(D.diodeDownCross,ones(1,numel(D.diodeDownCross)),'vk');
                h3=line([S.sessionStartTime;S.sessionEndTime]',[1.01;1.01]','linewidth',3);hold on;
                ylim([0.99 1.02]);
                l=legend([h1, h2, h3],{'diode up crossings','diode down crossings','stimulation session duration'});
            else
                disp('plotting triggers requires missing variables: run getDiodeTriggers, getSessionStartTime again');
            end
        end

    end
end