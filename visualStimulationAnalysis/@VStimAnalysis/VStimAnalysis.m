classdef (Abstract) VStimAnalysis < handle

    properties
        diodeUpCross % times [ms] of diode up crossing
        diodeDownCross % times [ms] of diode down crossing
        startSessionTrigger = [1 2] % Trigger number for start and end of visual stimulation session
        sessionNumInRecording = [] % the number of the relevant session in case a few sessions exist in the same recording
        sessionStartTime % the start time [ms] of the stimulation session
        sessionEndTime % the end time [ms] of the stimulation session
        visualStimulationFolder
        visualStimulationFile
        stimName %the name of the visual stimulation - extracted by removing analysis from the class name
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
        end

        %extract visual stimulation parameters from the file saved by VStim classes when running visualStimGUI
        function obj = getStimParams(obj)
           
            % Future implementaion - if finding name strategy does not work try to find the correct file according to its save date which should correspond to recording dates
            %find(dateNumber>obj.currentDataObj.startDate & dateNumber<obj.currentDataObj.endDate);

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

        function isAnalysis = checkAnalysisExistance(obj)
            db=dbstack;
            savedFileName=[obj.visualStimulationFolder,filesep,'Analysis',filesep,db(end).name,'.mat'];
            isAnalysis=isfile(savedFileName)
        end

        function saveMethodAnalysis(obj,varargin)
            db=dbstack;
            savedFileName=[obj.visualStimulationFolder,filesep,'Analysis',filesep,db(end).name,'.mat'];
            save(savedFileName,varargin{2:end});
        end

        function obj = syncDiodeTriggers(obj)

        end

        %set the visual stimulation folder as soon as a data recording object is populated
        function setVisualStimFolder(obj) %the events need to be here otherwise the code does not work
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
                    error('Visual stimulation folder was not found!!! Notice the the name of the folder should be visualStimulation');
                end
            else
                VSFileLocation=[obj.dataObj.recordingDir filesep tmpDir.name];
            end
            fprintf('Visual stimulation folder set as:\n %s\n',VSFileLocation);
            %find visual stimulation file according to recording file name
            VSFiles=dir([VSFileLocation filesep '*.mat']);
            dateTime=datetime({VSFiles.date},'InputFormat','dd-MMM-yyyy HH:mm:ss');
            VSFiles={VSFiles.name}; %do not switch with line above
        
            for i=1:numel(VSFiles)
                if contains(VSFiles{i},class(obj),'IgnoreCase',true)
                    return;
                end
            end
            
            obj.visualStimulationFile=VSFiles{i};

            obj.visualStimulationFolder = [obj.dataObj.recordingDir filesep 'visualStimlation'];
            if ~isfolder(obj.visualStimulationFolder)
                mkdir(obj.visualStimulationFolder);
            end
        end

        function obj=getSessionTime(obj,startEndChannel)
            if nargin == 2
                startSessionTrigger=startEndChannel;
            end
            T=obj.dataObj.getTrigger;
            if numel(T{obj.startSessionTrigger(1)})==1 && numel(T{obj.startSessionTrigger(2)})==1
                obj.sessionStartTime=T{obj.startSessionTrigger(1)};
                obj.sessionEndTime=T{obj.startSessionTrigger(2)};
            else
                fprintf('There are %d session in the recording. Assuming the session # %d. If not please modify the sessionNumInRecording property\n',numel(T{obj.startSessionTrigger(1)}),obj.sessionNumInRecording);
                obj.sessionStartTime=T{obj.startSessionTrigger(1)}(obj.sessionNumInRecording);
                obj.sessionEndTime=T{obj.startSessionTrigger(2)}(obj.sessionNumInRecording);
            end
            
        end

        %Extract the frame flips from the diode signal
        function [obj]=getDiodeTriggers(obj,dataObj,analogDataCh,NameValueArgs)
            arguments (Input)
                obj
                dataObj
                analogDataCh
                %extractionMethod (1,1) string {mustBeMember(extractionMethod,{'diodeThreshold','digitalTriggerDiode'})} = 'diodeThreshold';
                NameValueArgs.extractionMethod string = 'diodeThreshold'
                NameValueArgs.overwrite logical = False
            end

            if checkAnalysisExistance(obj,overwrite) && ~overwrite
                disp('Analysis already exists, use overwrite option to recalculate');
                return;
            end

            switch extractionMethod
                case "diodeThreshold"

                    if ~any(isempty([obj.sessionStartTime,obj.sessionEndTime]))
                        [A,t_ms]=dataObj.getAnalogData(analogDataCh,obj.sessionStartTime,obj.sessionEndTime-obj.sessionStartTime); %extract diode data for entire recording

                        Th=mean(A(1:100:end));
                        obj.diodeUpCross=t_ms(A(1:end-1)<Th & A(2:end)>=Th)+obj.sessionStartTime;
                        obj.diodeDownCross=t_ms(A(1:end-1)>Th & A(2:end)<=Th)+obj.sessionStartTime;
                    else
                        disp('Missing start and end times!!! Please run getSessionTime before extracting triggers');
                    end

                case "digitalTriggerDiode"

            end
            saveMethodAnalysis(obj,analogDataCh,extractionMethod);

        end


        function f=plotDiodeTriggers(obj)
            if ~any(isempty([obj.diodeUpCross,obj.diodeDownCross,obj.sessionStartTime]))
                f=figure('Position',[100,100,1200,300]);
                h1=plot(obj.diodeUpCross,ones(1,numel(obj.diodeUpCross)),'^k');hold on;
                h2=plot(obj.diodeDownCross,ones(1,numel(obj.diodeDownCross)),'vk');
                h3=line([obj.sessionStartTime;obj.sessionEndTime]',[1.01;1.01]','linewidth',3);hold on;
                ylim([0.99 1.02]);
                l=legend([h1, h2, h3],{'diode up crossings','diode down crossings','stimulation session duration'});
            else
                disp('plotting triggers requires running the methods: getDiodeTriggers, getSessionStartTime');
            end
        end

    end
end