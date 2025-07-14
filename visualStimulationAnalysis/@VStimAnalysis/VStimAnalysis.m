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

    properties (Constant, Abstract)
        trialType % The type of trials in terms of flips 'imageTrials' have on flip per trial and 'videoTrials' have many flips per trial 
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

        function results = getStimLFP(obj,params)

            arguments (Input)
                obj
                params.win = [500,500] % duration [1,2] [ms] (for on and off) for LFP analysis
                params.channelSkip = 5 %includes every 5th channel
                params.overwrite logical = false
                params.analysisTime = datetime('now')
            end

            %load previous results if analysis was previuosly performed and there is no need to overwrite
            results=[];
            if isfile(obj.getAnalysisFileName) && ~params.overwrite
                fprintf('Analysis already exists (use overwrite option to recalculate).\n');
                if nargout==1
                    fprintf('Loading saved results from file.\n');
                    results=load(obj.getAnalysisFileName);
                end
                return;
            end

            stimTimes=obj.getSyncedDiodeTriggers;

            %Design decimation Filter
            F=filterData(vs.dataObj.samplingFrequencyAP(1));
            F.downSamplingFactor=vs.dataObj.samplingFrequencyAP(1)/250;
            F=F.designDownSample;
            F.padding=true;
            samplingFreqLFP=F.filteredSamplingFrequency;

            %To add: for cases of low memory use a loop for calculating over groups of 10 trials and merge 
            [LFP_on,t_ms]=vs.dataObj.getData(1:params.channelSkip:end,stimTimes.diodeOnFlipTimes,params.win(1));
            LFP_on=F.getFilteredData(LFP_on);

            [LFP_off,t_ms]=vs.dataObj.getData(1:params.channelSkip:end,stimTimes.diodeOffFlipTimes,params.win(2));
            LFP_off=F.getFilteredData(LFP_off);

            fprintf('Saving results to file.\n');
            save(obj.getAnalysisFileName,'params','LFP_on','LFP_off','samplingFreqLFP');
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
            db=dbstack;currentMethod=strsplit(db(2).name,'.');
            analysisFile=[obj.visualStimAnalysisFolder,filesep,currentMethod{2},'.mat'];
        end

        function results = getSyncedDiodeTriggers(obj,params)
            arguments
                obj
                params.overwrite logical = false
                params.analysisTime = datetime('now')
            end

            %load previous results if analysis was previuosly performed and there is no need to overwrite
            results=[];
            if isfile(obj.getAnalysisFileName) && ~params.overwrite
                fprintf('Analysis already exists (use overwrite option to recalculate).\n');
                if nargout==1
                    fprintf('Loading saved results from file.\n');
                    results=load(obj.getAnalysisFileName);
                end
                return;
            end
            diode=obj.getDiodeTriggers;

            allDiodeFlips=sort([diode.diodeUpCross,diode.diodeDownCross]);
            allDiodeFlips(1+find(diff(allDiodeFlips)<obj.VST.ifi*1000/2))=[]; %remove double diode flip detections assuming intervals can not be faster than frame rate
            measuredFlips=numel(diode.diodeDownCross)+numel(diode.diodeUpCross);

            if isfield(obj.VST,'on_Flip')
                allFlips=[obj.VST.on_Flip;obj.VST.off_Flip];allFlips=allFlips(:);
            elseif isfield(obj.VST,'flip')
                allFlips=obj.VST.flip';allFlips=allFlips(:)*1000;
            end
            expectedFlips=numel(allFlips);
            fprintf('%d flips expected, %d found (diff=%d). Linking existing flip times with stimuli...\n',expectedFlips,measuredFlips,expectedFlips-measuredFlips);
            if (expectedFlips-measuredFlips)>0.1*expectedFlips
                fprintf('There are more than 10% mismatch in the number of diode and vStim expected flips. Cant continue!!! Please check diode extration!');
                return;
            end
            switch obj.trialType
                case 'videoTrials'
                    pTrialEnds=[find(diff(allDiodeFlips)>obj.VST.interTrialDelay*0.9*1000) numel(allDiodeFlips)];
                    pTrialStarts=[1 1+find(diff(allDiodeFlips)>obj.VST.interTrialDelay*0.9*1000)];
                    stimOnFlipTimes=allDiodeFlips(pTrialStarts);
                    stimOffFlipTimes=allDiodeFlips(pTrialEnds);

                    if numel(pTrialEnds)~=obj.VST.nTotTrials || numel(pTrialStarts)~=obj.VST.nTotTrials
                        disp('The total number of trials does not equal the number of inter trial delay gaps! Could not perform trial association');
                    end

                    trialDiodeFlips=cell(1,obj.VST.nTotTrials);
                    diodeFrameFlipTimes=nan(size(obj.VST.flip));
                    for i=1:obj.VST.nTotTrials
                        pFlips=~isnan(obj.VST.flip(i,:)); %not all trials have the same number of flips
                        currentDiodeFlipTimes=allDiodeFlips(pTrialStarts(i):pTrialEnds(i));
                        currentPCFlipTimes=obj.VST.flip(i,pFlips)*1000;
                        %plot(allFlips-allFlips(1),ones(1,numel(allFlips)),'or');hold on;plot(allDiodeFlips-allDiodeFlips(1),ones(1,numel(allDiodeFlips)),'.k');
                        %check for cases in which there was a missed frame in diode signal - this could be from a missed frame or a delayed frame
                        pDelayed=diff(currentDiodeFlipTimes)>obj.VST.ifi*1000*1.5;
                        frameMatch=numel(currentPCFlipTimes)-(numel(currentDiodeFlipTimes)+sum(pDelayed));
                        if frameMatch>=0 %all frames that were not present were missed
                            [in,p]=ismembertol(currentPCFlipTimes-currentPCFlipTimes(1),currentDiodeFlipTimes-currentDiodeFlipTimes(1),obj.VST.ifi*1000/2,'DataScale',1);
                        elseif (frameMatch+sum(pDelayed))==0 %at least some of the frames were delayed in presentation and not just missed
                            %look for a delay in diode flips which may explain a consistent delay
                            tmpDiode=currentDiodeFlipTimes+cumsum([zeros(1,-frameMatch) pDelayed])*(-obj.VST.ifi*1000);
                            [in,p]=ismembertol(currentPCFlipTimes-currentPCFlipTimes(1),tmpDiode-tmpDiode(1),obj.VST.ifi*1000/2,'DataScale',1);
                        elseif frameMatch<0 && (numel(currentPCFlipTimes)-numel(currentDiodeFlipTimes))>=0 %identifies irregularities in timing
                            in=ones(1,numel(currentDiodeFlipTimes));
                            p=1:numel(currentDiodeFlipTimes);
                        elseif frameMatch<0 %the match is negative and there was no compensation due the previous condition
                            currentDiodeFlipTimes=currentDiodeFlipTimes(1:numel(currentPCFlipTimes)); %remove excess frames at end and match one to one
                            in=ones(1,numel(currentPCFlipTimes));
                            p=1:numel(currentPCFlipTimes);
                        else
                            error('The case which there are more diode flips than stimulated frames was not addressed in the algorithm! Please add to code.')
                        end
                        diodeFrameFlipTimes(i,in)=currentDiodeFlipTimes(p(p~=0));

                        %{
                    plot(currentPCFlipTimes-currentPCFlipTimes(1),ones(1,numel(currentPCFlipTimes)),'.k');hold on;
                    plot(currentDiodeFlipTimes-currentDiodeFlipTimes(1),ones(1,numel(currentDiodeFlipTimes)),'or');legend({'trig','diode'})
                        %}
                    end
                    fprintf('Saving results to file.\n');
                    save(obj.getAnalysisFileName,'params','diodeFrameFlipTimes','stimOnFlipTimes','stimOffFlipTimes');
                case 'imageTrials'
                    if expectedFlips==measuredFlips
                        stimOnFlipTimes=allDiodeFlips(1:2:end);
                        stimOffFlipTimes=allDiodeFlips(2:2:end);
                    else
                        error('This case of unequal triggers for image type trials was not addressed in the code!');
                        return;
                    end
                    fprintf('Saving results to file.\n');
                    save(obj.getAnalysisFileName,'params','stimOnFlipTimes','stimOffFlipTimes');
            end
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
                    error('Visual stimulation folder was not found!!! Notice the the name of the folder should be visualStimulation');
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

            if isfolder(VSFileLocation)
                fprintf('Visual stimulation folder found:\n%s\n',VSFileLocation);
                obj.visualStimFolder = VSFileLocation;
            else
                fprintf('Visual stimulation folder not found!!!\nPlease place visual stimulations in a folder that starts with visualStimulation and is located in the same folder as the data folder and run again');
                return;
            end

            obj=setVisualStimulationFile(obj);
            %simon's code
            %{

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
                        cd(OldStimDir)
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
            %}

        end

        function obj=setVisualStimulationFile(obj,visualStimulationfile)
            %find visual stimulation file according to recording file names and the name of the visual stimulation analysis class
            if nargin==1
                VSFiles=dir([obj.visualStimFolder filesep '*.mat']);
                try
                    dateTime=datetime({VSFiles.date},'InputFormat','dd-MMM-yyyy HH:mm:ss');
                catch
                    %In case pc regional setting is Israel and hebrew
                    dateTime=datetime({VSFiles.date},'InputFormat','dd-MMM-yyyy HH:mm:ss','Locale', 'he_IL');
                end
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
                    fprintf('No matchings visual stimulation files found!!!\n Please check the names of visual stimulation files or run setVisualStimulationFile(file) with the filename as input.\n');
                    return;
                else
                    obj.visualStimulationFile=VSFiles{pSession};
                    [~,order]=sort(tmpDateTime);
                    obj.sessionOrderInRecording=find(order==pSession);
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

        function obj=getSessionTime(obj,params)
            arguments (Input)
                obj
                params.startEndChannel = []
                params.overwrite logical = false
                params.analysisTime = datetime('now')
            end


            %load previous results if analysis was previuosly performed and there is no need to overwrite
            if isfile(obj.getAnalysisFileName) && ~params.overwrite
                fprintf('Analysis already exists (use overwrite option to recalculate).\n');
                results=load(obj.getAnalysisFileName);
                obj.sessionStartTime=results.sessionStartTime;
                obj.sessionEndTime=results.sessionEndTime;
                obj.startSessionTrigger=results.startSessionTrigger;
                return;
            end

            if nargin == 2
                obj.startSessionTrigger=params.startEndChannel;
            end
            T=obj.dataObj.getTrigger;
            if isscalar(T{obj.startSessionTrigger(1)}) && isscalar(T{obj.startSessionTrigger(2)}) %only 1 visual stimulation session in the recording
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
            fprintf('Saving results to file.\n');
            save(obj.getAnalysisFileName,'sessionStartTime','sessionEndTime','startSessionTrigger');
        end

        %Extract the frame flips from the diode signal
        function results=getDiodeTriggers(obj,params)
            arguments (Input)
                obj
                %extractionMethod (1,1) string {mustBeMember(extractionMethod,{'diodeThreshold','digitalTriggerDiode'})} = 'diodeThreshold';
                params.extractionMethod string = 'diodeThreshold'
                params.overwrite logical = false
                params.analogDataCh = 1
                params.analysisTime = datetime('now')
            end

            %load previous results if analysis was previuosly performed and there is no need to overwrite
            results=[];
            if isfile(obj.getAnalysisFileName) && ~params.overwrite
                fprintf('Analysis already exists (use overwrite option to recalculate).\n');
                if nargout==1
                    fprintf('Loading saved results from file.\n');
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
                     if ~any(isempty([obj.sessionStartTime,obj.sessionEndTime]))
                        [A,t_ms]=obj.dataObj.getAnalogData(params.analogDataCh,obj.sessionStartTime,obj.sessionEndTime-obj.sessionStartTime); %extract diode data for entire recording               
                        
                        
                        
                        diodeUpCross=;
                        diodeDownCross=;
                    else
                        disp('Missing start and end times!!! Please run getSessionTime before extracting triggers');
                    end
            end
            results.diodeUpCross = diodeUpCross;
            results.diodeDownCross = diodeDownCross;

            %save results in the right file
            fprintf('Saving results to file.\n');
            save(obj.getAnalysisFileName,'params','diodeUpCross','diodeDownCross','Th');
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