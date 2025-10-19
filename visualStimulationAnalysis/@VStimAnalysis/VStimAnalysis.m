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
        visualStimPlotsFolder %the folder with results of visual stimulation analysis plots (under visualStimFolder)
        spikeSortingFolder %the folder with spike sorting data
        visualStimulationFile %the name of the visual stimulation mat file
        stimName %the name of the visual stimulation - extracted by removing analysis from the class name
        VST %all visual stimulation properties and values
    end

    properties (SetObservable, AbortSet = true, SetAccess=public)
        dataObj %data recording object
    end

    properties (Constant, Abstract)
        trialType % The type of trials in terms of flips 'imageTrials' have one flip per trial and 'videoTrials' have many flips per trial
    end

    methods (Hidden)
        %class constructor - gets name and adds listener to update initialization every time the dataRecording object is changed
        %If
        function obj=VStimAnalysis(dataObj, params)
            arguments (Input) %ResponseWindow.mat
                dataObj
                params.Session = 1;
            end
            StimSession =  params.Session;
            obj.stimName=class(obj);obj.stimName=obj.stimName(1:end-8); %
            addlistener(obj, 'dataObj', 'PostSet',@(src,evnt)obj.initialize(StimSession));
            if nargin==0 || isempty(dataObj)
                fprintf('No dataRecording object entered as input!!! Most functions will not work!!!\n Please manually populate the dataObj property.\n');
            else
                obj.dataObj=dataObj;
            end
           
        end
    end

    methods

        function obj = initialize(obj, StimSession)
            %Initialization - extraction of folders and visual parameters subclass name should match the visual stimulation
            %extract the visual stimulation parameters from parameter file
            obj.visualStimFolder=obj.findFolderInExperiment(obj.dataObj.recordingDir,'visualStimulation');
            obj=setVisualStimulationFile(obj,'Session',StimSession);
            obj=getStimParams(obj);
            obj.spikeSortingFolder=obj.findFolderInExperiment(obj.dataObj.recordingDir,'kilosort');
        end

        function printFig(obj,f,figName)
            %Prints plots in the relevant folder
            %f - figure handle
            %figName - the name of the figure in the figure folder
            set(f,'PaperPositionMode','auto');
            disp(['Printing fig: ',obj.visualStimPlotsFolder,filesep,figName]);
            print([obj.visualStimPlotsFolder,filesep,figName],'-djpeg','-vector','-r300');
        end

        function results = getStimLFP(obj,params)
            %Extracts the LFP for all visual stimulation times from the row data.
            arguments (Input)
                obj
                params.win = [500,500] % duration [1,2] [ms] (for on and off) for LFP analysis
                params.channelSkip = 5 %includes every 5th channel
                params.getWinFromStimDuration = false %if this option is used, only the on response is calculated
                params.overwrite logical = false %if true overwrites results
                params.analysisTime = datetime('now') %extract the time at which analysis was performed
                params.inputParams = false %if true - prints out the iput parameters so that it is clear what can be manipulated in the method
            end
            if params.inputParams,disp(params),return,end

            %load previous results if analysis was previuosly performed and there is no need to overwrite otherwise continue
            results = obj.isOutputAnalysis(obj.getAnalysisFileName,params.overwrite,nargout==1);
            if ~isempty(results), return, end

            stimTimes=obj.getSyncedDiodeTriggers;

            %Design decimation Filter
            F=filterData(obj.dataObj.samplingFrequencyAP(1));
            F.downSamplingFactor=obj.dataObj.samplingFrequencyAP(1)/250;
            F=F.designDownSample;
            F.padding=true;
            samplingFreqLFP=F.filteredSamplingFrequency;

            %To add: for cases of low memory use a loop for calculating over groups of 10 trials and merge
            if params.getWinFromStimDuration
                params.win=[obj.VST.stimDuration*1000, 500];
            end

            [LFP_on,~]=obj.dataObj.getData(1:params.channelSkip:obj.dataObj.channelNumbers(end),stimTimes.stimOnFlipTimes,params.win(1));
            LFP_on=F.getFilteredData(LFP_on);

            [LFP_off,~]=obj.dataObj.getData(1:params.channelSkip:obj.dataObj.channelNumbers(end),stimTimes.stimOffFlipTimes,params.win(2));
            LFP_off=F.getFilteredData(LFP_off);

            fprintf('Saving results to file.\n');
            save(obj.getAnalysisFileName,'params','LFP_on','LFP_off','samplingFreqLFP');
        end

        function obj = getStimParams(obj)
            %extract visual stimulation parameters from the file saved by VStim classes when running visualStimGUI
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
            %extract currently running analysis method name and use it to create a unique file name for saving analysis results
            db=dbstack;currentMethod=strsplit(db(2).name,'.');
            analysisFile=[obj.visualStimAnalysisFolder,filesep,currentMethod{end},'.mat'];
            analysisFile=[obj.visualStimAnalysisFolder,filesep,currentMethod{end},'.mat'];
        end

        %check if spike sorting data exists and converts to t,ic format if needed.
        function [s]=getSpikeTIcData(obj)
            %check that sorting exist
            if isdir(obj.spikeSortingFolder)
                if isfile([obj.spikeSortingFolder filesep 'sorting_tIc.mat'])
                    s=load([obj.spikeSortingFolder filesep 'sorting_tIc.mat']);
                else
                    fprintf('Did not find <strong> sorting_tIc.mat </strong> (t,ic format) in spike sorting folder!\ntrying to convert from Phy format, please wait...\n');
                    obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
                    s=load([obj.spikeSortingFolder filesep 'sorting_tIc.mat']);
                end
            else
                fprintf('Did not find spike sorting folder, put phy results into a folder starting with <strong> kilosort </strong> and try again,\n');
            end
        end

        function [results] = getCorrSpikePattern(obj,T,trialCat,params)
            arguments
                obj
                T = []%the start times of for the trials included in the analysis and its categories
                trialCat = []% the category of each of the trials.
                params.win = 1000 %the window post stimTimes times
                params.bin = 10 %[ms] - bins size for the generated rasters
                params.gaussConvSamples = 5 % one standard deviation of the Gaussian for smoothing rasters in units of bin
                params.overwrite logical = false %if true overwrites results
                params.analysisTime = datetime('now') %extract the time at which analysis was performed
                params.inputParams = false %if true - prints out the iput parameters so that it is clear what can be manipulated in the method
            end
            if params.inputParams,disp(params),return,end

            %load previous results if analysis was previuosly performed and there is no need to overwrite otherwise continue
            results = obj.isOutputAnalysis(obj.getAnalysisFileName,params.overwrite,nargout==1);
            if ~isempty(results), return, end

            s=obj.getSpikeTIcData;

            M=BuildBurstMatrix(s.ic,round(s.t/params.bin),round(T/params.bin),round(params.win/params.bin));
            M=ConvBurstMatrix(M,fspecial('gaussian',[1 params.gaussConvSamples*3],params.gaussConvSamples),'same');
            [nTrials,nNeu,nBins]=size(M);

            [CI]=CalcCorrAlfaB_tmp2(M);

            M=reshape(M,[nTrials,nNeu*nBins])';
            C=corrcoef(M);

            fprintf('Saving results to file.\n');
            save(obj.getAnalysisFileName,'params','C','CI','trialCat');
        end

        function plotCorrSpikePattern(obj,params)
            arguments
                obj
                params.overwrite logical = true %if true overwrites the existing plots. If false only generates the plots without saving
                params.analysisTime = datetime('now') %extract the time at which analysis was performed
                params.inputParams = false %if true - prints out the iput parameters so that it is clear what can be manipulated in the method
            end
            if params.inputParams,disp(params),return,end

            res=obj.getCorrSpikePattern;

            nTrials=size(res.C,1);
            uniqCat=unique(res.trialCat,'stable');
            nCat=numel(uniqCat);
            trialsPerCat=nTrials/nCat;
            [x,y]=meshgrid(0:trialsPerCat:nTrials);

            f1=figure;
            imagesc(res.C);
            hold on;
            line(x,y,'Color','k');line(y,x,'Color','k');
            axis(f1.Children,'square');
            set(f1.Children,'XTick',trialsPerCat/2:trialsPerCat:nTrials,'YTick',trialsPerCat/2:trialsPerCat:nTrials,'XTickLabel',uniqCat,'YTickLabel',uniqCat)
            if params.overwrite,obj.printFig(f1,'trialCorrMatrix'),end

            f2=figure;
            imagesc(res.CI);
            hold on;
            line(x,y,'Color','k');line(y,x,'Color','k');
            axis(f2.Children,'square');
            set(f2.Children,'XTick',trialsPerCat/2:trialsPerCat:nTrials,'YTick',trialsPerCat/2:trialsPerCat:nTrials,'XTickLabel',uniqCat,'YTickLabel',uniqCat)
            if params.overwrite,obj.printFig(f2,'timeInvariantTrialCorrMatrix'),end

            f3=figure;
            [DC,order]=DendrogramMatrix(res.CI,'toPlotBinaryTree',1,'maxClusters',8,'figureHandle',f3);
            if params.overwrite,obj.printFig(f3,'timeInvariantDendrogramedTrialCorrMatrix'),end

            f4=figure;
            text(1:numel(order),order,res.trialCat,'FontSize',8);hold on;plot(order,'.','MarkerSize',10);
            xlabel('Trial');ylabel('Reordered trial');
            if params.overwrite,obj.printFig(f4,'timeInvariantDendrogramedOrdering'),end

            f5=figure;
            [DC,order]=DendrogramMatrix(res.C,'toPlotBinaryTree',1,'maxClusters',8,'figureHandle',f5);
            if params.overwrite,obj.printFig(f5,'dendrogramedTrialCorrMatrix'),end

            f6=figure;
            text(1:numel(order),order,res.trialCat,'FontSize',8);hold on;plot(order,'.','MarkerSize',10);
            xlabel('Trial');ylabel('Reordered trial');
            if params.overwrite,obj.printFig(f6,'dendrogramedOrdering'),end
        end

        function results = getSyncedDiodeTriggers(obj,params)
            %Sychronize the times for each stimulation onet or frame between the diode analog data and the time stamps saved by the visual stimulation class used for presenting the stimulations
            arguments
                obj
                params.minDiodeInterval = 0.5 %removes diode frames shorter than minDiodeInterval fraction of the inter frame interval
                params.analyzeOnlyOnFlips = false %in some stimuli - the off flips exist but not analyzed.
                params.ignoreNLastFlips = 0 %some stimuli have residual flips that do not need to be analyzed.
                params.overwrite logical = false %if true overwrites results
                params.analysisTime = datetime('now') %extract the time at which analysis was performed
                params.inputParams = false %if true - prints out the iput parameters so that it is clear what can be manipulated in the method
            end
            if params.inputParams,disp(params),return,end

            %load previous results if analysis was previuosly performed and there is no need to overwrite otherwise continue
            results = obj.isOutputAnalysis(obj.getAnalysisFileName,params.overwrite,nargout==1);
            if ~isempty(results), return, end

            diode=obj.getDiodeTriggers;

            allDiodeFlips=sort([diode.diodeUpCross,diode.diodeDownCross]);
            allDiodeFlips(1+find(diff(allDiodeFlips)<obj.VST.ifi*1000*params.minDiodeInterval))=[]; %remove double diode flip detections assuming intervals can not be faster than frame rate
            fprintf('%d Trigger removed due to adjuscent intervals.\n',numel(diode.diodeUpCross)+numel(diode.diodeDownCross)-numel(allDiodeFlips));
            measuredFlips=numel(allDiodeFlips);


            if isequal(obj.stimName,'StaticDriftingGrating')
                params.analyzeOnlyOnFlips =true; %Weird case where flip onset and flip offset are the same and each of them contain both flips
            end

            if isfield(obj.VST,'on_Flip')
                if ~params.analyzeOnlyOnFlips
                    allFlips=[obj.VST.on_Flip;obj.VST.off_Flip];
                else
                    allFlips=[obj.VST.on_Flip];
                end
            elseif isfield(obj.VST,'flip')
                allFlips=obj.VST.flip';
            elseif isfield(obj.VST,'flipOnsetTimeStamp')
                 if ~params.analyzeOnlyOnFlips
                    allFlips=[obj.VST.flipOnsetTimeStamp;obj.VST.flipOffsetTimeStamp];
                else
                    allFlips=[obj.VST.flipOnsetTimeStamp];
                 end
            end


            %remove flips with NaN - which are just the result of fixed matrix size
            if all(isnan(allFlips(:)))
                fprintf('<strong>All flip times in the visual stimulation meta data were NaNs!!!</strong> Please check that your vStim is valid\nTrying to use the start times of diode signals and expected frame rates\n');
                pTrialEnds=[find(diff(allDiodeFlips)>obj.VST.interTrialDelay*0.9*1000) numel(allDiodeFlips)];
                pTrialStarts=[1 1+find(diff(allDiodeFlips)>obj.VST.interTrialDelay*0.9*1000)];
                allFlips=bsxfun(@plus, (0:size(allFlips,1)-1)*obj.VST.ifi*1000,allDiodeFlips(pTrialStarts)');
            else
                allFlips=allFlips(~isnan(allFlips))*1000;
            end

            if isequal(obj.stimName,'StaticDriftingGrating')
                params.ignoreNLastFlips =2;
            end

            %ignore a few last flips - needed for some stimuli
            if params.ignoreNLastFlips~=0
                allFlips(end-params.ignoreNLastFlips+1:end)=[];
            end

            expectedFlips=numel(allFlips);
            fprintf('%d flips expected, %d found (diff=%d). Linking existing flip times with stimuli...\n',expectedFlips,measuredFlips,expectedFlips-measuredFlips);
            if (expectedFlips-measuredFlips)>0.1*expectedFlips
                fprintf('There are more than 10 percent mismatch in the number of diode and vStim expected flips. Cant continue!!! Please check diode extraction!\n');
                return;
            end
            switch obj.trialType
                case 'videoTrials'
                    pTrialEnds=[find(diff(allDiodeFlips)>obj.VST.interTrialDelay*0.9*1000) numel(allDiodeFlips)];
                    pTrialStarts=[1 1+find(diff(allDiodeFlips)>obj.VST.interTrialDelay*0.9*1000)];
                    stimOnFlipTimes=allDiodeFlips(pTrialStarts);
                    stimOffFlipTimes=allDiodeFlips(pTrialEnds);

                    if isequal(obj.stimName,'StaticDriftingGrating') %Change because static time could be equal to intertrial delat
                      
                        % Compute differences
                        dv_prev = [NaN diff(allDiodeFlips)];      % difference from previous element
                        dv_next = [diff(allDiodeFlips) NaN];      % difference from next element
                        pTrialStarts = [1 find(abs(dv_prev) >= obj.VST.static_time*0.7*1000 & abs(dv_next) >= obj.VST.interTrialDelay*0.7*1000)];
                        pTrialEnds = [find(abs(dv_prev) <= (1/obj.VST.fps)*10*1000 & abs(dv_next) >= obj.VST.interTrialDelay*0.7*1000) numel(allDiodeFlips)];

                        stimOnFlipTimes=allDiodeFlips(pTrialStarts);
                        stimOffFlipTimes=allDiodeFlips(pTrialEnds);
                    end

                    if numel(pTrialEnds)~=obj.VST.nTotTrials || numel(pTrialStarts)~=obj.VST.nTotTrials
                        disp('The total number of trials does not equal the number of inter trial delay gaps! Could not perform trial association');
                    end

                    trialDiodeFlips=cell(1,obj.VST.nTotTrials);
                    
                    try
                        diodeFrameFlipTimes=nan(size(obj.VST.flip));
                    catch
                        obj.VST.flip = obj.VST.flipOnsetTimeStamp;
                        diodeFrameFlipTimes=nan(size(obj.VST.flip));
                    end

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
                        elseif (frameMatch+sum(pDelayed))==0 && ~isequal(obj.stimName,'StaticDriftingGrating')%at least some of the frames were delayed in presentation and not just missed
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
                    results = load(obj.getAnalysisFileName);
                case 'imageTrials'
                    if expectedFlips==measuredFlips
                        if params.analyzeOnlyOnFlips
                            stimOnFlipTimes=allDiodeFlips;
                            stimOffFlipTimes=[];
                        else
                            stimOnFlipTimes=allDiodeFlips(1:2:end);
                            stimOffFlipTimes=allDiodeFlips(2:2:end);
                        end
                    else
                        error('This case of unequal triggers for image type trials was not addressed in the code!');
                        return;
                    end

                    fprintf('Saving results to file.\n');
                    save(obj.getAnalysisFileName,'params','stimOnFlipTimes','stimOffFlipTimes');
                    results=load(obj.getAnalysisFileName);
            end
        end

        function copyFilesFromRecordingFolder(obj)
            %searches visual stimulation files and copies them to a dedicated visual stimulation folder
            numberOfParentFolders=2; %How many parent folders to go in looking for the visual stimulation files

            tmpFolder=obj.dataObj.recordingDir;
            filesFound=false;
            for f=1:numberOfParentFolders
                matFiles=dir([tmpFolder,filesep,'*.mat']);
                if ~isempty(matFiles)
                    for i=1:numel(matFiles)
                        tmpVars=whos('-file',[tmpFolder filesep matFiles(i).name]);
                        if strcmp(tmpVars.name,'VSMetaData')
                            copyfile([tmpFolder filesep matFiles(i).name], obj.visualStimFolder);
                            fprintf('%s was copied to the visual stimulation folder: %s\n',matFiles(i).name,obj.visualStimFolder);
                            filesFound=true;
                        end
                    end
                end
                if ~filesFound
                    if tmpFolder(end)==filesep
                        [tmpFolder, ~, ~] = fileparts(tmpFolder(1:end-1));
                    else
                        [tmpFolder, ~, ~] = fileparts(tmpFolder(1:end));
                    end
                end
            end

            if ~filesFound
                error('No visual stimulation mat files found!!! Please copy stimulation files to the visual stimulation folder')
            end
        end

        function obj=setVisualStimulationFile(obj,params)
            %find visual stimulation file according to recording file names and the name of the visual stimulation analysis class
            arguments (Input) %ResponseWindow.mat
                obj
                params.visualStimulationfile = [];
                params.Session = 1;
            end

            if isempty(params.visualStimulationfile)
                VSFiles=dir([obj.visualStimFolder filesep '*.mat']);
                if isempty(VSFiles)
                    obj.copyFilesFromRecordingFolder;
                    VSFiles=dir([obj.visualStimFolder filesep '*.mat']);
                end

                try
                    dateTime=datetime({VSFiles.date},'InputFormat','dd-MMM-yyyy HH:mm:ss');
                catch
                    %In case pc regional setting is Israel and hebrew
                    dateTime=datetime({VSFiles.date},'InputFormat','dd-MMM-yyyy HH:mm:ss','Locale', 'he_IL');
                end
                [~,pDate]=sort(dateTime);
                VSFiles={VSFiles.name}; %do not switch with line above
                VSFiles = VSFiles(~contains(lower(VSFiles), 'metadata')); %exclude metadata
                recordingsFound=0;
                tmpDateTime = datetime.empty(0,numel(VSFiles));
                pSession = [];
                for i=1:numel(VSFiles)
                    if contains(VSFiles{i},obj.stimName,'IgnoreCase',true)
                        recordingsFound=recordingsFound+1;
                        pSession=[pSession i];
                    end
                    try
                        vStimIdentifiers=split(VSFiles{i},["_","."]);
                        tmpDateTime(i)=datetime(join(vStimIdentifiers(2:8), "-"),'InputFormat','yyyy-MM-dd-HH-mm-ss-SSS');
                    catch
                        fprintf('<strong>!!!!Important!!!!!</strong>\nUnable to extract date and time from a visual stimulation file name!!!!\nPlease correct the format or remove the file and run again!!!!\n')
                    end
                end

                if recordingsFound<1
                    fprintf('No matchings visual stimulation files found!!!\n Please check the names of visual stimulation files or run setVisualStimulationFile(file) with the filename as input.\n');
                    return;
                else
                    obj.visualStimulationFile=VSFiles{pSession(params.Session)};
                    [~,order]=sort(tmpDateTime);
                    obj.sessionOrderInRecording=find(order==pSession(params.Session));
                end
            else
                obj.visualStimulationFile=visualStimulationfile;
            end

            %populate properties and create folders for analysis if needed
            [~,fileWithoutExtension]=fileparts(obj.visualStimulationFile);
            obj.visualStimAnalysisFolder=[obj.visualStimFolder filesep fileWithoutExtension '_Analysis'];
            if ~isfolder(obj.visualStimAnalysisFolder)
                mkdir(obj.visualStimAnalysisFolder);
                fprintf('Visual stimulation Analysis folders created:\n%s\n',obj.visualStimAnalysisFolder);
            end
            obj.visualStimPlotsFolder=[obj.visualStimFolder filesep fileWithoutExtension '_Plots'];
            if ~isfolder(obj.visualStimPlotsFolder)
                mkdir(obj.visualStimPlotsFolder);
                fprintf('Visual stimulation analysis Plots folders created:\n%s\n',obj.visualStimAnalysisFolder);
            end
        end

        function results=getSessionTime(obj,params)
            %obj=getSessionTime(obj,params) - Gets the start times of each visual stimulation session from digital triggers in the recoring
            arguments (Input)
                obj
                params.startEndChannel = [] %[1,2] - The digital triger channel for stim onset and offset
                params.analysisTime = datetime('now') %extract the time at which analysis was performed
                params.inputParams = false %if true - prints out the iput parameters so that it is clear what can be manipulated in the method
                params.overwrite logical = false %if true overwrites results %if true overwrites results
            end
            if params.inputParams,disp(params),return,end

            %In future versions remove these properties and dont use them for analysis results
            %load previous results if analysis was previuosly performed and there is no need to overwrite
            if isfile(obj.getAnalysisFileName) && ~params.overwrite
                fprintf('Analysis already exists (use overwrite option to recalculate).\n');
                results=load(obj.getAnalysisFileName);
                obj.sessionStartTime=results.sessionStartTime;
                obj.sessionEndTime=results.sessionEndTime;
                obj.startSessionTrigger=results.startSessionTrigger;
                return;
            end

            %load previous results if analysis was previuosly performed and there is no need to overwrite otherwise continue
            results = obj.isOutputAnalysis(obj.getAnalysisFileName,params.overwrite,nargout==1);
            if ~isempty(results), return, end

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
                params.extractionMethod string = 'diodeThreshold' %the method used to extract frame flipes-{'diodeThreshold','digitalTriggerDiode'}
                params.analogDataCh = 1
                params.overwrite logical = false %if true overwrites results
                params.analysisTime = datetime('now') %extract the time at which analysis was performed
                params.inputParams = false %if true - prints out the iput parameters so that it is clear what can be manipulated in the method
            end
            if params.inputParams,disp(params),return,end

            %load previous results if analysis was previuosly performed and there is no need to overwrite otherwise continue
            results = obj.isOutputAnalysis(obj.getAnalysisFileName,params.overwrite,nargout==1);
            if ~isempty(results), return, end

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
                        return;
                    end
                case "digitalTriggerDiode"

                switch obj.trialType

                    case 'videoTrials'

                    if ~any(isempty([obj.sessionStartTime,obj.sessionEndTime]))

                        if all(obj.trialType == 'videoTrials') && (isequal(obj.stimName,'linearlyMovingBall')...
                                || isequal(obj.stimName, 'linearlyMovingBar'))
                            expectedFlipsperTrial = unique(obj.VST.nFrames);
                            speeds = obj.VST.speeds;
                            framesNspeed = zeros(2,length(speeds));
                            framesNspeed(1,:) =  speeds;
                            %Works for two speeds
                            framesNspeed(2,speeds == min(speeds)) = max(expectedFlipsperTrial);
                            framesNspeed(2,speeds ==  max(speeds)) = min(expectedFlipsperTrial);
                        elseif isequal(obj.stimName,'StaticDriftingGrating') || isequal(obj.stimName,'movie')
                            framesNspeed = ones(2,1);
                            framesNspeed(2,1) = round(obj.VST.actualStimDuration*obj.VST.fps); 

                            if isequal(obj.stimName,'movie')
                                framesNspeed(2,1) = round(obj.VST.movFrameCount);
                                framesNspeed = repmat(framesNspeed,1,numel(obj.VST.movieSequence));
                            end
                            
                            if isequal(obj.stimName,'StaticDriftingGrating')
                                framesNspeed(2,1) = framesNspeed(2,1)+1;%adding static time. 
                                framesNspeed = repmat(framesNspeed,1,numel(obj.VST.angleSequence));
                            end             

                        else 
                            framesNspeed = ones(2,1);
                        end

                        t = obj.dataObj.getTrigger;
                        trialOn = t{3}(t{3} > obj.sessionStartTime & t{3} < obj.sessionEndTime);
                        trialOff = t{4}(t{4} > obj.sessionStartTime & t{4} < obj.sessionEndTime);
                        interDelayMs = obj.VST.interTrialDelay*1000;

                        %[A,t_ms]=obj.dataObj.getAnalogData(params.analogDataCh,trialOn(1)-interDelayMs/2,trialOff(end)-trialOn(1)+interDelayMs); %extract diode data for entire recording
                        [A,t_ms]=obj.dataObj.getAnalogData(params.analogDataCh,trialOn(1),trialOff(end)-trialOn(1)+interDelayMs); %extract diode data for entire recording

                        DiodeCrosses = cell(2,numel(trialOn));
                        moreCross =0;
                        trialMostcross=inf;
                        intTrials =[];
                        iMC =0;
                        intrialsNum = 0;
                        trialFail =0;
                        failedTrials =[];
                        for i =1:length(trialOff)

                            startSnip  = round((trialOn(i)-trialOn(1))*(obj.dataObj.samplingFrequencyNI/1000))+1;
                            endSnip  = round((trialOff(i)-trialOn(1)+100)*(obj.dataObj.samplingFrequencyNI/1000));

                            if endSnip>length(A)
                                signal = squeeze(A(startSnip:end));
                                t_msS = t_ms(startSnip:end);
                            else
                                signal =squeeze(A(startSnip:endSnip));
                                t_msS = t_ms(startSnip:endSnip);
                            end
                            fDat=medfilt1(signal,15);
                            Th=mean(fDat(1:100:end));
                            stdS = std(fDat(1:100:end));
                            sdK = 0.1;
                            upTimes=t_msS(fDat(1:end-1)<Th-sdK*stdS & fDat(2:end)>=Th-sdK*stdS)+trialOn(1);%+interDelayMs/2; %get real recording times
                            downTimes=t_msS(fDat(1:end-1)>Th-sdK*stdS  & fDat(2:end)<=Th-sdK*stdS )+trialOn(1);%+interDelayMs/2;

                            % Filter crossings: Remove those too close together (e.g., < 50 ms)
                            minISI = 2*floor(1000/obj.VST.fps);  % ms
                            filterISI = @(x) x([true, diff(x) > minISI]);

                            try
                                DiodeCrosses{1,i} = filterISI(upTimes);
                                DiodeCrosses{2,i} = filterISI(downTimes);
                            catch
                                DiodeCrosses{1,i} = upTimes;
                                DiodeCrosses{2,i} = downTimes;
                            end

                            if (length(DiodeCrosses{1,i}) + length(DiodeCrosses{2,i}))*1.1 < framesNspeed(2,i)
                                %if the number of calculated frames is less than 10%
                                %then perform an interpolation with the
                                %first and last cross

                                if (length(DiodeCrosses{1,i}) + length(DiodeCrosses{2,i})) < framesNspeed(2,i)*0.5
                                    %%Diode failure. Use digital triggers
                                    %%and interpolate.
                                    diodeAll = zeros(1,2);
                                    diodeAll(1) = trialOn(i);
                                    diodeAll(2) = trialOff(i);
                                    ind = 2;
                                    trialFail = trialFail +1;
                                    failedTrials = [failedTrials i];
                                else
                                    [~, ind]=min([DiodeCrosses{1,i}(1)  DiodeCrosses{2,i}(1)]); %check if trial starts with up or down cross
                                    diodeAll = sort([DiodeCrosses{1,i} DiodeCrosses{2,i}]);
                                end

                                %DiodeInterp = linspace(diodeAll(1),diodeAll(end),framesNspeed(2,i));
                                DiodeInterp = diodeAll(1):1000/obj.VST.fps: min([diodeAll(1) + (framesNspeed(2,i)-1)*(1000/obj.VST.fps), trialOff(i)]);
                                if ind == 2 %Trial starts with down cross
                                    DiodeCrosses{2,i} = DiodeInterp(1:2:end);
                                    DiodeCrosses{1,i} = DiodeInterp(2:2:end);
                                else
                                    DiodeCrosses{2,i} = DiodeInterp(2:2:end);
                                    DiodeCrosses{1,i} = DiodeInterp(1:2:end);
                                end

                                intTrials = [intTrials i];
                                intrialsNum = intrialsNum+1;

                            end

                            if (length(DiodeCrosses{1,i}) + length(DiodeCrosses{2,i}))>framesNspeed(2,i)
                                %if there are more crosses than there
                                %should be
                                moreCross = moreCross+1;
                                if trialMostcross>(length(DiodeCrosses{1,i}) + length(DiodeCrosses{2,i})) - framesNspeed(2,i)
                                    trialMostcross = (length(DiodeCrosses{1,i}) + length(DiodeCrosses{2,i})) - framesNspeed(2,i);
                                    iMC = i;
                                end
                            end

                            if isequal(obj.stimName,'StaticDriftingGrating') %Make sure there is only one start frame. 
                                
                                [firstCross, idx]= min([DiodeCrosses{1,i}(1),DiodeCrosses{2,i}(1)]);

                                if firstCross > t_msS(1) + trialOn(1) + (obj.VST.static_time*1000)/2 %Add first frame that might not be read brcsuse of no diode change
                                    if idx ==1
                                        DiodeCrosses{2,i} = [50+trialOn(i) DiodeCrosses{2,i}];
                                    else
                                         DiodeCrosses{1,i} = [50+trialOn(i) DiodeCrosses{1,i}];
                                    end
                                    [firstCross, idx]= min([DiodeCrosses{1,i}(1),DiodeCrosses{2,i}(1)]);
                                end

                                ups = DiodeCrosses{1,i};
                                downs = DiodeCrosses{2,i};
                                ups = ups(ups == firstCross | ups>firstCross+obj.VST.static_time*1000*0.9);
                                downs = downs(downs == firstCross | downs>firstCross+obj.VST.static_time*1000*0.9);

                                DiodeCrosses{1,i} =  ups;
                                DiodeCrosses{2,i} = downs;
                            end


                        end

                        diodeUpCross=cell2mat(DiodeCrosses(1,:));
                        diodeDownCross=cell2mat(DiodeCrosses(2,:));

                        fprintf('%d trials out of %d have little or no diode signal, assuming diode failure but correct fliping in trials:',trialFail,length(trialOff))
                        fprintf('%.0f ', failedTrials);
                        fprintf('\n');
                        fprintf('%d trials have excess crossings out of %d; trial %d has the most excess crossings: %d',moreCross,length(trialOff),iMC,trialMostcross)
                        fprintf('\n');
                        fprintf('%d Interpolated trials (out of %d) with more than 10%% of crosses missing: ',intrialsNum,length(trialOff));
                        fprintf('%.0f ', intTrials);
                        fprintf('\n');
                        %Test
                        figure;plot(squeeze(fDat));
                        hold on;xline((DiodeCrosses{1,i} - trialOn(i))*(obj.dataObj.samplingFrequencyNI/1000))
                        xline((DiodeCrosses{2,i}  - trialOn(i))*(obj.dataObj.samplingFrequencyNI/1000),'r')
                        yline(Th)
                        hold on;xline((upTimes-trialOn(i))*(obj.dataObj.samplingFrequencyNI/1000));xline((50)*(obj.dataObj.samplingFrequencyNI/1000))
                        % %xline((trialOff(i)-trialOn(1))*(obj.dataObj.samplingFrequencyNI/1000),'b')
                        % figure;plot(1:length(DiodeCrosses{1,i}) + length(DiodeCrosses{2,i})-1,diff(sort([DiodeCrosses{1,i} DiodeCrosses{2,i}])))
                    else
                        disp('Missing start and end times!!! Please run getSessionTime before extracting triggers');
                    end

                    case 'imageTrials'

                        t = obj.dataObj.getTrigger;
                        trialOn = t{3}(t{3} > obj.sessionStartTime & t{3} < obj.sessionEndTime);
                        trialOff = t{4}(t{4} > obj.sessionStartTime & t{4} < obj.sessionEndTime);
                        interDelayMs = obj.VST.interTrialDelay*1000;

                        %[A,t_ms]=obj.dataObj.getAnalogData(params.analogDataCh,trialOn(1)-interDelayMs/2,trialOff(end)-trialOn(1)+interDelayMs); %extract diode data for entire recording
                        [A,t_ms]=obj.dataObj.getAnalogData(params.analogDataCh,trialOn(1),trialOff(end)-trialOn(1)+interDelayMs); %extract diode data for entire recording

                        DiodeCrosses = cell(2,numel(trialOn));
                        moreCross =0;
                        trialMostcross=inf;
                        intTrials =[];
                        iMC =0;
                        intrialsNum = 0;
                        trialFail =0;
                        failedTrials =[];
                        for i =1:length(trialOff)

                            startSnip  = round((trialOn(i)-trialOn(1))*(obj.dataObj.samplingFrequencyNI/1000))+1;
                            endSnip  = round((trialOff(i)-trialOn(1)+interDelayMs/2)*(obj.dataObj.samplingFrequencyNI/1000));

                            if endSnip>length(A)
                                signal = squeeze(A(startSnip:end));
                                t_msS = t_ms(startSnip:end);
                            elseif startSnip<1
                                signal =squeeze(A(1:endSnip));
                                t_msS = t_ms(1:endSnip);
                            else
                                signal =squeeze(A(startSnip:endSnip));
                                t_msS = t_ms(startSnip:endSnip);
                            end

                            fDat=-1*medfilt1(signal,(obj.VST.stimDuration/4)*1000);
                            Th=mean(fDat(1:100:end));
                            stdS = std(fDat(1:100:end));
                            sdK = 0;
                            upTimes=t_msS(fDat(1:end-1)<Th-sdK*stdS & fDat(2:end)>=Th-sdK*stdS)+trialOn(1);%+interDelayMs/2; %get real recording times
                            downTimes=t_msS(fDat(1:end-1)>Th-sdK*stdS  & fDat(2:end)<=Th-sdK*stdS )+trialOn(1);%+interDelayMs/2;

                            if length(upTimes) >1 || length(downTimes)>1
                                upTimes=upTimes(1);
                                downTimes = downTimes(2);
                            end

                            if length(upTimes) <1 || length(downTimes)<1
                                
                                    upTimes = trialOn(i)+10;
                                    downTimes = trialOff(i)+10;
                            end

                            DiodeCrosses{1,i} = upTimes;
                            DiodeCrosses{2,i} = downTimes;

                            
                        end

                        figure;plot(squeeze(fDat));
                        hold on;xline((DiodeCrosses{1,i} - trialOn(i))*(obj.dataObj.samplingFrequencyNI/1000))
                        xline((DiodeCrosses{2,i}  - trialOn(i))*(obj.dataObj.samplingFrequencyNI/1000),'r')

                        diodeUpCross=cell2mat(DiodeCrosses(1,:));
                        diodeDownCross=cell2mat(DiodeCrosses(2,:));

                end
            end
            results.diodeUpCross = diodeUpCross;
            results.diodeDownCross = diodeDownCross;

            %save results in the right file
            fprintf('Saving results to file.\n');
            save(obj.getAnalysisFileName,'params','diodeUpCross','diodeDownCross','Th');
        end

        function f=plotDiodeTriggers(obj)
            %Plots the position of detected diode triggers
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

    methods (Static)
        %find a specific folder in the experiment
        %folderLocation=findFolderInExperiment(rootFolder,folderNamePart,params)
        folderLocation=findFolderInExperiment()
        Fig = PlotZScoreComparison()

        function results=isOutputAnalysis(analysisFileName,overwrite,isOutput)
            %load previous results if analysis was previuosly performed and there is no need to overwrite
            results=[];
            if ~overwrite
                if isOutput
                    if isfile(analysisFileName)
                        fprintf('Loading saved results from file.\n');
                        results=load(analysisFileName);
                    else
                        fprintf('No results for this analyis!!! Running analysis first but will not be able to load and output the results.\n');
                    end
                else
                    if isfile(analysisFileName)
                        fprintf('Analysis already exists (use overwrite option to recalculate).\n');
                        results=false;
                    end
                end
            else
                if isOutput
                    fprintf('Cant not calculate and return results!\n <strong>Please run without output argument first and then run again to load the results.</strong>\n');
                    results=false;
                end
            end
        end

    end
end