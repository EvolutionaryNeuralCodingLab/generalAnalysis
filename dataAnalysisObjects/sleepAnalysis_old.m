classdef sleepAnalysis < handle
    properties
        par
        parFreqBandDetection
        parDBRatio
        excelRecordingDataFileName
        fieldNames
        nSessions
        filt
        currentDataObj
        currentPRec
        currentAnalysisFolder
        currentPlotFolder
        currentExpFolder
        currentVideosFolder
        files
        parPool4Batch=false;
        relevantFields={'Sleep','Awake','Eye','HR','Num','Animal','CheetahFolder','VideoFiles','MatroxTrigScheme','FrameRange','VideoType','AnimalState','DVRLFPCh','ElecLayout','Sorting','tStartAwake'};
        relevantFieldsXls
        nTotalRecordings
        currentMEAFiles
        gridSorterObj
        currentRecName
    end
    properties (Constant)
        xlsSheet=1;
        titleLine=1;
        startCol=1;
        figResJPG=400;
        figResEPS=1000; %is not important
        defaultDir='\\storage.laur.corp.brain.mpg.de\Data_2\chronicLizardExpp'
        defaultXlsFile='\\storage.laur.corp.brain.mpg.de\Data_2\chronicLizardExp\SleepExp.xlsx';
    end
    methods
        
        %% getDayTimeInRecTime
        function data=getSleepVsLights(obj,varargin)
            %% parameter and settings
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'referenceClock','19:00:00'); %reference for lights on/off
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'clockStartTime',[]); %cell array with the format 'HH:MM:SS'
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
            
            %make parameter structure
            parDayTimeOnRecTime=parseObj.Results;
            
            %check if analysis was already done done
            obj.files.dayTimeOnRecTime=[obj.currentAnalysisFolder filesep 'dayTimeOnRecTime.mat'];
            if exist(obj.files.dayTimeOnRecTime,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.dayTimeOnRecTime);
                else
                    disp('dayTimeOnRecTime file already exists');
                end
                return;
            end
            
            dbRatioFile=[obj.currentAnalysisFolder filesep 'dbRatio_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(dbRatioFile,'Delta to beta analysis missing, please first run getDBRatio');
            dataDB=load(dbRatioFile,'t_ms'); %load data 
            
            dbAutocorrFile=[obj.currentAnalysisFolder filesep 'dbAutocorr_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(dbAutocorrFile,'Delta to beta autocorr analysis missing, please first run getDBRatioAC');
            dataAC=load(dbAutocorrFile,'pSleepDBRatio','period'); %load data

            if strcmp(obj.currentDataObj.startDate(1),'(m/d/y):')
                obj.currentDataObj=obj.currentDataObj.getStartRecordingTime;
            end
            recordingStartTimeClock=obj.currentDataObj.startDate;
            
            pStartSleep=find(dataAC.pSleepDBRatio==1,1,'first');
            pEndSleep=find(dataAC.pSleepDBRatio(pStartSleep:end)==0,1,'first')+pStartSleep;

            if pStartSleep==1
                sleepStartEnd=[0 dataDB.t_ms(pEndSleep)];
            else
                sleepStartEnd=dataDB.t_ms([pStartSleep pEndSleep]);
            end
            
            if exist([obj.currentAnalysisFolder filesep 'light.mat'],'file')
                l=load([obj.currentAnalysisFolder filesep 'light.mat']);
                startSleepFromRef_h=(sleepStartEnd-l.light(1))/1000/60/60;
                manualLightAnnotation=true;
            else
                manualLightAnnotation=false;
                tmpDV=datevec(datenum(referenceClock,'HH:MM:SS')-datenum(recordingStartTimeClock,'HH:MM:SS') );
                if tmpDV(1)<0
                    startSleepFromRef_h=sleepStartEnd/1000/60/60+(24-tmpDV(:,4)+(60-tmpDV(:,5))/60+(60-tmpDV(:,6))/3600);
                    disp('Interval between start recording and reference time was too large -> assuming recording started after reference time');
                else
                    startSleepFromRef_h=sleepStartEnd/1000/60/60-(tmpDV(:,4)+tmpDV(:,5)/60+tmpDV(:,6)/3600);
                end
            end
            
            save(obj.files.dayTimeOnRecTime,'sleepStartEnd','recordingStartTimeClock','referenceClock','startSleepFromRef_h','parDayTimeOnRecTime','manualLightAnnotation');
        end
        
        %% getSpikeSTAs
        function data=getSpikeSTAs(obj,varargin)
            %% parameter and settings
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'nCycles',10,@isnumeric); 
            addParameter(parseObj,'cycleSelection','first',@(x) any(strcmp(x,{'first','rand'})));
            addParameter(parseObj,'binSW',10,@isnumeric);
            addParameter(parseObj,'preSW',1000,@isnumeric);
            addParameter(parseObj,'winSW',2000,@isnumeric);
            addParameter(parseObj,'binSO',1000,@isnumeric);
            addParameter(parseObj,'preSO',40000,@isnumeric);
            addParameter(parseObj,'winSO',80000,@isnumeric);
            addParameter(parseObj,'minSpikeRate',0.05,@isnumeric);
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
            
            %make parameter structure
            parSpikeSTA=parseObj.Results;
            
            %check if analysis was already done done
            obj.files.spikeSTA=[obj.currentAnalysisFolder filesep 'spikeSTA.mat'];
            if exist(obj.files.spikeSTA,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.spikeSTA);
                else
                    disp('Spike STA file already exists');
                end
                return;
            end
            
            slowCyclesFile=[obj.currentAnalysisFolder filesep 'slowCycles_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(slowCyclesFile,'slow cycles file missing, please first run getSlowCycles');
            load(slowCyclesFile,'TcycleMid','TcycleOffset','TcycleOnset'); %load data
            
            sharpWavesFile=[obj.currentAnalysisFolder filesep 'sharpWaves_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(sharpWavesFile,'Sharp wave file missing, please run getSharpWaves');
            load(sharpWavesFile);
            
            AIfile=[obj.currentAnalysisFolder filesep 'AI_2ms_4std.mat'];
            obj.checkFileRecording(sharpWavesFile,'AI file missing, please run getAI');
            load(AIfile);
            
            nNeurons=size(icA,2);
            if strcmp(cycleSelection,'first')
                pCycles=1:nCycles;
            elseif strcmp(cycleSelection,'rand')
                pCycles=randperm(numel(TcycleMid),nCycles);
            end
            
            %get high freq channel correlation matrix
            C=zeros(nCycles,nNeurons,nNeurons);
            for i=1:nCycles
                MSW=squeeze(BuildBurstMatrixA(icA,round(tA/2),IA,round(TcycleMid(pCycles(i))/2),round(TcycleOffset(pCycles(i))/2)))';
                tmpC=corrcoef(MSW);
                C(i,:,:)=tmpC;
            end
            ch=icA(1,:);
            avgCrossChCorr=squeeze(mean(C,1));
            clear IA tA icA;
            
            %cluster high freq channel correlation matrix
            [~,orderCtxDVR,ctxDVRClass]=DendrogramMatrix(avgCrossChCorr,'toPlotBinaryTree',0,'linkMethod','average','linkMetric','spearman','maxClusters',2);
            
            if mean(find(ctxDVRClass==1))<mean(find(ctxDVRClass==2))
                disp('Warning: cortex and DVR possition flipped in clustering, flipping position');
                ctxDVRClass(ctxDVRClass==1)=0;
                ctxDVRClass(ctxDVRClass==2)=1;
                ctxDVRClass(ctxDVRClass==0)=2;
            end
            %load spike sorting data
            load([obj.currentDataObj.recordingDir filesep obj.currentDataObj.recordingName '_spikeSort' filesep 'spikeSorting.mat']);

            MSO=BuildBurstMatrix(ic,round(t/binSO),round((TcycleOnset-preSO)/binSO),round(winSO/binSO));
            firingRateSO=squeeze(mean(MSO,1))';
            tMSO=(-preSO+binSO/2):binSO:(winSO-preSO+binSO/2);
            
            MSW=BuildBurstMatrix(ic,round(t/binSW),round((tSW'-preSW)/binSW),round(winSW/binSW));
            tMSW=(-preSW+binSW/2):binSW:(winSW-preSW-binSW/2);
            firingRateSW=squeeze(mean(MSW,1))';

            if numel(manQual)==size(ic,2)
                MUALabel=manQual;
            else
                MUALabel=2*ones(1,size(ic,2));
            end

            avg.CtxSWSU=[];
            avg.CtxSOSU=[];
            avg.CtxSWMU=[];
            avg.CtxSOMU=[];
            avg.DvrSWSU=[];
            avg.DvrSOSU=[];
            avg.DvrSWMU=[];
            avg.DvrSOMU=[];
            for i=1:numel(MUALabel)
                neuSpikeTemp=squeeze(allWaveforms(:,i,:));
                minV=min(neuSpikeTemp);
                [~,peakElec(i)]=min(minV);
                
                if ctxDVRClass(peakElec(i))==2
                    if MUALabel(i)==1
                        avg.CtxSWSU=[avg.CtxSWSU firingRateSW(:,i)];
                        avg.CtxSOSU=[avg.CtxSOSU firingRateSO(:,i)];
                    elseif MUALabel(i)==2
                        avg.CtxSWMU=[avg.CtxSWMU firingRateSW(:,i)];
                        avg.CtxSOMU=[avg.CtxSOMU firingRateSO(:,i)];
                    end
                else
                    if MUALabel(i)==1
                        avg.DvrSWSU=[avg.DvrSWSU firingRateSW(:,i)];
                        avg.DvrSOSU=[avg.DvrSOSU firingRateSO(:,i)];
                    elseif MUALabel(i)==2
                        avg.DvrSWMU=[avg.DvrSWMU firingRateSW(:,i)];
                        avg.DvrSOMU=[avg.DvrSOMU firingRateSO(:,i)];
                    end
                end
            end

            save(obj.files.spikeSTA,'parSpikeSTA','ctxDVRClass','avg','avgCrossChCorr','ch','orderCtxDVR','tMSW','tMSO','MUALabel');
        end
        
        %% getAI
        function data=getAI(obj,varargin)
            %% parameter and settings
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'selectedChannels',obj.currentDataObj.channelNumbers,@isnumeric);
            addParameter(parseObj,'stdAbsNoiseConstant',4,@isnumeric);
            addParameter(parseObj,'bin_ms',2,@isnumeric);
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
            
            %make parameter structure
            parAI=parseObj.Results;
            
            %check if analysis was already done done
            obj.files.AIfile=[obj.currentAnalysisFolder filesep 'AI_' num2str(bin_ms) 'ms_' num2str(stdAbsNoiseConstant) 'std.mat'];
            if exist(obj.files.AIfile,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.AIfile);
                else
                    disp('AI file already exists');
                end
                return;
            end

            [IA,tA,icA]=mcdIntensityExtractorUnfiltered(obj.currentDataObj,'Bin_ms',bin_ms,'SelectedChannels',selectedChannels,'StdAbsNoiseConstant',stdAbsNoiseConstant);
            
            save(obj.files.AIfile,'IA','tA','icA','parAI');
        end
        
        %% getHPSegments 
        function data=getHPSegments(obj,varargin)
            %% parameter and settings
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'win',80*1000,@isnumeric); %median filter window for extracting optic flow baseline
            addParameter(parseObj,'bin',1000,@isnumeric); %MAD (std) threshold for
            addParameter(parseObj,'artifactTreshHP',500,@isnumeric); %threshold in uV
            addParameter(parseObj,'maxSegments',500,@isnumeric); %MAD (std) threshold for
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'saveFileName',[]);
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
            
            %make parameter structure
            parHPSegments=parseObj.Results;
            
            %check if analysis was already done done
            if isempty(saveFileName)
                obj.files.HPSegments=[obj.currentAnalysisFolder filesep 'HPSegments_ch' num2str(ch) '.mat'];
            else
                obj.files.HPSegments=saveFileName;
            end
            if exist(obj.files.HPSegments,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.HPSegments);
                else
                    disp(['High pass segments file already exists']);
                end
                return;
            end
            
            animalStates=strsplit(obj.par.AnimalState{obj.currentPRec},'/');
            awakeStartTimeSec=obj.par.tStartAwake{obj.currentPRec};
            
            downSamplingFactor=obj.currentDataObj.samplingFrequency(1)*(bin/1000);
            for i=1:numel(animalStates)
                if strcmp(animalStates{i},'Awake') || strcmp(animalStates{i},'Running') || strcmp(animalStates{i},'Resting')
                    recDuration=obj.currentDataObj.recordingDuration_ms;
                    if ~isnan(awakeStartTimeSec)
                        allSegments=(awakeStartTimeSec*1000+win/2):win:(recDuration-win/2);
                    else
                        allSegments=(win/2):win:(recDuration-win/2);
                    end
                    nSeg=numel(allSegments);
                    tSeg=sort(allSegments(randperm(nSeg,min(maxSegments,nSeg))));
                    
                elseif strcmp(animalStates{i},'Sleep')
                    slowCyclesFile=[obj.currentAnalysisFolder filesep 'slowCycles_ch' num2str(ch) '.mat'];
                    obj.checkFileRecording(slowCyclesFile,'slow cycles file missing, please first run getSlowCycles');
                    load(slowCyclesFile); %load data
                    
                    tSeg=TcycleOnset;
                    nSeg=numel(tSeg);
                    tSeg=sort(tSeg(randperm(nSeg,min(maxSegments,nSeg))));
                end
                allAI{i}=zeros(numel(tSeg),win/bin);
                for j=1:numel(tSeg)
                    MF=obj.filt.FH2.getFilteredData(obj.currentDataObj.getData(ch,tSeg(j)-win/2,win));
                    if all(abs(MF)<artifactTreshHP)
                        AI=squeeze(mean(abs(reshape(MF,[downSamplingFactor  size(MF,3)/downSamplingFactor 1])),1))';
                        allAI{i}(j,:)=AI;
                    else
                        allAI{i}(j,:)=NaN;
                    end
                end
                allSeg{i}=tSeg;
                allStates{i}=animalStates{i};
            end
            
            save(obj.files.HPSegments,'allAI','allSeg','allStates','parHPSegments');
        end
        
        
        %% getHeartRateDBCorr
        function data=getHeartRateDBCorr(obj,varargin)
            %% parameter and settings
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'heartRateFile',[obj.currentAnalysisFolder filesep 'HR.mat'],@(x) exist(x,'file'));
            addParameter(parseObj,'interpTimeBin',1000);
            addParameter(parseObj,'stdWin',30*1000,@isnumeric); %median filter window for extracting optic flow baseline
            addParameter(parseObj,'plotWin',200*1000,@isnumeric); %MAD (std) threshold for 
            addParameter(parseObj,'plotBin',2*1000,@isnumeric);
            addParameter(parseObj,'InterpSmoothness',0.9,@isnumeric);
            addParameter(parseObj,'plotResults',[],@isnumeric);
            addParameter(parseObj,'hAxes',[],@isnumeric);
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
            
            %make parameter structure
            parheartRate=parseObj.Results;
            
            %check if analysis was already done done
            obj.files.heartRateAnal=[obj.currentAnalysisFolder filesep 'heartRateAnal.mat'];
            if exist(obj.files.heartRateAnal,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.heartRateAnal);
                else
                    disp(['Syncing DB with heart rate file already exists']);
                end
                return;
            end
            
            dbRatioFile=[obj.currentAnalysisFolder filesep 'dbRatio_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(dbRatioFile,'delta to beta file missing, please first run getDBRatio');
            load(dbRatioFile); %load data
            
            slowCyclesFile=[obj.currentAnalysisFolder filesep 'slowCycles_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(slowCyclesFile,'slow cycles file missing, please first run getSlowCycles');
            load(slowCyclesFile); %load data
            
            tStart=TcycleOnset(1);
            tEnd=TcycleOnset(end);
            
            tInterp=tStart:interpTimeBin:tEnd;
            
            HR=load([obj.currentAnalysisFolder filesep 'HR.mat']);
            
            interHR=csaps(HR.time,HR.bpm,InterpSmoothness,tInterp);
            interDB=csaps(t_ms,bufferedBetaRatio,InterpSmoothness,tInterp);
            interpStdHR = movingstd(interHR,stdWin/interpTimeBin);
            %plot(tInterp,interHR);hold on;plot(HR.time,HR.bpm,'r')
            
            HRStdRaster=BuildBurstMatrixA([1;1;1;numel(tInterp)],tInterp/plotBin,interpStdHR,(TcycleOnset-(plotWin/2))/plotBin,plotWin/plotBin);
            DBRaster=BuildBurstMatrixA([1;1;1;numel(tInterp)],tInterp/plotBin,interDB,(TcycleOnset-(plotWin/2))/plotBin,plotWin/plotBin);
            tRaster=(-plotWin/2+plotBin/2):plotBin:(plotWin/2);
            
            [~,pSmallRates]=sort(max(abs(HRStdRaster),[],3)); %sort according to hr ampitude
            nEvents=round(numel(pSmallRates)/3); %take only a third of events
            
            save(obj.files.heartRateAnal,'HRStdRaster','DBRaster','tRaster','tInterp','interHR','interDB','interpStdHR','parheartRate');

            if plotResults
                f=figure;
                if isempty(hAxes)
                    hAxes(1)=subaxis(f,2,2,1,'S',0.01,'MR',0.2);
                    hAxes(2)=subaxis(f,2,2,2,'S',0.01,'MR',0.2);
                    hAxes(3)=subaxis(f,2,2,3,'S',0.01,'MR',0.2);
                    hAxes(4)=subaxis(f,2,2,4,'S',0.01,'MR',0.2);
                    hold(hAxes(1),'on');
                    hold(hAxes(3),'on');
                end
                
                plot(tRaster/1000,normZeroOne(mean(squeeze(HRStdRaster))),'Parent',hAxes(1));
                plot(tRaster/1000,normZeroOne(mean(squeeze(DBRaster))),'r','Parent',hAxes(1));
                [hl,hO]=legend(hAxes(1),{'norm. HR variability','norm. \delta/\beta'},'Box','off');
                horizontalLegend(hO);
                hl.Position=[0.4887    0.8825    0.2875    0.0905];
                line([0 0],[0 1],'color',[0.8 0.8 0.8],'Parent',hAxes(1));
                
                plot(tRaster/1000,normZeroOne(mean(squeeze(HRStdRaster(pSmallRates(1:nEvents),1,:)))),'Parent',hAxes(3));
                plot(tRaster/1000,normZeroOne(mean(squeeze(DBRaster(pSmallRates(1:nEvents),1,:)))),'r','Parent',hAxes(3));
                [hl,hO]=legend(hAxes(3),{'norm. HR variability','norm. \delta/\beta'},'Box','off');
                horizontalLegend(hO);
                hl.Position=[0.4887    0.8825    0.2875    0.0905];
                line([0 0],[0 1],'color',[0.8 0.8 0.8],'Parent',hAxes(3));
                
                imagesc(tRaster/1000,1:size(DBRaster(pSmallRates(1:nEvents),1,:),1),squeeze(DBRaster(pSmallRates(1:nEvents),1,:)),'Parent',hAxes(2));
                ylabel('Cycle #','Parent',hAxes(2));
                set(hAxes(2),'XTickLabel',[]);
                xlabel('Time [s]');
                
                imagesc(tRaster/1000,1:size(HRStdRaster(pSmallRates(1:nEvents),1,:),1),squeeze(HRStdRaster(pSmallRates(1:nEvents),1,:)),'Parent',hAxes(4));
                ylabel('Cycle #','Parent',hAxes(4));
                set(hAxes(4),'XTickLabel',[]);
                xlabel('Time [s]');
                
                %subplot(2,2,2);crosscorr(interHR,mStdHR,500);
            end
            
        end
        
        %% plotSyncedDBEyeMovements
        function hOut=plotSyncedDBEyeMovements(obj,varargin)
            %% parameter and settings
            obj.checkFileRecording;
            
            parseObj = inputParser;
            parseObj.FunctionName='sleepAnalysis\plotSyncedDBEyeMovements';
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'videoFile',[obj.currentVideosFolder filesep obj.par.VideoFiles{obj.currentPRec}],@(x) exist(x,'file'));
            addParameter(parseObj,'saveFigures',1,@isnumeric);
            addParameter(parseObj,'rLim4Rose',[],@isnumeric);
            addParameter(parseObj,'RoseAlpha',0.9,@isnumeric);
            addParameter(parseObj,'noBackground',0,@isnumeric);
            addParameter(parseObj,'printLocalCopy',0,@isnumeric);
            addParameter(parseObj,'h',0,@ishandle);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            addParameter(parseObj,'plotRandomDist',1,@isnumeric);
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
            parPlotSyncDBEye=parseObj.Results;
            
            [~,videoFileName]=fileparts(videoFile);

            syncDBEyeFile=[obj.currentAnalysisFolder filesep 'syncDBEye_' videoFileName '.mat'];
            obj.checkFileRecording(syncDBEyeFile,'Eye sync analysis missing, please first run getSyncedDBEyeMovements');
            load(syncDBEyeFile); %load data
            
            mResampledTemplate=mean(resampledTemplate);
            
            if h==0
                fH=figure;
                h=axes;
            else
                saveFigures=0;
                axes(h);
            end
            cMap=lines(8);
            
            if ~isempty(rLim4Rose)
                hTmp = polarTight(0, rLim4Rose);
                delete(hTmp)
                set(h, 'Nextplot','add');hold on;
            end

            hOut.hRose=rose(phaseMov*2*pi-mPhaseDB,parSyncDBEye.nBins);
            hOut.hRose.Color=[0.9 0.078 0.184];
            XdataRose = get(hOut.hRose,'Xdata');XdataRose=reshape(XdataRose,[4,numel(XdataRose)/4]);
            YdataRose = get(hOut.hRose,'Ydata');YdataRose=reshape(YdataRose,[4,numel(YdataRose)/4]);
            hOut.hPatch=patch(XdataRose,YdataRose,[0.9 0.078 0.184]);
            set(hOut.hPatch,'FaceAlpha',RoseAlpha);
            %set(h,'color','k');
            maxSamplesInBin=max(max(sqrt(XdataRose.^2+YdataRose.^2)));hold on;
            
            hOut.hPolar=polar([0 (1:parSyncDBEye.nBins)/parSyncDBEye.nBins]*pi*2-mPhaseDB,[mResampledTemplate(end) mResampledTemplate]/(max(mResampledTemplate/maxSamplesInBin)));
            hOut.hPolar.LineWidth=2;
            hOut.hPolar.Color=cMap(1,:,:);
            
            uistack(hOut.hPatch, 'top');
            
            delete(findall(h, 'String', '30', '-or','String','60', '-or','String','120', '-or','String','150', '-or','String','210', '-or','String','240', '-or','String','300', '-or','String','330'));
            
            hOut.l=legend([hOut.hRose hOut.hPolar],'SWC','OF');
            hOut.l.Color=[1 1 1];
            hOut.l.Box='off';
                        
            if plotRandomDist
                hOut.hRose2=rose(phaseRand*2*pi-mPhaseDB,parSyncDBEye.nBins);
                hOut.hRose2.Color=[0.5 0.5 0.5];
            end
            
            %if ~isempty(rLim4Rose)
            %    set(h_fake,'Visible','off');
            %end
            
            if saveFigures
                set(fH,'PaperPositionMode','auto');
                fileName=[obj.currentPlotFolder filesep 'syncEye_' videoFileName];
                print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                if printLocalCopy
                    fileName=[cd filesep obj.par.Animal{obj.currentPRec} '_Rec' num2str(obj.currentPRec) '_syncEyeDB_' videoFileName];
                    print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                    print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                end
            end
            
        end
        
        %% plotSyncedDBEyeMovementsRaster
        function hOut=plotSyncedDBEyeMovementsRaster(obj,varargin)
            %% parameter and settings
            hOut=[];
            obj.checkFileRecording;
            
            parseObj = inputParser;
            parseObj.FunctionName='sleepAnalysis\plotSyncedDBEyeMovementsRaster';
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'videoFile',[obj.currentVideosFolder filesep obj.par.VideoFiles{obj.currentPRec}],@(x) exist(x,'file'));
            addParameter(parseObj,'saveFigures',1,@isnumeric);
            addParameter(parseObj,'printLocalCopy',0,@isnumeric);
            addParameter(parseObj,'h',0,@ishandle);
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
            
            %make parameter structure
            parPlotSyncDBEye=parseObj.Results;
            
            [~,videoFileName]=fileparts(videoFile);

            syncDBEyeFile=[obj.currentAnalysisFolder filesep 'syncDBEye_' videoFileName '.mat'];
            obj.checkFileRecording(syncDBEyeFile,'Eye sync analysis missing, please first run getSyncedDBEyeMovements');
            load(syncDBEyeFile); %load data
            
            if h==0
                fH=figure;
                h(1)=subaxis(fH,4,2,1,1,1,3,'S',0.01);
                h(2)=subaxis(fH,4,2,1,4,1,1,'S',0.01);
                h(3)=subaxis(fH,4,2,2,1,1,3,'S',0.01);
                h(4)=subaxis(fH,4,2,2,4,1,1,'S',0.01);
            else
                saveFigures=0;
            end
            
            if numel(h)>=1
                edges=(0:parSyncDBEye.nBins)/(parSyncDBEye.nBins-0.0000001);
                middles=(edges(1:end-1)+edges(2:end))/2;
                
                axes(h(1));
                hOut.imagesc=imagesc(0:parSyncDBEye.nBins,1:size(resampledTemplate,1),resampledTemplate);hold on;
                set(h(1),'XTickLabel',[]);
                ylabel('# cycle');
                p=cell2mat(cellfun(@(x) ~isempty(x),phaseAll,'UniformOutput',0));
                for i=find(p)
                    hOut.hP(i)=plot(phaseAll{i}*parSyncDBEye.nBins,i*ones(size(phaseAll{i})),'.r');
                end

                I=histc(phaseMov,edges);
                
                if numel(h)==1
                    xlabel('Phase');
                    set(h,'XTick',[-0.5 parSyncDBEye.nBins+0.5],'XTickLabel',{'0','2\pi'});
                    xlim([-0.5 parSyncDBEye.nBins+0.5])
                end
            end
            
            if numel(h)>=2
                axes(h(2));
                hOut.p1=plot(middles,normZeroOne(mean(resampledTemplate)),'lineWidth',2);hold on;
                hOut.p2=plot(middles,normZeroOne(I(1:end-1)),'r','lineWidth',2);
                xlim([0 1]);
                xlabel('Phase');
                set(h(2),'XTick',[0 1],'XTickLabel',{'0','2\pi'});
                hOut.l=legend('norm. \delta/\beta','norm. OF counts');
                hOut.l.Box='off';
                hOut.l.Position=[0.6434    0.9061    0.2596    0.0812];
            end
            
            if numel(h)>=3
                axes(h(3));
                imagesc(resampledTemplate);hold on;
                set(h(3),'YTickLabel',[],'XTickLabel',[]);
                p=cell2mat(cellfun(@(x) ~isempty(x),phaseAllRand,'UniformOutput',0));
                for i=find(p)
                    plot(phaseAllRand{i}*parSyncDBEye.nBins,i*ones(size(phaseAllRand{i})),'*r');
                end
                I=histc(phaseRand,edges);
            end
            
            if numel(h)>=4
                axes(h(4));
                hOut.p3=plot(middles,normZeroOne(mean(resampledTemplate)),'lineWidth',2);hold on;
                hOut.p4=plot(middles,normZeroOne(I(1:end-1)),'r','lineWidth',2);
                xlim([0 1]);
                xlabel('Phase');
                set(h(4),'XTick',[0 1],'XTickLabel',{'0','2\pi'},'YTickLabel',[])
            end
            
            if saveFigures
                set(fH,'PaperPositionMode','auto');
                fileName=[obj.currentPlotFolder filesep 'syncEyeDBRaster_' videoFileName];
                print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                if printLocalCopy
                    fileName=[cd filesep obj.par.Animal{obj.currentPRec} '_Rec' num2str(obj.currentPRec) '_syncEyeDBRaster_' videoFileName];
                    print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                    print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                end
            end
            
        end
        
        %timeBinDB=(parDBRatio.movWin-parDBRatio.movOLWin);

            %{
            figure;
            subplot(1,2,1);
            imagesc(resampledTemplate);hold on;
            p=cell2mat(cellfun(@(x) ~isempty(x),phaseAll,'UniformOutput',0));
            for i=find(p)
                plot(phaseAll{i}*nBins,i*ones(size(phaseAll{i})),'or');
            end
            
            subplot(1,2,2);
            imagesc(resampledTemplate);hold on;
            p=cell2mat(cellfun(@(x) ~isempty(x),phaseAllRand,'UniformOutput',0));
            for i=find(p)
                plot(phaseAllRand{i}*nBins,i*ones(size(phaseAllRand{i})),'or');
            end
            %}
            
        %% plotEyeVideoOFDB
        function obj=plotEyeVideoOFDB(obj,varargin)
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'winFrame',100,@isnumeric);
            addParameter(parseObj,'tStart',[],@isnumeric);
            addParameter(parseObj,'tEnd',[],@isnumeric);
            addParameter(parseObj,'outputVideo',[]);
            addParameter(parseObj,'outputFrameRate',10);
            addParameter(parseObj,'opticFlowFile',[]);
            addParameter(parseObj,'OFlineColor','black');
            addParameter(parseObj,'ampOFLine',50);
            addParameter(parseObj,'showOnlyEye',true);
            addParameter(parseObj,'videoCompressor','DV Video Encoder');
            addParameter(parseObj,'saveVideo',false,@isnumeric);
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'videoFile',[obj.currentVideosFolder filesep obj.par.VideoFiles{obj.currentPRec}],@(x) exist(x,'file'));
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
                            
            dbRatioFile=[obj.currentAnalysisFolder filesep 'dbRatio_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(dbRatioFile,'delta to beta file missing, please first run getDBRatio');
            DB=load(dbRatioFile); %load data
            
            [~,videoFileName]=fileparts(videoFile);
            eyeTrackingFile=[obj.currentAnalysisFolder filesep 'eyeTracking_' videoFileName '.mat'];
            obj.checkFileRecording(eyeTrackingFile,'Eye tracking analysis missing, please first run getEyeMovement');
            ET=load(eyeTrackingFile,'parEyeTracking','pFrames','mOF','pbboxUpdate','bboxCenterAll'); %load data
            
            syncDBEyeFile=[obj.currentAnalysisFolder filesep 'syncDBEye_' videoFileName '.mat'];
            obj.checkFileRecording(syncDBEyeFile,'sync eye to beta 2 delta file missing, please first run getSyncedDBEyeMovements');
            sync=load(syncDBEyeFile); %load data
            
            if ~showOnlyEye
                videoReader = VideoReader(videoFile); %initiate video obj since number of frames was already read (not allowed by matlab)
            end
            %videoPlayer  = vision.VideoPlayer('Position',[100 100 [videoReader.Width, videoReader.Height]+30]);
            if saveVideo
                if isempty(outputVideo)
                    outputVideo=[obj.currentAnalysisFolder filesep videoFileName 'OFDB.avi'];
                end
                videoWriter = vision.VideoFileWriter([outputVideo '.avi'],'FrameRate',outputFrameRate);
                videoWriter.VideoCompressor=videoCompressor;
            end
            
            cMapLines=lines(8);
            dPix=5;
            if ~isempty(opticFlowFile)
                OF=load(opticFlowFile);
                tFramesOF=sync.tVideoFrames(OF.pFrames);
                plotFlowField=true;
                
                rX=round(dPix/2):dPix:size(OF.allOF,2)-round(dPix/2);
                rY=round(dPix/2):dPix:size(OF.allOF,1)-round(dPix/2);
                [Y, X] = meshgrid(rX,rY);
                if ~showOnlyEye
                    X=X+OF.initialFrameSubregion(2);
                    Y=Y+OF.initialFrameSubregion(1);
                end
                
                shapes = vision.ShapeInserter;
                shapes.Shape = 'Lines';
                shapes.BorderColor = OFlineColor;
            else
                plotFlowField=false;
            end
            
            
            smoothDB=2e-11;
            interpDBOF=csaps(DB.t_ms,DB.bufferedBetaRatio,smoothDB,tFramesOF);

            %plot(tAllFrames(pFrames),interpOF);hold on;plot(tAnalyzedFrames,validmOF);
            if isempty(tStart)
                tStart=tFramesOF(1);
            end
            if isempty(tEnd)
                tEnd=tFramesOF(end);
            end
            
            pFrames=find(tFramesOF>=tStart & tFramesOF<=tEnd);
            tFrames=OF.pFrames(pFrames);
            
            interpDB=interpDBOF(pFrames);
            interpOF=OF.mOF(pFrames);
            pFramesOrig=OF.pFrames(pFrames);
            %plot(normZeroOne(interpDB));hold on;plot(normZeroOne(interpOF));
            
            eInterpDB=[zeros(1,winFrame) interpDB zeros(1,winFrame)]./std([zeros(1,winFrame) interpDB zeros(1,winFrame)]);
            eInterpOF=[zeros(1,winFrame) interpOF zeros(1,winFrame)]./2/std([zeros(1,winFrame) interpOF zeros(1,winFrame)]);
            %plot(eInterpDB);hold on;plot(eInterpOF)
            %plot(tAllFrames(pFrames),interpDB);hold on;plot(t_ms,bufferedBetaRatio);xlim([16027110.8963877          18570778.2072246])
            %plot(tAllFrames(pFrames),interpOF);hold on;plot(tAnalyzedFrames,validmOF);xlim([16027110.8963877          18570778.2072246])

            %set scaling parameters for curves
            if showOnlyEye
                f=figure('position',[100 100 350 600]);
                set(gcf,'PaperPositionMode','auto');
                
                videoFrame=squeeze(OF.allIm(:,:,1));
                h(1)=subaxis(f,2,1,1,'M',0.03,'S',0.07);
                h(2)=subaxis(f,2,1,2,'M',0.03,'S',0.07);
                axis(h(1),'off');
                xlim(h(1),[0 2*winFrame]);
                yl=[0 4];
                ylim(h(1),yl);
                imshow(videoFrame,'Parent',h(2));
                set(h(2),'nextplot','replacechildren');
                hold(h(1),'on');
                pH=[];
            else
                f=figure('position',[100 100 500 500]);
                videoFrame=rgb2gray(videoReader.read(pFramesOrig(i))); %read will be replaced by readFrame in future versions but it is not possible to skip frames with readframes
                
                W=videoReader.Width;
                H=videoReader.Height;
                pX=(1:(2*winFrame+1))/(2*winFrame+1)*W;
                yStartPixDB=H*0.4; %from top down
                yPixDB=H*0.05;
                yStartPixOF=H*0.5;
                yPixOF=H*0.05;
                ylineLim=H*0.55;
                imshow(videoFrame);
                set(gca,'nextplot','replacechildren');
            end
            set(f,'Renderer','zbuffer');
            
            for i=1:numel(pFramesOrig)
                tmpDB=eInterpDB(i:i+2*winFrame);
                tmpOF=eInterpOF(i:i+2*winFrame);
                
                if showOnlyEye
                    videoFrame=squeeze(OF.allIm(:,:,pFrames(i)));
                else
                    videoFrame=rgb2gray(videoReader.read(pFramesOrig(i))); %read will be replaced by readFrame in future versions but it is not possible to skip frames with readframes
                end

                if plotFlowField
                    currentFrame=squeeze(OF.allOF(:,:,pFrames(i)));
                    currentFrame=currentFrame(rY,rX);
                    Hor = imag(currentFrame)*ampOFLine;
                    Ver = real(currentFrame)*ampOFLine;
                    
                    OFlines = [Y(:)'; X(:)'; Y(:)'+Ver(:)'; X(:)'+Hor(:)'];
                    videoFrame = step(shapes, videoFrame,  int32(OFlines)');
                    % Draw lines on top of image
                    %pHL=line([X(:) X(:)+Hor(:)]',[Y(:) Y(:)+Ver(:)]','color',OFlineColor);
                end
                
                if showOnlyEye
                    delete(pH);
                    imshow(videoFrame,'Parent',h(2));
                    pH(1)=plot(h(1),tmpDB,'lineWidth',1,'color',cMapLines(1,:));
                    pH(2)=plot(h(1),tmpOF,'lineWidth',1,'color',cMapLines(2,:));
                    pH(3)=line([winFrame+1 winFrame+1],yl,'color',cMapLines(5,:),'Parent',h(1));
                    %pH(4)=text(170,-0.4,[num2str(i-1) 's'],'Parent',h(1),'FontSize',16,'FontWeight','Bold');
                else
                    imshow(videoFrame);hold on;
                    pH(1)=plot(pX,-tmpDB*yPixDB+yStartPixDB,'lineWidth',3,'color',cMapLines(1,:));
                    pH(2)=plot(pX,-tmpOF*yPixOF+yStartPixOF,'lineWidth',3,'color',cMapLines(2,:));
                    pH(3)=line([pX(winFrame+1) pX(winFrame+1)],[ylineLim 0],'color',cMapLines(5,:));
                end
                
                if saveVideo %save tracked video
                    frame=getframe(f);
                    step(videoWriter, frame.cdata);
                end
                
            end
            
            if saveVideo
                release(videoWriter);
            end
            if ~showOnlyEye
                delete(videoReader);
            end
            delete(f);
            
        end
        
        %% getSyncedDBEyeMovements
        function data=getSyncedDBEyeMovements(obj,varargin)
            %% parameter and settings
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'videoFile',[obj.currentVideosFolder filesep obj.par.VideoFiles{obj.currentPRec}],@(x) exist(x,'file'));
            addParameter(parseObj,'matroxTrigScheme',obj.par.MatroxTrigScheme{obj.currentPRec});
            addParameter(parseObj,'win',180*1000,@isnumeric); %median filter window for extracting optic flow baseline
            addParameter(parseObj,'nStd',6,@isnumeric); %MAD (std) threshold for 
            addParameter(parseObj,'nBins',18,@isnumeric);
            addParameter(parseObj,'pixelMoveThresh',10,@isnumeric);
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'nFramesRemoveAfterROIShift',3,@isnumeric);
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
            
            %make parameter structure
            parSyncDBEye=parseObj.Results;
            
            [~,videoFileName]=fileparts(videoFile);
            %check if analysis was already done done
            obj.files.syncDBEye=[obj.currentAnalysisFolder filesep 'syncDBEye_' videoFileName '.mat'];
            if exist(obj.files.syncDBEye,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.syncDBEye);
                else
                    disp(['Syncing DB with eye tracking file for video: ' videoFileName ' already exists']);
                end
                return;
            end
            
            eyeTrackingFile=[obj.currentAnalysisFolder filesep 'eyeTracking_' videoFileName '.mat'];
            obj.checkFileRecording(eyeTrackingFile,'Eye tracking analysis missing, please first run getEyeMovement');
            load(eyeTrackingFile,'parEyeTracking','pFrames','mOF','pbboxUpdate','bboxCenterAll'); %load data
            
            digiTrigFile=[obj.currentAnalysisFolder filesep 'digiTrig.mat'];
            obj.checkFileRecording(digiTrigFile,'digital trigger file missing, please first run getDigiData');
            load(digiTrigFile); %load data
            
            dbRatioFile=[obj.currentAnalysisFolder filesep 'dbRatio_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(dbRatioFile,'delta to beta file missing, please first run getDBRatio');
            load(dbRatioFile); %load data
            
            slowCyclesFile=[obj.currentAnalysisFolder filesep 'slowCycles_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(slowCyclesFile,'slow cycles file missing, please first run getSlowCycles');
            load(slowCyclesFile); %load data
            
            %remove frames that are close to a ROI shift and frames with large shifts
            p2RemoveShifts=find(sqrt(diff(bboxCenterAll(:,1)).^2+diff(bboxCenterAll(:,2)).^2)>pixelMoveThresh)+1;
            pFramesValid=pFrames;
            validmOF=mOF;
            if ~isempty(pbboxUpdate) || ~isempty(p2RemoveShifts)
                p2Remove=union(p2RemoveShifts,pbboxUpdate)';
                p2Remove=bsxfun(@plus,p2Remove,(0:nFramesRemoveAfterROIShift-1)');
                p2Remove=unique(p2Remove(:));
                
                validmOF(p2Remove)=[];
                pFramesValid(p2Remove)=[];
                bboxCenterAll(p2Remove,:)=[];
            end

            startFrameInVideo=cellfun(@(x) str2num(x),strsplit(obj.par.FrameRange{obj.currentPRec},'-'),'UniformOutput', 1);
            if strcmp(matroxTrigScheme,'startD2FramesD1NoEndNoMissed') %correct start trigger (on digi ch 2), no stop trigger, no missed frames (on digi ch 1)
                pMatrox=1;
                pStart=2;
                tAllFrames=tTrig{pMatrox}(1:2:end);
                tAllFrames(tAllFrames<tTrig{pStart}(1))=[];
                nFrames=pFramesValid(end);
                tAllFrames=tAllFrames(1:nFrames);
                tVideoFrames=tAllFrames(startFrameInVideo(1):startFrameInVideo(2));
                tAnalyzedFrames=tAllFrames(pFramesValid+startFrameInVideo(1)-1);
            elseif strcmp(matroxTrigScheme,'D10FromMatroxNoMissed')
                pMatrox=1;
                pMatrox10=10;
                tAllFrames=tTrig{pMatrox}(1:2:end);
                tAllFrames(tAllFrames<tTrig{pMatrox10}(1))=[];
                lastStr=strsplit(tTrig_string{pMatrox10}{end},';');
                lastFrameFromStr=str2num(lastStr{1}(8:end));
                trigNumDiff=find(tAllFrames<tTrig{pMatrox10}(end),1,'last')~=lastFrameFromStr;
                if trigNumDiff~=0
                    warning(['Number of triggers in T1 differs by ' num2str(trigNumDiff) ' trig than what is found from 10 strings-check for missed frames!!!!']);
                end
                tVideoFrames=tAllFrames(startFrameInVideo(1):startFrameInVideo(2));
                tAnalyzedFrames=tAllFrames(pFramesValid+startFrameInVideo(1)-1);
            elseif strcmp(matroxTrigScheme(1:11),'Lizard4Case')
                stringNames=tTrig_string{10};
                tmp=cellfun(@(x) x(8:15),stringNames,'Uniformoutput',0);
                tmp1=cellfun(@(x) find(x==';'),tmp,'Uniformoutput',0);
                frameNumber=cell2mat(cellfun(@(x,y) str2double(x(1:(y-1))),tmp,tmp1,'Uniformoutput',0));
                pStartFrame=find(frameNumber==0);
                pEndFrame=pStartFrame(3)-1;
                pStartFrame=pStartFrame(2);
                tAllFrames=interp1(frameNumber(pStartFrame:pEndFrame),tTrig{10}(pStartFrame:pEndFrame),0:frameNumber(pEndFrame),'linear');
                tVideoFrames=tAllFrames(startFrameInVideo(1):startFrameInVideo(2));
                tAnalyzedFrames=tAllFrames(pFramesValid+startFrameInVideo(1)-1);
            elseif strcmp(matroxTrigScheme,'interpFrom10')
                stringNames=tTrig_string{10};
                tmp=cellfun(@(x) x(8:15),stringNames,'Uniformoutput',0);
                tmp1=cellfun(@(x) find(x==';'),tmp,'Uniformoutput',0);
                frameNumber=cell2mat(cellfun(@(x,y) str2double(x(1:(y-1))),tmp,tmp1,'Uniformoutput',0));
                tAllFrames=interp1(frameNumber,tTrig{10},0:frameNumber(end),'linear');
                tVideoFrames=tAllFrames(startFrameInVideo(1):startFrameInVideo(2));
                tAnalyzedFrames=tAllFrames(pFramesValid+startFrameInVideo(1)-1);
            elseif strcmp(matroxTrigScheme,'xxxx')
                
                pCameraStart=tTrig{10}(1);
                pCameraLast=tTrig{10}(end);
                
                frameNumberStartString=tTrig_string{10}{1};
                frameNumberEndString=tTrig_string{10}{end};
                
                frameNumberStart=str2num(frameNumberStartString(8:find(frameNumberStartString==';',1,'first')-1));
                frameNumberEnd=str2num(frameNumberEndString(8:find(frameNumberEndString==';',1,'first')-1));
                
                frameShutterTimes=tTrig{1}(find(tTrig{1}>pCameraStart & tTrig{1}<pCameraLast));
                frameShutterTimes(diff(frameShutterTimes)<2)=[]; %remove off trigger - every trigger is 1ms long and is recorded twice (once for on and once for off)
                
                frameShutterTimesAfterLastTrigger=Ts{1}(find(Ts{1}>pCameraLast));
                frameShutterTimesAfterLastTrigger(diff(frameShutterTimesAfterLastTrigger)<2)=[]; %remove off trigger - every trigger is 1ms long and is recorded twice (once for on and once for off)
                frameShutterTimesAfterLastTrigger=frameShutterTimesAfterLastTrigger(1:(nFrames-frameNumberEnd));
                frameShutterTimes=[frameShutterTimes frameShutterTimesAfterLastTrigger];
            end
            %plot(tAnalyzedFrames/3600000,normZeroOne(validmOF));hold on;plot(tAnalyzedFrames/3600000,normZeroOne(bboxCenterAll));
            winSamples=round(win/1000*(parEyeTracking.frameRate/parEyeTracking.skipFrames));
            mOFmed=fastmedfilt1d(validmOF,winSamples)';
            mOFMAD=fastmedfilt1d(abs(validmOF-mOFmed),winSamples)'*1.4826;
            tMovement=tAnalyzedFrames(validmOF>(mOFmed+nStd*mOFMAD));
            %plot(tAnalyzedFrames/3600000,validmOF);hold on;plot(tAnalyzedFrames/3600000,mOFmed+nStd*mOFMAD);
            for i=1:numel(TcycleOnset)
                cycleDuration=TcycleOffset(i)-TcycleOnset(i);
                pTmp=find(tMovement>(TcycleMid(i)-cycleDuration/2) & tMovement<(TcycleMid(i)+cycleDuration/2));
                phaseAll{i}=(tMovement(pTmp)-(TcycleMid(i)-cycleDuration/2))/cycleDuration;
                
                shufTimes=rand(1,numel(pTmp))*cycleDuration;
                phaseAllRand{i}=shufTimes/cycleDuration;
                
                pTmp=find(t_ms>(TcycleMid(i)-cycleDuration/2) & t_ms<(TcycleMid(i)+cycleDuration/2));
                resampledTemplate(i,:) = interp1((0:(numel(pTmp)-1))./(numel(pTmp)-1),bufferedBetaRatio(pTmp)',(0:(nBins-1))/(nBins-1),'spline');

                %{
                cycleDuration=TcycleOffset(i)-TcycleOnset(i);
                pTmp=find(tMovement>TcycleOnset(i) & tMovement<TcycleOffset(i));
                phaseAll{i}=(tMovement(pTmp)-TcycleOffset(i))/(TcycleOffset(i)-TcycleOnset(i));
                
                shufTimes=rand(1,numel(pTmp))*(TcycleOffset(i)-TcycleOnset(i));
                phaseAllRand{i}=shufTimes/(TcycleOffset(i)-TcycleOnset(i));
                
                pTmp=find(t_ms>TcycleOnset(i) & t_ms<TcycleOffset(i));
                resampledTemplate(i,:) = interp1((0:(numel(pTmp)-1))./(numel(pTmp)-1),bufferedBetaRatio(pTmp)',(0:(nBins-1))/(nBins-1),'spline');

                %}
            end
            phaseMov=cell2mat(phaseAll);
            phaseRand=cell2mat(phaseAllRand);
            
            mPhaseMov=angle(mean(exp(1i*phaseMov*2*pi))); %Mean of circular quantities - wiki
            binCenters=(0:(nBins))/(nBins);binCenters=(binCenters(1:end-1)+binCenters(2:end))/2;
            mPhaseDB=angle(mean(mean(resampledTemplate).*exp(1i.*binCenters*2*pi))); %Mean of circular quantities - wiki
            
            save(obj.files.syncDBEye,'phaseMov','mPhaseDB','mPhaseMov','phaseAll','phaseAllRand','phaseRand','resampledTemplate','validmOF','tAnalyzedFrames','tVideoFrames','tAllFrames','tMovement','parSyncDBEye','pFramesValid','mOFmed','mOFMAD');
           
        end
        
        %% getEyeMovements
        function [data]=getEyeMovements(obj,varargin)
            %% parameter and settings
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'videoFile',[obj.currentVideosFolder filesep obj.par.VideoFiles{obj.currentPRec}],@(x) exist(x,'file'));
            addParameter(parseObj,'startFrame',1,@isnumeric); %max freq. to examine
            addParameter(parseObj,'endFrame',Inf,@isnumeric);
            addParameter(parseObj,'initialFrameSubregion',[],@isnumeric);
            addParameter(parseObj,'frameForEyePosExtraction',[],@isnumeric);
            addParameter(parseObj,'fractionOfBoxJumpThreshold',0.25,@isnumeric);
            addParameter(parseObj,'saveFullOFMatrices',false,@isnumeric);
            addParameter(parseObj,'loadInitialConditions',true,@isnumeric);
            addParameter(parseObj,'skipFramesBoundingBox',30,@isnumeric);
            addParameter(parseObj,'removeBorderOF',true);
            addParameter(parseObj,'borderPix',5);
            addParameter(parseObj,'skipFrames',10,@isnumeric);
            addParameter(parseObj,'plotTracking',false,@isnumeric);
            addParameter(parseObj,'saveTrackingVideo',false,@isnumeric);
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'savedFileName',[]);
            addParameter(parseObj,'minTrackingPoints',50,@isnumeric);
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

            %make parameter structure
            parEyeTracking=parseObj.Results;
            
            [~,videoFileName]=fileparts(videoFile);
            if isempty(savedFileName)
                obj.files.eyeTracking=[obj.currentAnalysisFolder filesep 'eyeTracking_' videoFileName '.mat'];
                trackingFileName=[obj.currentAnalysisFolder filesep videoFileName '_EyeTracking.avi'];
            else
                obj.files.eyeTracking=[savedFileName '.mat'];
                trackingFileName=[savedFileName '_EyeTracking.avi'];
            end
            
            %check if analysis was already done done
            if exist(obj.files.eyeTracking,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.eyeTracking);
                else
                    disp(['Eye tracking file for video: ' videoFileName 'already exists']);
                end
                return;
            end
            
            if exist(obj.files.eyeTracking,'file') & loadInitialConditions
                load(obj.files.eyeTracking,'startFrame','initialFrameSubregion');
                parEyeTracking.startFrame=startFrame;
                parEyeTracking.initialFrameSubregion=initialFrameSubregion;
            end
            
            if saveTrackingVideo
                plotTracking=1;
            end
            %% Pre processing
            videoReader = VideoReader(videoFile); %initiate video obj
            nFramesVideo=videoReader.NumberOfFrames;
            frameWidth=videoReader.Width;
            frameHeight=videoReader.Height;
            frameRate=videoReader.FrameRate;
            
            parEyeTracking.nFramesVideo=nFramesVideo;
            parEyeTracking.frameWidth=frameWidth;
            parEyeTracking.frameHeight=frameHeight;
            parEyeTracking.frameRate=frameRate;
            
            if isinf(endFrame) %analyze the complete video
                endFrame=nFramesVideo;
            end
            pFrames=startFrame:skipFrames:endFrame;
            nFrames=numel(pFrames);
            
            if isempty(frameForEyePosExtraction)
                frameForEyePosExtraction=pFrames(1);
            end
            
            initFrame = rgb2gray(read(videoReader,frameForEyePosExtraction));% videoReader.CurrentTime=(1/videoReader.FrameRate)*(pFrames(1)-1);
            if isempty(initialFrameSubregion) %to manually select region for extracting eye movements
                f=figure('position',[100 100 1200 600]);
                subplot(1,3,1:2);imshow(initFrame);
                h = imrect(gca);
                initialFrameSubregion=h.getPosition;
                xInd=round(initialFrameSubregion(1):(initialFrameSubregion(1)+initialFrameSubregion(3)));
                yInd=round(initialFrameSubregion(2):(initialFrameSubregion(2)+initialFrameSubregion(4)));
                subplot(1,3,3);imshow(initFrame(yInd,xInd,:));
                title('Selected region - press any key');
                pause;
                close(f);
                initialFrameSubregion=round(initialFrameSubregion);
            else
                xInd=round(initialFrameSubregion(1):(initialFrameSubregion(1)+initialFrameSubregion(3)));
                yInd=round(initialFrameSubregion(2):(initialFrameSubregion(2)+initialFrameSubregion(4)));
            end
            delete(videoReader);
            
            parEyeTracking.initialFrameSubregion=initialFrameSubregion;
            
            if saveFullOFMatrices %if save all OP matrices initialize array
                allOF=zeros(numel(yInd),numel(xInd),nFramesVideo,'like',complex(zeros(1,'single'),zeros(1,'single')));
            end
            %determine the position of border pixels to remove from OF analysis
            if removeBorderOF
                tmp=zeros([numel(yInd),numel(xInd)]);
                tmp(:)=1:numel(tmp);
                pSizeBorder=tmp([1:borderPix end-borderPix+1:end],:);
                pUpDownBorder=tmp(:,[1:borderPix end-borderPix+1:end]);
                pBorder=unique([pSizeBorder(:);pUpDownBorder(:)]);
            end
            
            %defition of optic flow and video reader/converter objects
            if skipFrames~=1 || pFrames(1)>1
                videoReader = VideoReader(videoFile); %initiate video obj since number of frames was already read (not allowed by matlab)
                %videoReader.CurrentTime=(1/videoReader.FrameRate)*(pFrames(1)-1);
                nonConsecutiveVideo=true;
            else
                videoReader = vision.VideoFileReader(videoFile,'ImageColorSpace','Intensity','VideoOutputDataType','uint8'); % create required video objects
                nonConsecutiveVideo=false;
            end
            
            % optic flow definitions
            converter = vision.ImageDataTypeConverter;
            opticalFlow = vision.OpticalFlow(...
                'Method','Lucas-Kanade',...
                'ReferenceFrameDelay', 1,... 
                'OutputValue' ,'Horizontal and vertical components in complex form');% use of the Lucas-Kanade method for optic flow determination
            
            bboxPoints=[initialFrameSubregion(1) initialFrameSubregion(2);initialFrameSubregion(1) initialFrameSubregion(2)+initialFrameSubregion(4);initialFrameSubregion(1)+initialFrameSubregion(3) initialFrameSubregion(2)+initialFrameSubregion(4);initialFrameSubregion(1)+initialFrameSubregion(3) initialFrameSubregion(2)];                
            bboxCenter=[(bboxPoints(3,1)+bboxPoints(1,1))/2 (bboxPoints(3,2)+bboxPoints(1,2))/2];
            bboxCenterOld=[(bboxPoints(3,1)+bboxPoints(1,1))/2 (bboxPoints(3,2)+bboxPoints(1,2))/2]; 
            bboxPointsOld=bboxPoints;

            bboxShiftDistanceThreshold=round(min(initialFrameSubregion(3)*fractionOfBoxJumpThreshold,initialFrameSubregion(4)*fractionOfBoxJumpThreshold));
            
            % Detect feature points in the face region.
            points = detectMinEigenFeatures(initFrame, 'ROI', round(initialFrameSubregion));

            %Display the detected points.
            %figure, imshow(videoFrame), hold on, title('Detected features');
            %plot(points);

            % Create a point tracker and enable the bidirectional error constraint to make it more robust in the presence of noise and clutter.
            pointTracker = vision.PointTracker('MaxBidirectionalError', 2);

            % Initialize the tracker with the initial point locations and the initial video frame.
            points = points.Location;
            initialize(pointTracker, points, initFrame);
            
            if plotTracking
                videoPlayer  = vision.VideoPlayer('Position',[100 100 [size(initFrame, 2), size(initFrame, 1)]+30]);
            end
            if saveTrackingVideo
               videoWriter = vision.VideoFileWriter(trackingFileName,'FrameRate',30);
            end
            %savePlottedTracking
            
            % Make a copy of the points to be used for computing the geometric transformation between the points in the previous and the current frames
            oldPoints = points;

            if saveFullOFMatrices %if to save all optic flow data
                allOF=zeros(numel(yInd),numel(xInd),nFrames,'like',complex(zeros(1,'single'),zeros(1,'single')));
                allIm=zeros(numel(yInd),numel(xInd),nFrames,'single');
            else
                allOF=[];
                allIm=[];
            end
            
            %% main loop
            pbboxUpdate=[];
            bboxCenterAll=zeros(nFrames,2);
            mOF=zeros(1,nFrames);
            skipBoundingBoxInSkip=round(skipFramesBoundingBox/skipFrames);
            parEyeTracking.skipBoundingBoxInSkip=skipBoundingBoxInSkip;

            hWB=waitbar(0,'Calculating optic flow');
            for i=1:nFrames
                %frame = step(videoReader); this is faster but cant start from an arbitrary frame or jump frames
                if nonConsecutiveVideo
                    videoFrame=rgb2gray(videoReader.read(pFrames(i))); %read will be replaced by readFrame in future versions but it is not possible to skip frames with readframes
                else
                    videoFrame = step(videoReader);
                    for j=1:numel(pFrames(i+1)-pFrames(i)-1)
                        step(videoReader);
                    end
                end
                
                if mod(i,skipBoundingBoxInSkip)==0
                    waitbar(i/nFrames,hWB);
                    
                    % Track the points. Note that some points may be lost.
                    [points, isFound] = step(pointTracker, videoFrame);
                    visiblePoints = points(isFound, :);
                    oldInliers = oldPoints(isFound, :);
                    
                    if size(visiblePoints, 1) >= 2 % need at least 2 points to ensure we are still reliably tracking the object
                        
                        % Estimate the geometric transformation between the old points and the new points and eliminate outliers
                        [xform, oldInliers, visiblePoints] = estimateGeometricTransform(oldInliers, visiblePoints, 'similarity', 'MaxDistance', 4);
                        
                        % Apply the transformation to the bounding box points
                        bboxPoints = transformPointsForward(xform, bboxPoints);
                        
                        % Reset the points
                        if size(oldInliers,1)<minTrackingPoints
                            newBox=round([min(bboxPoints(:,1)) min(bboxPoints(:,2))  max(bboxPoints(:,1))-min(bboxPoints(:,1)) max(bboxPoints(:,2))-min(bboxPoints(:,2))]);
                            newPoints = detectMinEigenFeatures(videoFrame, 'ROI', newBox );
                            newPoints = newPoints.Location;
                            in = inpolygon(newPoints(:,1),newPoints(:,2),bboxPoints(:,1),bboxPoints(:,2));
                            points=newPoints(in,:);
                            setPoints(pointTracker,points);
                            %initialize(pointTracker, points, initFrame);
                            oldPoints = points;
                        else
                            oldPoints = visiblePoints;
                            setPoints(pointTracker, oldPoints);
                        end
                        
                        bboxCenter=[(bboxPoints(3,1)+bboxPoints(1,1))/2 (bboxPoints(3,2)+bboxPoints(1,2))/2];
                        if sqrt((bboxCenter(1)-bboxCenterOld(1)).^2+(bboxCenter(2)-bboxCenterOld(2)).^2) > bboxShiftDistanceThreshold
                            %xInd=round(bboxPoints(1,1):bboxPoints(3,1));
                            %yInd=round(bboxPoints(1,2):bboxPoints(3,2));
                            
                            xyShift=round(bboxCenter-bboxCenterOld);
                            yInd=round(yInd-(xyShift(1)));
                            xInd=xInd-(xyShift(2));
                            if any(yInd<1)
                                yInd=1:numel(yInd);
                            end
                            if any(xInd<1)
                                xInd=1:numel(xInd);
                            end
                            if any(yInd>frameHeight)
                                yInd=(frameHeight-numel(yInd)+1):frameHeight;
                            end
                            if any(xInd>frameWidth)
                                xInd=(frameWidth-numel(xInd)+1):frameWidth;
                            end
                            
                            bboxPointsOld=bboxPoints;
                            bboxCenterOld=bboxCenter;
                            pbboxUpdate=[pbboxUpdate i];
                        end
                        
                        if plotTracking
                            
                            % Insert a bounding box around the object being tracked
                            bboxPolygon = reshape(bboxPoints', 1, []);
                            bboxPolygonOld = reshape(bboxPointsOld', 1, []);
                            
                            videoFramePlot = insertShape(videoFrame, 'Polygon', bboxPolygon,'LineWidth', 2);
                            videoFramePlot = insertShape(videoFramePlot, 'Polygon', bboxPolygonOld,'LineWidth', 2,'color','r');
                            
                            % Display tracked points
                            %videoFramePlot = insertMarker(videoFramePlot, visiblePoints, '+','Color', 'white');
                            
                            % Display the annotated video frame using the video player object
                            step(videoPlayer, videoFramePlot);
                            
                            if saveTrackingVideo %save tracked video
                                step(videoWriter, videoFramePlot);
                            end
                        end
                    else
                        disp(['Tracking analysis stopped at ' num2str(i) '/' num2str(nFrames) ' since all tracking points were lost']);
                        parEyeTracking.pStopDue2LostPoints=i;
                        mOF(i:end)=[];
                        bboxCenterAll(i:end,:)=[];
                        pFrames(i:end)=[];
                        break; %stop for loop
                    end
                    
                end
                im = step(converter, videoFrame(yInd,xInd));
                tmpOF=step(opticalFlow, im);
                
                if removeBorderOF
                    tmpOF(pBorder)=0;
                end
                
                if saveFullOFMatrices
                    allOF(:,:,i) = tmpOF;
                    allIm(:,:,i) = im;
                end
                
                mOF(i)=mean(mean(abs(tmpOF))); %mean velocity for every pixel
                bboxCenterAll(i,:)=bboxCenter;
                
            end
            close(hWB);
            
            save(obj.files.eyeTracking,'mOF','allOF','allIm','pbboxUpdate','parEyeTracking','pFrames','bboxCenterAll','initialFrameSubregion');
            
            % Clean uprelease(videoReader);
            release(opticalFlow);
            release(converter);
            release(pointTracker);
            if nonConsecutiveVideo
                delete(videoReader);
            else
                release(videoReader);
            end
            
            if saveTrackingVideo %save tracked video
                release(videoWriter);
            end
            if plotTracking
                release(videoPlayer);
            end
            
        end
        
        %% getDelta2BetaRatio
        function data=getDelta2BetaRatio(obj,varargin)
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'movLongWin',1000*60*30,@isnumeric); %max freq. to examine
            addParameter(parseObj,'movWin',10000,@isnumeric);
            addParameter(parseObj,'movOLWin',9000,@isnumeric);
            addParameter(parseObj,'segmentWelch',1000,@isnumeric);
            addParameter(parseObj,'dftPointsWelch',2^10,@isnumeric);
            addParameter(parseObj,'OLWelch',0.5);
            addParameter(parseObj,'tStart',0,@isnumeric);
            addParameter(parseObj,'win',0,@isnumeric); %if 0 uses the whole recording duration
            addParameter(parseObj,'deltaBandCutoff',4,@isnumeric);
            addParameter(parseObj,'betaBandLowCutoff',10,@isnumeric);
            addParameter(parseObj,'betaBandHighCutoff',40,@isnumeric);
            addParameter(parseObj,'applyNotch',0,@isnumeric);
            addParameter(parseObj,'maxVoltage',1500,@isnumeric);
            addParameter(parseObj,'overwrite',0,@isnumeric);
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
            
            %make parameter structure
            parDBRatio=parseObj.Results;
            
            if isnan(ch)
                disp('Error: no reference channel for Delta 2 Beta extraction');
                return;
            end
            %check if analysis was already done done
            obj.files.dbRatio=[obj.currentAnalysisFolder filesep 'dbRatio_ch' num2str(ch) '.mat'];
            if exist(obj.files.dbRatio,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.dbRatio);
                else
                    disp('DB analysis already exists for this recording');
                end
                return;
            end
            
            movWinSamples=movWin/1000*obj.filt.FFs;%obj.filt.FFs in Hz, movWin in samples
            movOLWinSamples=movOLWin/1000*obj.filt.FFs;
            timeBin=(movWin-movOLWin); %ms
            
            segmentWelchSamples = round(segmentWelch/1000*obj.filt.FFs);
            samplesOLWelch = round(segmentWelchSamples*OLWelch);
            
            %run welch once to get frequencies for every bin (f) determine frequency bands
            [~,f] = pwelch(randn(1,movWinSamples),segmentWelchSamples,samplesOLWelch,dftPointsWelch,obj.filt.FFs);
            pfLowBand=find(f<=deltaBandCutoff);
            pfHighBand=find(f>=betaBandLowCutoff & f<betaBandHighCutoff);
            
            %if obj.currentDataObj.recordingDuration_ms<movLongWin
            %    movLongWin=obj.currentDataObj.recordingDuration_ms;
            %end
            
            if win==0
                win=obj.currentDataObj.recordingDuration_ms-tStart;
                endTime=obj.currentDataObj.recordingDuration_ms;
            else
                endTime=min(win+tStart,obj.currentDataObj.recordingDuration_ms);
            end
            startTimes=tStart:(movLongWin-movOLWin):endTime;
            nChunks=numel(startTimes);
            deltaBetaRatioAll=cell(1,nChunks);
            t_ms=cell(1,nChunks);
            %deltaBetaRatioAllLow=cell(1,nChunks);;deltaBetaRatioAllHigh=cell(1,nChunks);
            fprintf('\nDelta2Beta extraction (%d chunks)-',nChunks);
            for i=1:nChunks
                fprintf('%d,',i);
                MLong=obj.currentDataObj.getData(ch,startTimes(i),movLongWin);
                if ~applyNotch
                    FMLong=obj.filt.F.getFilteredData(MLong);
                else
                    obj.filt.FN=filterData(Fs);
                    obj.filt.FN.filterDesign='cheby1';
                    obj.filt.FN.padding=true;
                    obj.filt.FN=obj.filt.FN.designNotch;
                    
                    FMLong=obj.filt.FN.getFilteredData(MLong); %for 50Hz noise
                    FMLong=obj.filt.F.getFilteredData(FMLong);
                end
                
                FMLong(FMLong<-maxVoltage | FMLong>maxVoltage)=nan; %remove high voltage movement artifacts
                
                FMLongB = buffer(FMLong,movWinSamples,movOLWinSamples,'nodelay');
                pValid=all(~isnan(FMLongB));
                
                [pxx,f] = pwelch(FMLongB(:,pValid),segmentWelchSamples,samplesOLWelch,dftPointsWelch,obj.filt.FFs);
                
                deltaBetaRatioAll{i}=zeros(1,numel(pValid));
                deltaBetaRatioAll{i}(pValid)=(mean(pxx(pfLowBand,:))./mean(pxx(pfHighBand,:)))';
                
                deltaRatioAll{i}=zeros(1,numel(pValid));
                deltaRatioAll{i}(pValid)=mean(pxx(pfLowBand,:))';
                
                betaRatioAll{i}=zeros(1,numel(pValid));
                betaRatioAll{i}(pValid)=mean(pxx(pfHighBand,:))';
                
                t_ms{i}=startTimes(i)+((movWin/2):timeBin:(movLongWin-movWin/2));
            end
            fprintf('\n');
            deltaBetaRatioAll{end}(t_ms{end}>(endTime-movWin/2))=NaN; 
            deltaRatioAll{end}(t_ms{end}>(endTime-movWin/2))=NaN; 
            betaRatioAll{end}(t_ms{end}>(endTime-movWin/2))=NaN; 
            
            bufferedBetaRatio=cell2mat(deltaBetaRatioAll);bufferedBetaRatio=bufferedBetaRatio(:);
            bufferedDeltaRatio=cell2mat(deltaRatioAll);bufferedDeltaRatio=bufferedDeltaRatio(:);
            bufferedOnlyBetaRatio=cell2mat(betaRatioAll);bufferedOnlyBetaRatio=bufferedOnlyBetaRatio(:);
            
            t_ms=cell2mat(t_ms);
            
            save(obj.files.dbRatio,'t_ms','bufferedBetaRatio','parDBRatio','bufferedOnlyBetaRatio','bufferedDeltaRatio');
        end
        
        %% getSlowCycles
        function data=getSlowCycles(obj,varargin)
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'medianFiltWin',1000*20,@isnumeric);
            addParameter(parseObj,'longOrdFiltWin',1000*1000,@isnumeric);
            addParameter(parseObj,'longOrdFiltOrd',0.6,@isnumeric);
            addParameter(parseObj,'maxCycle',1000*140,@isnumeric); %for excluding cycles which do not have a regular duration
            addParameter(parseObj,'minCycle',1000*5,@isnumeric); %for merging cycles with gaps less than these
            addParameter(parseObj,'overwrite',0,@isnumeric);
            
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
            
            %make parameter structure
            parSlowCycles=parseObj.Results;
            
            %check if analysis was already done done
            obj.files.slowCycles=[obj.currentAnalysisFolder filesep 'slowCycles_ch' num2str(ch) '.mat'];
            if exist(obj.files.slowCycles,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.slowCycles);
                else
                    disp('Slow cycle analysis already exists for this recording');
                end
                return;
            end
            
            dbRatioFile=[obj.currentAnalysisFolder filesep 'dbRatio_ch' num2str(ch) '.mat'];
            dbAutocorrFile=[obj.currentAnalysisFolder filesep 'dbAutocorr_ch' num2str(ch) '.mat'];
            
            obj.checkFileRecording(dbRatioFile,'Delta to beta analysis missing, please first run getDBRatio');
            load(dbRatioFile,'t_ms','bufferedBetaRatio','parDBRatio'); %load data  
            
            obj.checkFileRecording(dbAutocorrFile,'Delta to beta autocorr analysis missing, please first run getDBRatioAC');
            load(dbAutocorrFile,'pSleepDBRatio'); %load data
                
            timeBin=(parDBRatio.movWin-parDBRatio.movOLWin);
            
            %smooth with median filter
            medianFiltSamples=medianFiltWin/timeBin;
            DBRatioMedFilt = fastmedfilt1d(bufferedBetaRatio, medianFiltSamples);
            
            %long order filter to determine edges of DB fluctuation 
            longOrdFiltSamples=longOrdFiltWin/timeBin;
            longOrdFiltOrdSamples=longOrdFiltOrd*longOrdFiltSamples;
            DBLongOrdFilt = ordfilt2(DBRatioMedFilt, longOrdFiltOrdSamples, ones(longOrdFiltSamples,1));

            %plot(t_ms/1000/60/60,DBRatioMedFilt);hold on;plot(t_ms/1000/60/60,DBLongOrdFilt);
            
            edgeSamples=100;
            sortDBLongOrdFilt=sort(DBLongOrdFilt);
            sortDBLongOrdFilt(isnan(sortDBLongOrdFilt))=[];
            Th=mean(sortDBLongOrdFilt(1:edgeSamples))+(mean(sortDBLongOrdFilt((end-edgeSamples):end))-mean(sortDBLongOrdFilt(1:edgeSamples)))/2;
            %Th=min(DBLongOrdFilt)+(max(DBLongOrdFilt)-min(DBLongOrdFilt))/2;

            maxCycleSamples=maxCycle/timeBin;
            minCycleSamples=minCycle/timeBin;
            
            pTcycleOnset=find((DBRatioMedFilt(2:end)>=Th & DBRatioMedFilt(1:end-1)<Th) & pSleepDBRatio(1:end-1));
            pTcycleOnset(1+diff(pTcycleOnset)<minCycleSamples)=[];
            pTcycleOffset=pTcycleOnset(2:end);
            pTcycleOnset=pTcycleOnset(1:end-1);
            ppRemove=(pTcycleOffset-pTcycleOnset)>maxCycleSamples;
            
            pTcycleOffset(ppRemove)=[];
            pTcycleOnset(ppRemove)=[];
            
            TcycleOnset=t_ms(pTcycleOnset);
            TcycleOffset=t_ms(pTcycleOffset);
            
            edgesSamples=10;
            pTcycleMid=zeros(numel(TcycleOnset),1);
            for i=1:numel(TcycleOnset)
                pTmp=find(DBRatioMedFilt((pTcycleOnset(i)+edgesSamples):(pTcycleOffset(i)-edgesSamples))<Th,1,'first');
                if ~isempty(pTmp)
                    pTcycleMid(i)=pTmp;
                end
            end
            pTcycleMid=pTcycleMid+pTcycleOnset+edgesSamples+1;
            TcycleMid=t_ms(pTcycleMid);
            
            save(obj.files.slowCycles,'parSlowCycles','TcycleOnset','TcycleOffset','TcycleMid','pSleepDBRatio','DBRatioMedFilt');
        end
        
        %% plotDelta2BetaRatio
        function [h]=plotSlowCycles(obj,varargin)
            
            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'saveFigures',1,@isnumeric);
            addParameter(parseObj,'h',0,@ishandle);
            addParameter(parseObj,'printLocalCopy',0,@isnumeric);
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
            
            slowCyclesFile=[obj.currentAnalysisFolder filesep 'slowCycles_ch' num2str(ch) '.mat'];
            dbRatioFile=[obj.currentAnalysisFolder filesep 'dbRatio_ch' num2str(ch) '.mat'];
            
            obj.checkFileRecording(slowCyclesFile,'Delta to beta analysis missing, please first run getSlowCycles');
            obj.checkFileRecording(dbRatioFile,'Delta to beta analysis missing, please first run getDBRatio');
            load(slowCyclesFile); %load data
            load(dbRatioFile); %load data

            
            if h==0
                f=figure('Position',[200 200 900 500]);
                h(1)=subaxis(7,1,1,'S',0.01);
                h(2)=subaxis(7,1,2,'S',0.01);
                h(3)=subaxis(7,1,1,3,1,5,'S',0.01);
            else
                saveFigures=0;
            end
            
            %'parSlowCycles','TcycleOnset','TcycleOffset','TcycleMid','pSleep'
            
            axes(h(1));
            plot(t_ms(pSleepDBRatio')/1000/60/60,ones(1,numel(find(pSleepDBRatio))),'.k','MarkerSize',10);hold on;
            ylim([0 1]);
            axis off;
            l=legend('sleep');
            l.Box='off';l.Location='northeastoutside';l.Position=[0.8881    0.8805    0.1015    0.0979];
            
            axes(h(2));
            plot(TcycleOnset/1000/60/60,ones(1,numel(TcycleOnset)),'.b','MarkerSize',10);ylim([-1 2]);hold on;
            plot(TcycleOffset/1000/60/60,ones(1,numel(TcycleOffset)),'or','MarkerSize',10);ylim([-1 2]);
            plot(TcycleMid/1000/60/60,ones(1,numel(TcycleOffset)),'.g','MarkerSize',10);ylim([-1 2]);
            l=legend('onset','offset','middle');
            ylim([0.5 1.5]);
            axis off;
            l.Box='off';l.Location='northeastoutside';l.Position=[0.8881    0.6787    0.1015    0.0979];
            
            axes(h(3));
            plot(t_ms/1000/60/60,bufferedBetaRatio);
            ylabel('\delta/\beta ratio');
            xlabel('Time [h]');
            
            linkaxes(h,'x');
            
            if saveFigures
                set(f,'PaperPositionMode','auto');
                fileName=[obj.currentPlotFolder filesep 'slowCycles_ch' num2str(ch)];
                print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                if printLocalCopy
                    fileName=[cd filesep obj.par.Animal{obj.currentPRec} '_Rec' num2str(obj.currentPRec) '_slowCycles_ch' num2str(ch)];
                    print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                    print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                end
            end
        end
        
        %% getSharpWaves
        function [data]=getSharpWavesAnalysis(obj,varargin)
            
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'preSW',1000,@isnumeric);
            addParameter(parseObj,'winSW',2500,@isnumeric);
            addParameter(parseObj,'overwrite',0,@isnumeric);
                    
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
            
            %make parameter structure
            parSharpWavesAnalysis=parseObj.Results;

            %check if analysis was already done
            obj.files.sharpWaveAnalysis=[obj.currentAnalysisFolder filesep 'sharpWavesAnalysis_ch' num2str(ch) '.mat'];
            if exist(obj.files.sharpWaveAnalysis,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.sharpWaveAnalysis);
                else
                    disp('Sharp wave analysis already exists for this recording');
                end
                return;
            end
            
            sharpWavesFile=[obj.currentAnalysisFolder filesep 'sharpWaves_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(sharpWavesFile,'Sharp wave file missing, please run getSharpWaves');
            load(sharpWavesFile);
            
            %determine new length based on downsampling factor
            winSW=ceil(winSW/obj.filt.DS4Hz.downSamplingFactor)*obj.filt.DS4Hz.downSamplingFactor;
            
            nSW=min(1000,numel(tSW));
            [allRaw,tRaw]=obj.currentDataObj.getData(ch,tSW(1:nSW)-preSW,winSW);
            [allDS4Hz,tDS4Hz]=obj.filt.DS4Hz.getFilteredData(allRaw);
            [allFHR,tFHR]=obj.filt.FHR.getFilteredData(allRaw);
            
            allSWAbsAI=squeeze(mean(abs(reshape(permute(allFHR,[3,1,2]),[obj.filt.DS4Hz.downSamplingFactor  size(allFHR,3)/obj.filt.DS4Hz.downSamplingFactor nSW])),1))';
            meanProfiles.mSWAbsAI=mean(allSWAbsAI);
            meanProfiles.mSWAbsHP=squeeze(mean(abs(allFHR),2));
            meanProfiles.mSWRaw=squeeze(mean(allRaw,2));
            meanProfiles.mSWLP=squeeze(mean(allDS4Hz,2));
            meanProfiles.tSWAbsAI=tDS4Hz;
            meanProfiles.tSWAbsHP=tFHR;
            meanProfiles.tSWRaw=tRaw;
            meanProfiles.tSWLP=tDS4Hz;
            
            %polar((1:numel(mSWLP))/numel(mSWLP)*2*pi,mSWLP'/max(abs(mSWLP)),'b');hold on;
            %polar((1:numel(mSWAbsAI))/numel(mSWAbsAI)*2*pi,mSWAbsAI/max(mSWAbsAI),'r');
            
            %{
            FAall=zeros(size(FAtmp));
            FAallLog=zeros(size(FAtmp));
            for i=1:500
                fprintf('%d ',i);
                [~,~,~,FAtmp]=spectrogram(squeeze(allSW(1,i,:))',2^13,round(0.9*(2^13)),2^13,obj.filt.FHR.samplingFrequency);
                FAall=FAtmp/500+FAall;
                FAallLog=10*log10(FAtmp)/500+FAallLog;
            end
            imagesc(10*log10(FAall));
            imagesc(flipud(FAallLog));
            %}
            save(obj.files.sharpWaveAnalysis,'allSWAbsAI','meanProfiles','parSharpWavesAnalysis');
        end
        
        
        %% plotDelta2BetaRatio
        function [h]=plotSharpWavesAnalysis(obj,varargin)
            
            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'saveFigures',1,@isnumeric);
            addParameter(parseObj,'h',0,@ishandle);
            addParameter(parseObj,'printLocalCopy',0,@isnumeric);
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
            
            sharpWaveAnalysisFile=[obj.currentAnalysisFolder filesep 'sharpWavesAnalysis_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(sharpWaveAnalysisFile,'Delta to beta analysis missing, please first run getSharpWaveAnalysis');
            load(sharpWaveAnalysisFile); %load data
            
            if h==0
                f=figure('Position',[100 100 500 900]);
                h(1)=subaxis(f,2,1,1,'S',0.01,'ML',0.15);
                h(2)=subaxis(f,2,1,2,'S',0.01,'ML',0.15);
            else
                saveFigures=0;
            end
            
            plot(meanProfiles.tSWAbsHP,meanProfiles.mSWAbsHP,'Parent',h(1));
            ylabel('High pass amp','Parent',h(1));
            set(h(1),'XTickLabel',[]);
            plot(meanProfiles.tSWLP,meanProfiles.mSWLP,'r','Parent',h(2));
            ylabel('LFP [\muV]','Parent',h(2));
            xlabel('Time [ms]','Parent',h(2));
            
            linkaxes(h,'x');
            axis tight;

            if saveFigures
                set(f,'PaperPositionMode','auto');
                fileName=[obj.currentPlotFolder filesep 'sharpWaveAnalysis_ch' num2str(ch)];
                print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                if printLocalCopy
                    fileName=[cd filesep obj.par.Animal{obj.currentPRec} '_Rec' num2str(obj.currentPRec) '_sharpWaveAnalysis_ch' num2str(ch)];
                    print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                    print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                end
            end
        end
        
        
        %% getSharpWaves
        function data=getSharpWaves(obj,varargin)
            
            obj.checkFileRecording;

            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'nTestSegments',20,@isnumeric);
            addParameter(parseObj,'minPeakWidth',200,@isnumeric);
            addParameter(parseObj,'minPeakInterval',1000,@isnumeric);
            addParameter(parseObj,'detectOnlyDuringSWS',true);
            addParameter(parseObj,'preTemplate',400,@isnumeric);
            addParameter(parseObj,'winTemplate',1500,@isnumeric);
            addParameter(parseObj,'resultsFileName',[],@isstr);
            addParameter(parseObj,'percentile4ScaleEstimation',5,@isnumeric);
            addParameter(parseObj,'overwrite',0,@isnumeric);
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
            
            %make parameter structure
            parSharpWaves=parseObj.Results;
            
            %check if analysis was already done
            if isempty(resultsFileName)
                obj.files.sharpWaves=[obj.currentAnalysisFolder filesep 'sharpWaves_ch' num2str(ch) '.mat'];
            else
                obj.files.sharpWaves=[resultsFileName '.mat'];
            end
            if exist(obj.files.sharpWaves,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.sharpWaves);
                else
                    disp('Sharp wave analysis already exists for this recording');
                end
                return;
            end

            slowCyclesFile=[obj.currentAnalysisFolder filesep 'slowCycles_ch' num2str(ch) '.mat'];

            obj.checkFileRecording(slowCyclesFile,'slow cycle analysis missing, please first run getSlowCycles');
            load(slowCyclesFile); %load data
            
            nCycles=numel(TcycleOnset);
            pCycle=sort(randperm(nCycles,nTestSegments));
            
            Mtest=cell(nTestSegments,1);
            tTest=cell(nTestSegments,1);
            for i=1:numel(pCycle)
                MTmp=obj.currentDataObj.getData(ch,TcycleOnset(pCycle(i)),TcycleMid(pCycle(i))-TcycleOnset(pCycle(i)));
                [Mtest{i},tTest{i}]=obj.filt.DS4Hz.getFilteredData(MTmp);
                tTest{i}=tTest{i}'+TcycleOnset(pCycle(i));
                Mtest{i}=squeeze(Mtest{i});
            end
            Mtest=cell2mat(Mtest);
            tTest=cell2mat(tTest);

            sortedMtest=sort(Mtest);
            scaleEstimator=sortedMtest(round(percentile4ScaleEstimation/100*numel(sortedMtest)));
            tmpFs=obj.filt.DS4Hz.filteredSamplingFrequency;
            
            [peakVal,peakTime,peakWidth,peakProminance]=findpeaks(-Mtest,...
                'MinPeakHeight',-scaleEstimator,'MinPeakDistance',minPeakInterval/1000*tmpFs,'MinPeakProminence',-scaleEstimator/2,'MinPeakWidth',minPeakWidth/1000*tmpFs,'WidthReference','halfprom');
            
            [allSW,tSW]=obj.currentDataObj.getData(ch,tTest(peakTime)-preTemplate,winTemplate);
            [FLallSW,tFLallSW]=obj.filt.DS4Hz.getFilteredData(allSW);
            
            template=squeeze(median(FLallSW,2));
            nTemplate=numel(template);
            ccEdge=floor(nTemplate/2);
            [~,pTemplatePeak]=min(template);
            
            if detectOnlyDuringSWS
                TOn=TcycleOnset;
                TWin=TcycleMid-TcycleOnset;
            else
                seg=60000;
                TOn=0:seg:(obj.currentDataObj.recordingDuration_ms-seg);
                TWin=seg*ones(1,numel(TOn));
                nCycles=numel(TOn);
            end
            
            absolutePeakTimes=cell(nCycles,1);
            for i=1:nCycles
                [tmpM,tmpT]=obj.currentDataObj.getData(ch,TOn(i),TWin(i));
                [tmpFM,tmpFT]=obj.filt.DS4Hz.getFilteredData(tmpM);
                
                [C]=xcorrmat(squeeze(tmpFM),template);
                C=C(numel(tmpFM)-ccEdge:end-ccEdge);
                %C=xcorr(squeeze(tmpFM),template,'coeff');
                
                [~,peakTime]=findpeaks(C,'MinPeakHeight',0.1,'MinPeakProminence',0.2,'WidthReference','halfprom');
                peakTime(peakTime<=pTemplatePeak)=[]; %remove peaks at the edges where templates is not complete
                absolutePeakTimes{i}=tmpFT(peakTime-pTemplatePeak)'+TOn(i);
                
                %h(1)=subplot(2,1,1);plot(squeeze(tmpFM));h(2)=subplot(2,1,2);plot((1:numel(C))-pTemplatePeak,C);linkaxes(h,'x');
            end
            tSW=cell2mat(absolutePeakTimes);
            
            save(obj.files.sharpWaves,'tSW','parSharpWaves');
        end
        
        %% plotDelta2BetaRatio
        function [h]=plotDelta2BetaRatio(obj,varargin)
            
            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'saveFigures',1,@isnumeric);
            addParameter(parseObj,'chunksLength',1000*60*30,@isnumeric);
            addParameter(parseObj,'h',0,@ishandle);
            addParameter(parseObj,'printLocalCopy',0,@isnumeric);
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
            
            dbRatioFile=[obj.currentAnalysisFolder filesep 'dbRatio_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(dbRatioFile,'Delta to beta analysis missing, please first run getDelta2BetaRatio');
            load(dbRatioFile); %load data
            
            timeBin=(parDBRatio.movWin-parDBRatio.movOLWin);
            nSamples=numel(bufferedBetaRatio);
            tmov=(1:nSamples)*timeBin;
            
            movWinSamples=round(chunksLength/timeBin);
            chunks=buffer(bufferedBetaRatio,movWinSamples);
            tLong=t_ms(round(movWinSamples/2):movWinSamples:nSamples)/1000/60/60;
            
            sortedBetaRatio=sort(bufferedBetaRatio);
            estimateColorMapMax=round(sortedBetaRatio(round(numel(sortedBetaRatio)*0.95))/100)*100;
            
            if h==0
                fDB=figure('Position',[100 100 900 500]);
                h=axes;
            else
                saveFigures=0;
                axes(h);
            end
            imagesc((1:size(chunks,1))*timeBin/1000/60,tLong,chunks',[0 estimateColorMapMax]);
            xlabel('Time [min]');ylabel('Time [hour]');
            
            h(2)=colorbar;
            %set(cb,'position',[0.9115    0.7820    0.0096    0.1440]);
            ylabel(h(2),'\delta/\beta');
            
            if saveFigures
                set(fDB,'PaperPositionMode','auto');
                fileName=[obj.currentPlotFolder filesep 'dbRatio_ch' num2str(ch)];
                print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                if printLocalCopy
                    fileName=[cd filesep obj.par.Animal{obj.currentPRec} '_Rec' num2str(obj.currentPRec) '_dbRatio_ch' num2str(ch)];
                    print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                    print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                end
            end
        end
        
        %% plotDelta2BetaAC
        function h=plotDelta2BetaAC(obj,varargin)
            %sleepAnalysis.getDelta2BetaAC - input parameters: 
            parseObj = inputParser;
            parseObj.FunctionName='sleepAnalysis\plotDelta2BetaAC';
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'saveFigures',1,@isnumeric);
            addParameter(parseObj,'printLocalCopy',0,@isnumeric);
            addParameter(parseObj,'h',0,@ishandle);
            
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
            
            dbAutocorrFile=[obj.currentAnalysisFolder filesep 'dbAutocorr_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(dbAutocorrFile,'Autocorr analysis missing, please first run getDelta2BetaAC');
            load(dbAutocorrFile);
            
            if h==0
                fAC=figure;
                h=axes;
            else
                saveFigures=0;
                axes(h);
            end
            
            lineHandles = stem(xcf_lags/1000,real(xcf),'filled','r-o');
            set(lineHandles(1),'MarkerSize',4);
            grid('on');
            xlabel('Period [s]');
            ylabel('Auto corr.');
            hold('on');
            
            plot(period/1000,real(xcf(pPeriod)),'o','MarkerSize',5,'color','k');
            
            a = axis;
            plot([a(1) a(1); a(2) a(2)],[xcf_bounds([1 1]) xcf_bounds([2 2])],'-b');
            plot([a(1) a(2)],[0 0],'-k');
            hold('off');
            
            if saveFigures
                set(fAC,'PaperPositionMode','auto');
                fileName=[obj.currentPlotFolder filesep 'dbAC_ch' num2str(parDbAutocorr.ch) '_t' num2str(parDbAutocorr.tStart) '_w' num2str(parDbAutocorr.win)];
                print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                if printLocalCopy
                    fileName=[cd filesep obj.par.Animal{obj.currentPRec} '_Rec' num2str(obj.currentPRec) '_dbAC_ch' num2str(parDbAutocorr.ch) '_t' num2str(parDbAutocorr.tStart) '_w' num2str(parDbAutocorr.win)];
                    print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                    print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                end
            end
            
        end
      
        %% plotDelta2BetaSlidingAC
        function h=plotDelta2BetaSlidingAC(obj,varargin)
            %sleepAnalysis.plotDelta2BetaSlidingAC - input parameters: 
            parseObj = inputParser;
            parseObj.FunctionName='sleepAnalysis\plotDelta2BetaSlidingAC';
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'saveFigures',1,@isnumeric);
            addParameter(parseObj,'printLocalCopy',0,@isnumeric);
            addParameter(parseObj,'h',0);
            
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
            
            dbAutocorrFile=[obj.currentAnalysisFolder filesep 'dbAutocorr_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(dbAutocorrFile,'Autocorr analysis missing, please first run getDelta2BetaAC');
            load(dbAutocorrFile);

            if h(1)==0
                fSAC=figure('position',[200 200 550 600]);
                h(1)=subaxis(fSAC,2,1,1,'S',0.05,'M',0.1);
                h(2)=subaxis(fSAC,2,1,2,'S',0.05,'M',0.1);
            else
                saveFigures=0;
            end
            
            axes(h(1));
            h(3)=imagesc(tSlidingAC/1000/60/60,autocorrTimes/1000,real(acf),[-0.5 0.5]);
            ylabel('Autocorr lag [s]');
            ylim(xcf_lags([1 end])/1000);%important for panel plots
            yl=ylim;
            xlim(tSlidingAC([1 end])/1000/60/60); %important for panel plots
            xl=xlim;
            set(h(1),'YDir','normal');
            set(h(1),'XTickLabel',[]);
            hold on;
            
            x=[tStartSleep/1000/60/60 tEndSleep/1000/60/60 tEndSleep/1000/60/60 tStartSleep/1000/60/60];
            W=0.03;
            y=yl(2)+W*[diff(yl) diff(yl) diff(yl)*3 diff(yl)*3];
            h(4)=patch(x,y,[0.2 0.2 0.2],'Clipping','off','lineStyle','none','FaceAlpha',0.5); 
            text((x(1)+x(2))/2,(y(1)+y(3))/2,'E-Sleep','HorizontalAlignment','center','VerticalAlignment','middle');
            
            axes(h(2));
            
            h(5)=scatter(tSlidingAC(pSleepSlidingAC)/1000/60/60,acfPeriodAll(pSleepSlidingAC)/1000,5,[0.8 0.8 1],'filled');hold on;
            h(6)=plot(tFilteredSlidingPeriod/1000/60/60,filteredSlidingPeriod/1000,'-','lineWidth',3);
            ylabel('Period [s]');
            xlabel('Time [h]');
            set(h(2),'Box','on');
            axis tight;
            xlim(xl);
            yl=ylim;
            marg=diff(yl)*0.02;
            ylim([yl(1)-marg,yl(2)+marg]);
            h(7)=line(xlim,[period/1000 period/1000],'color',[1 0.8 0.8]);
            h(8:9)=line([parDbAutocorr.tStart parDbAutocorr.tStart;parDbAutocorr.tStart+parDbAutocorr.win parDbAutocorr.tStart+parDbAutocorr.win]'/1000/60/60,[yl;yl]','color',[0.8 1 0.8]);


            if saveFigures
                set(fSAC,'PaperPositionMode','auto');
                fileName=[obj.currentPlotFolder filesep 'dbSAC_ch' num2str(parDbAutocorr.ch) '_t' num2str(parDbAutocorr.tStart) '_w' num2str(parDbAutocorr.win)];
                print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                if printLocalCopy
                    fileName=[cd filesep obj.par.Animal{obj.currentPRec} '_Rec' num2str(obj.currentPRec) '_dbSAC_ch' num2str(parDbAutocorr.ch) '_t' num2str(parDbAutocorr.tStart) '_w' num2str(parDbAutocorr.win)];
                    print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                    print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                    print(fileName,'-dpdf',['-r' num2str(obj.figResEPS)]);
                end
            end
        end
                
        %% getDelta2BetaAC
        function [data]=getDelta2BetaAC(obj,varargin)
            %sleepAnalysis.getDelta2BetaAC - input parameters: ch,tStart,win,movOLWin,XCFLag,movingAutoCorrWin,movingAutoCorrOL
            parseObj = inputParser;
            parseObj.FunctionName='sleepAnalysis\getDelta2BetaAC';
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'tStart',0,@isnumeric); %max freq. to examine
            addParameter(parseObj,'win',obj.currentDataObj.recordingDuration_ms,@isnumeric);
            addParameter(parseObj,'maxPeriodBand',20,@isnumeric);
            addParameter(parseObj,'movOLWin',4000,@isnumeric);
            addParameter(parseObj,'XCFLag',500000,@isnumeric);
            addParameter(parseObj,'movingAutoCorrWin',1000*1000,@isnumeric);
            addParameter(parseObj,'movingAutoCorrOL',900*1000,@isnumeric);
            addParameter(parseObj,'oscilDurationMovingWin',60*60*1000,@isnumeric);
            addParameter(parseObj,'smoothingDuration',1000*60*60,@isnumeric);
            addParameter(parseObj,'oscilDurationThresh',0.25,@isnumeric);
            addParameter(parseObj,'overwrite',0,@isnumeric);
            
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
            
            %make parameter structure
            parDbAutocorr=parseObj.Results;
            
            %check if analysis was already done done
            obj.files.dbAutocorr=[obj.currentAnalysisFolder filesep 'dbAutocorr_ch' num2str(ch) '.mat'];
            if exist(obj.files.dbAutocorr,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.dbAutocorr);
                else
                    disp('Autocorr analysis already exists for this recording');
                end
                return;
            end
            
            dbRatioFile=[obj.currentAnalysisFolder filesep 'dbRatio_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(dbRatioFile,'Delta to beta analysis missing, please first run getDelta2BetaRatio');
            load(dbRatioFile); %load data
            
            bufferedBetaRatio(isnan(bufferedBetaRatio))=0; %for cross-corr analysis nans result in imaginary values
            
            %cross correlation analysis
            pOscillation=find(t_ms>tStart & t_ms<(tStart+win));
            timeBin=(parDBRatio.movWin-parDBRatio.movOLWin);
            XCFLagSamples=XCFLag/timeBin;
            [xcf,xcf_lags,xcf_bounds]=crosscorr(bufferedBetaRatio(pOscillation),bufferedBetaRatio(pOscillation),XCFLagSamples);
            
            xcf_lags=xcf_lags*1000;
            %calculate periodicity
            
            %find first vally and peak in the autocorrelation function
            [~,pPeak] = findpeaks(xcf(XCFLagSamples+1:end),'MinPeakProminence',0.1);
            [~,pVally] = findpeaks(-xcf(XCFLagSamples+1:end),'MinPeakProminence',0.1);
            if isempty(pPeak) %if peak is weak, try a different value
                [~,pPeak] = findpeaks(xcf(XCFLagSamples+1:end),'MinPeakProminence',0.05);
                [~,pVally] = findpeaks(-xcf(XCFLagSamples+1:end),'MinPeakProminence',0.05);
                disp('Prominance for peak detection was reduced to 0.05 due to low periodicity values!!!');
            end
            
            
            
            if isempty(pPeak) | isempty(pVally)
                pPeriod=NaN;
                period=NaN;
                pVally=NaN;
                vallyPeriod=NaN;
                peak2VallyDiff=NaN;
            else
                pPeriod=pPeak(1)+XCFLagSamples;
                period=xcf_lags(pPeriod);
                pAntiPeriod=pVally(1)+XCFLagSamples;
                vallyPeriod=xcf_lags(pAntiPeriod);
                peak2VallyDiff=xcf(pPeriod)-xcf(pAntiPeriod);
            end
            
            
            %sliding autocorr analysis
            movingAutoCorrWinSamples=movingAutoCorrWin/timeBin;
            movingAutoCorrOLSamples=movingAutoCorrOL/timeBin;
            autoCorrTimeBin=(movingAutoCorrWin-movingAutoCorrOL);
            BetaRatioForSlidingAutocorr = buffer(bufferedBetaRatio,movingAutoCorrWinSamples,movingAutoCorrOLSamples,'nodelay');
            tSlidingAC=(movingAutoCorrWin/2):(movingAutoCorrWin-movingAutoCorrOL):(t_ms(end)-movingAutoCorrWin/2+movingAutoCorrWin-movingAutoCorrOL);
           
            %R=xcorrmat(BetaRatioForSlidingAutocorr,BetaRatioForSlidingAutocorr,autoCorrSamples);
            
            acfSamples=floor(movingAutoCorrWinSamples/2);
            acf=zeros(size(BetaRatioForSlidingAutocorr,1)+1,size(BetaRatioForSlidingAutocorr,2));
            peak2VallyDiff=zeros(1,size(BetaRatioForSlidingAutocorr,2));
            for i=1:size(BetaRatioForSlidingAutocorr,2)
                [acf(:,i),autoCorrSamples] = crosscorr(BetaRatioForSlidingAutocorr(:,i),BetaRatioForSlidingAutocorr(:,i),acfSamples);
                %calculate peak2VallyDiff for different times
                acf(:,i)=smooth(acf(:,i),10,'moving');
                
                [acfPeakAll(i),acfPeriodAll(i)]=max(acf((acfSamples+pPeak(1)-maxPeriodBand):(acfSamples+pPeak(1)+maxPeriodBand),i));
                [acfVallyAll(i),acfAntiPeriodAll(i)]=min(acf((acfSamples+pVally(1)-maxPeriodBand):(acfSamples+pVally(1)+maxPeriodBand),i));
                peak2VallyDiff(i)=acfPeakAll(i)-acfVallyAll(i);
                
                %{
                [~,pPeak] = findpeaks(acf(acfSamples+1:end,i),'MinPeakProminence',0.1);
                [~,pVally] = findpeaks(-acf(acfSamples+1:end,i),'MinPeakProminence',0.1);
                findpeaks(acf(acfSamples+1:end,i),'MinPeakProminence',0.1);
                disp([pPeak(1) pVally(1)])
                pause;
                
                if ~isempty(pPeak) & ~isempty(pVally)
                    pPeriodTmp=pPeak(1)+acfSamples;
                    pVally=pVally(1)+acfSamples;
                    peak2VallyDiff(i)=acf(pPeriodTmp,i)-acf(pVally,i);
                end
                %}
            end
            autocorrTimes=autoCorrSamples*timeBin;
            acfPeriodAll=autocorrTimes((acfPeriodAll+acfSamples+pPeak(1)-maxPeriodBand-1));
            
            oscilDurationMovingSamples=oscilDurationMovingWin/autoCorrTimeBin;
            tmpOscDuration=peak2VallyDiff>oscilDurationThresh;
            filtOscilDuration = medfilt1(double(tmpOscDuration),oscilDurationMovingSamples);
            pSleepSlidingAC=filtOscilDuration>=0.5;
            
            tmpBin=movingAutoCorrWinSamples-movingAutoCorrOLSamples;
            pSleepDBRatio=false(numel(bufferedBetaRatio),1);
            for i=1:numel(pSleepSlidingAC)
                pSleepDBRatio(((i-1)*tmpBin+1):(i*tmpBin))=pSleepSlidingAC(i);
            end
            
            pStartSleep=find(pSleepDBRatio==1,1,'first');
            tStartSleep=t_ms(pStartSleep);
            tEndSleep=t_ms(find(pSleepDBRatio(pStartSleep:end)==0,1,'first')+pStartSleep);
            
            pSleepSlidingAC=find(tSlidingAC>=tStartSleep & tSlidingAC<=tEndSleep & peak2VallyDiff>oscilDurationThresh);
            
            smoothingSamples=round(smoothingDuration/autoCorrTimeBin);
            filteredSlidingPeriod=smooth(tSlidingAC(pSleepSlidingAC),acfPeriodAll(pSleepSlidingAC),smoothingSamples,'moving');
            edgeSamples=tSlidingAC(pSleepSlidingAC)<=(tSlidingAC(pSleepSlidingAC(1))+smoothingDuration/2) | tSlidingAC(pSleepSlidingAC)>=(tSlidingAC(pSleepSlidingAC(end))-smoothingDuration/2);
            filteredSlidingPeriod(edgeSamples)=[];
            tFilteredSlidingPeriod=tSlidingAC(pSleepSlidingAC(~edgeSamples))';
            
            %save data
            save(obj.files.dbAutocorr,'parDbAutocorr','xcf','xcf_lags','xcf_bounds','BetaRatioForSlidingAutocorr','autoCorrTimeBin','autocorrTimes','timeBin',...
                'pPeriod','period','acf','vallyPeriod','peak2VallyDiff','pSleepDBRatio','pSleepSlidingAC','acfPeakAll','acfVallyAll','tSlidingAC','acfPeriodAll',...
                'tStartSleep','tEndSleep','filteredSlidingPeriod','tFilteredSlidingPeriod','pSleepSlidingAC');
        end
        
        %% plotFreqBandDetection
        function [h]=plotFreqBandDetection(obj,varargin)
            
            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'plotDendrogram',true);
            addParameter(parseObj,'plotSpectralBands',true);
            addParameter(parseObj,'savePlots',true);
            addParameter(parseObj,'freqBandFile',[]);
            addParameter(parseObj,'cLim',0);
            addParameter(parseObj,'printLocalCopy',0,@isnumeric);
            addParameter(parseObj,'hDendro',0);
            addParameter(parseObj,'hSpectra',0);
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
            
            if isempty(freqBandFile)
                spectralClusteringFile=[obj.currentAnalysisFolder filesep 'spectalClustering_ch' num2str(ch) '.mat'];
            else
                spectralClusteringFile=freqBandFile;
            end
            obj.checkFileRecording(spectralClusteringFile,'Spectral band analysis missing, please first run getDelta2BetaRatio');
            load(spectralClusteringFile);
            
            if plotDendrogram
                maxDendroClusters=2;
                
                if cLim==0
                    cLim=[];
                end
                if hDendro==0
                    hDendro=[];
                else
                    savePlots=[];
                end
                [DC,order,clusters,h]=DendrogramMatrix(corrMat,'linkMetric','euclidean','linkMethod','ward','maxClusters',maxDendroClusters,...
                    'toPlotBinaryTree',1,'cLim',cLim,'hDendro',hDendro,'plotOrderLabels',0);
                %h(3).Position=[0.9149    0.7595    0.0137    0.1667];
                ylabel(h(3),'Corr.');
                xlabel(h(2),'Segment');
                xlabel(h(1),'Distance');
                if savePlots
                    set(gcf,'PaperPositionMode','auto');
                    fileName=[obj.currentPlotFolder filesep 'dendrogram_ch' num2str(parFreqBandDetection.ch) '_t' num2str(parFreqBandDetection.tStart) '_w' num2str(parFreqBandDetection.win)];
                    print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                    print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                    if printLocalCopy
                        fileName=[cd filesep obj.par.Animal{obj.currentPRec} '_Rec' num2str(obj.currentPRec) '_dendrogram_ch' num2str(parFreqBandDetection.ch) '_t' num2str(parFreqBandDetection.tStart) '_w' num2str(parFreqBandDetection.win)];
                        print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                        print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                    end
                end
            end
            
            if plotSpectralBands
                if hSpectra==0
                    fTmp=figure('position',[680   678   658   420]);
                    hTmp=axes;
                    h=[h hTmp];
                else
                    axes(hSpectra);
                    h=[h hSpectra];
                    savePlots=0;
                end
                for i=1:maxDendroClusters
                    PS=mean(normsPxx(:,clusters==i),2);
                    plot(freqHz,PS,'lineWidth',2);hold on;
                end
                plot(crossFreq,PS(crossFreq==freqHz),'ok','MarkerSize',8,'LineWidth',2);
                text(crossFreq+(diff(xlim))*0.15,PS(crossFreq==freqHz),'F_{trans.}');
                xlabel('Frequency (Hz)');
                ylabel('nPSD');
                xlim([0 parFreqBandDetection.fMax]);
                
                if savePlots
                    set(fTmp,'PaperPositionMode','auto');
                    fileName=[obj.currentPlotFolder filesep 'spectralBands_ch' num2str(parFreqBandDetection.ch) '_t' num2str(parFreqBandDetection.tStart) '_w' num2str(parFreqBandDetection.win)];
                    print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                    print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                    if printLocalCopy
                        fileName=[cd filesep obj.par.Animal{obj.currentPRec} '_Rec' num2str(obj.currentPRec) '_spectralBands_ch' num2str(parFreqBandDetection.ch) '_t' num2str(parFreqBandDetection.tStart) '_w' num2str(parFreqBandDetection.win)];
                        print(fileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
                        print(fileName,'-depsc',['-r' num2str(obj.figResEPS)]);
                    end
                end
                
            end
            
        end
        
        
        %% getHPSegments 
        function data=getAwakeVsSleepFreq(obj,varargin)
            %% parameter and settings
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'sleepFreqBandFile',[]); %median filter window for extracting optic flow baseline
            addParameter(parseObj,'binDuration',10000);
            addParameter(parseObj,'fMax',30,@isnumeric); %max freq. to examine
            addParameter(parseObj,'maxWin',1000*60*60*2,@isnumeric);
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'saveFileName',[]);
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
            
            %make parameter structure
            parAwakeVsSleepFreq=parseObj.Results;
            
            if isempty(sleepFreqBandFile)
                spectralClusteringFile=[obj.currentAnalysisFolder filesep 'spectalClustering_ch' num2str(ch) '.mat'];
            else
                spectralClusteringFile=[sleepFreqBandFile '.mat'];
            end
            %obj.checkFileRecording(spectralClusteringFile,'Spectral band analysis missing, please first run getFreqBandDetection');
            %dataSleep=load(spectralClusteringFile,'sPxx','freqHz','normsPxx','clusters');
            
            %check if analysis was already done done
            if isempty(saveFileName)
                obj.files.AwakeVsSleepFreq=[obj.currentAnalysisFolder filesep 'AwakeVsSleepFreq_ch' num2str(ch)];
            else
                obj.files.AwakeVsSleepFreq=saveFileName;
            end
            
            if exist(obj.files.AwakeVsSleepFreq,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.AwakeVsSleepFreq);
                else
                    disp(['AwakeVsSleepFreq analysis file already exists']);
                end
                return;
            end
            
            animalStates=strsplit(obj.par.AnimalState{obj.currentPRec},'/');
            awakeStartTimeSec=obj.par.tStartAwake{obj.currentPRec};
            
            for i=1:numel(animalStates)
                if strcmp(animalStates{i},'Awake') || strcmp(animalStates{i},'Running') || strcmp(animalStates{i},'Resting')
                    recDuration=obj.currentDataObj.recordingDuration_ms;
                    if ~isnan(awakeStartTimeSec)
                        tStart=awakeStartTimeSec*1000;
                    else
                        tStart=0;
                    end
                    win=min(maxWin,floor((recDuration-tStart)/binDuration)*binDuration);
                    obj.getFreqBandDetection('tStart',tStart,'win',win,'binDuration',binDuration,'saveFile',obj.files.AwakeVsSleepFreq,'maxDendroClusters',1,'overwrite',overwrite,'fMax',fMax);
                    dataAwake=obj.getFreqBandDetection('saveFile',obj.files.AwakeVsSleepFreq);
                    
                    %{
                    for j=unique(dataSleep.clusters)'
                        PS(j,:)=mean(10*log10(dataSleep.normsPxx(:,dataSleep.clusters==j)),2);
                    end
                    PA=mean(bsxfun(@rdivide,dataAwake.sPxx,mean(dataSleep.sPxx,2)),2);
                    plot(dataSleep.freqHz,PS,'lineWidth',2);hold on;
                    plot(dataAwake.freqHz,PA,'lineWidth',2);
                    %}
                end
            end
            
            %save(obj.files.AwakeVsSleepFreq,'parAwakeVsSleepFreq');
        end
        
        %% getActivity4OpenVsClosedEyes
        function getActivity4OpenVsClosedEyes(obj,varargin)
            obj.checkFileRecording;

            parseObj = inputParser;
            addParameter(parseObj,'ch',[],@isnumeric);
            addParameter(parseObj,'saveFile',[]);
            addParameter(parseObj,'nZoomPanels',1);
            addParameter(parseObj,'overwrite',0,@isnumeric);
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
            
            %make parameter structure
            parActivity4OpenVsClosedEyes=parseObj.Results;
            
            %check if analysis was already done done
            if isempty(saveFile)
                obj.files.Activity4OpenVsClosedEyes=[obj.currentAnalysisFolder filesep 'Activity4OpenVsClosedEyes_ch' num2str(ch) '.mat'];
            else
                obj.files.Activity4OpenVsClosedEyes=[saveFile '.mat'];
            end
            
            if exist(obj.files.Activity4OpenVsClosedEyes,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.Activity4OpenVsClosedEyes);
                else
                    disp('Activity4OpenVsClosedEyes analysis already performed, use ''overwrite'' to recalculate');
                end
                return;
            end
            
            obj=obj.getFileNames;
            if exist(obj.files.dbRatio,'file')
                db=load(obj.files.dbRatio);
            end
            if exist(obj.files.dbRatio,'file')
                dbAC=load(obj.files.dbAutocorr,'tStartSleep','tEndSleep');
            end
            
            eye=load([obj.currentAnalysisFolder filesep 'eye.mat']);
            tmpEvents=min(numel(eye.eyeOpenStart),numel(eye.eyeCloseStart));
            if (eye.eyeCloseStart(1)-eye.eyeOpenStart(1)) >= 0 %first anotation segment is eye open
                openEyeStart=eye.eyeOpenStart(1:tmpEvents);
                openEyeEnd=eye.eyeCloseStart(1:tmpEvents);
                closedEyeStart=eye.eyeCloseStart(1:tmpEvents-1);
                closedEyeEnd=eye.eyeOpenStart(2:tmpEvents);
            else %first anotation segment is eyes close
                openEyeStart=eye.eyeOpenStart(1:tmpEvents-1);
                openEyeEnd=eye.eyeCloseStart(2:tmpEvents);
                closedEyeStart=eye.eyeCloseStart(1:tmpEvents);
                closedEyeEnd=eye.eyeOpenStart(1:tmpEvents);
            end
            
            pNonSleepOpen=find((openEyeStart<dbAC.tStartSleep & openEyeEnd<dbAC.tStartSleep) | (openEyeStart>dbAC.tEndSleep & openEyeEnd>dbAC.tEndSleep));
            pNonSleepClosed=find((closedEyeStart<dbAC.tStartSleep & closedEyeEnd<dbAC.tStartSleep) | (closedEyeStart>dbAC.tEndSleep & closedEyeEnd>dbAC.tEndSleep));
            
            %remove open closed segments that occur during the main sleep epoch
            openEyeStart=openEyeStart(pNonSleepOpen);
            openEyeEnd=openEyeEnd(pNonSleepOpen);
            closedEyeStart=closedEyeStart(pNonSleepClosed);
            closedEyeEnd=closedEyeEnd(pNonSleepClosed);
            
            allOpenDb={};allClosedDb={};
            for i=1:numel(openEyeStart)
                pTmp=find(db.t_ms>openEyeStart(i) & db.t_ms<=openEyeEnd(i));
                allOpenDb{i}=db.bufferedBetaRatio(pTmp);
            end
            for i=1:numel(closedEyeStart)
                pTmp=find(db.t_ms>closedEyeStart(i) & db.t_ms<=closedEyeEnd(i));
                allClosedDb{i}=db.bufferedBetaRatio(pTmp);
            end
            allOpenDb=cell2mat(allOpenDb');
            allClosedDb=cell2mat(allClosedDb');
            
            
            
            f=figure('position',[520   -97   560   895]);
            for i=1:nZoomPanels

                h(i)=subaxis(f,nZoomPanels+1,1,i,'s',0.03,'mt',0.01);
                plot(db.t_ms/1000/60,db.bufferedBetaRatio);hold on;
                yl=ylim;
                patch([openEyeStart;openEyeEnd;openEyeEnd;openEyeStart]/1000/60,(ones(numel(openEyeStart),1)*[yl(1) yl(1) yl(2) yl(2)])',[0 0 1],'FaceAlpha',0.2,'edgeColor','none');
                patch([closedEyeStart;closedEyeEnd;closedEyeEnd;closedEyeStart]/1000/60,(ones(numel(closedEyeStart),1)*[yl(1) yl(1) yl(2) yl(2)])',[1 0 0],'FaceAlpha',0.2,'edgeColor','none');
                axis tight;
                
                if i==1
                    [hl,ho]=legend({'\delta/\beta','open','closed'},'box','off','location','northwest');
                    %horizontalLegend(ho);
                end
            end
            
            xlabel('Time [min]');
            ylabel('\delta/\beta');
            
            
            
            h(nZoomPanels+1)=subaxis(f,nZoomPanels+1,1,nZoomPanels+1,'s',0.01);
            
            maxEdge=min(2000,6*std([allOpenDb;allClosedDb]));
            
            edges=[0:(maxEdge/10):maxEdge];
            [IOpen]=histc(allOpenDb,edges);
            [IClosed]=histc(allClosedDb,edges);
            bar(edges,[IOpen./sum(IOpen) IClosed/sum(IClosed)],1.1);
            xlabel('\delta/\beta');
            ylabel('Prob.');
            l=legend({'Open','Close'},'box','off');
        end
        
        %% getFreqBandDetectionEMG - function under construction
        function [data]=getFreqBandDetectionEMG(obj,varargin)
            obj.checkFileRecording;

            parseObj = inputParser;
            addParameter(parseObj,'ch',[],@isnumeric);
            addParameter(parseObj,'fMax',500,@isnumeric); %max freq. to examine
            addParameter(parseObj,'dftPoints',2^12,@isnumeric);
            addParameter(parseObj,'tStart',0,@isnumeric);
            addParameter(parseObj,'win',1000*60*60,@isnumeric);
            addParameter(parseObj,'maxDendroClusters',2,@isnumeric);
            addParameter(parseObj,'saveFile',[]);
            addParameter(parseObj,'overwrite',0,@isnumeric);
            addParameter(parseObj,'segmentLength',1000);
            addParameter(parseObj,'WelchOL',0.5);
            addParameter(parseObj,'binDuration',10000);
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
            
            %make parameter structure
            parFreqBandDetectionEMG=parseObj.Results;
            
            %check if analysis was already done done
            if isempty(saveFile)
                obj.files.spectralClusteringEMG=[obj.currentAnalysisFolder filesep 'spectalClusteringEMG_ch' num2str(ch) '.mat'];
            else
                obj.files.spectralClusteringEMG=[saveFile '.mat'];
            end
            
            if exist(obj.files.spectralClusteringEMG,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.spectralClusteringEMG);
                else
                    disp('Spectral clustering EMG analysis already exists for this recording');
                end
                return;
            end
            
            obj.filt.EMG1=filterData(obj.currentDataObj.samplingFrequency(1));
            obj.filt.EMG1.downSamplingFactor=32;
            obj.filt.EMG1.padding=true;
            obj.filt.EMG1=obj.filt.EMG1.designDownSample;
            
            MLong=obj.currentDataObj.getData(ch,tStart,win);
            FMLong=obj.filt.EMG1.getFilteredData(MLong);
            
            %calculate initial parameters
            segmentSamples = round(segmentLength/1000*obj.filt.EMG1.filteredSamplingFrequency);
            samplesOL = round(segmentSamples*WelchOL);
            samplesBin = binDuration/1000*obj.filt.EMG1.filteredSamplingFrequency;
            
            nBins=numel(FMLong)/samplesBin;

            FMLongB=reshape(FMLong,[samplesBin,nBins]);
            
            if (numel(FMLong)/samplesBin)~=round(numel(FMLong)/samplesBin)
                nBins=nBins-1;
                FMLong=FMLong(1:(samplesBin*nBins));
                disp('Last bin in recording not included due to a missmatch between recording duration and binDuration');
            end
                
            [pxx,f] = pwelch(FMLongB,segmentSamples,samplesOL,dftPoints,obj.filt.EMG1.filteredSamplingFrequency);
            %plot(10*log10(pxx))
            p=find(f<fMax);
            pp=find(sum(pxx(p,:))<500);
            
            %{
            [data]=obj.getDelta2BetaRatio;
            
            
            pDB=find(data.t_ms>tStart & data.t_ms<tStart+win);
            tEMG=2*1000*60*60+((5*1000):(10*1000):(2*60*60*1000));
            
            h(1)=subplot(2,1,1);
            plot(data.t_ms(pDB),data.bufferedBetaRatio(pDB)); hold on;
            h(2)=subplot(2,1,2);
            plot(tEMG,sum(pxx(p,:)),'r');
            plot(tEMG,mean(abs(FMLongB),1),'r');
            linkaxes(h,'x');
            
            figure;
            plot(data.t_ms(pDB),data.bufferedBetaRatio(pDB)); hold on;
            plot(tEMG(order(285:387)),200*ones(1,numel(order(285:387))),'.')
            %}
            
            sPxx=pxx(p,pp);
            freqHz=f(p);
            normsPxx=bsxfun(@rdivide,sPxx,mean(sPxx,2));
            corrMat=corrcoef(normsPxx);
            
            if maxDendroClusters==2
                
                [DC,order,clusters]=DendrogramMatrix(corrMat,'linkMetric','euclidean','linkMethod','ward','maxClusters',maxDendroClusters);
                
                S1=mean(normsPxx(:,clusters==1),2);
                S2=mean(normsPxx(:,clusters==2),2);
                if mean(S1(1:3))>mean(S2(1:3))
                    crossFreq=freqHz(1+find(S2-S1>=0,1,'last'));
                else
                    crossFreq=freqHz(1+find(S1-S2>=0,1,'last'));
                end
            else
                crossFreq=[];order=[];clusters=[];
            end
            
            save(obj.files.spectralClustering,'corrMat','sPxx','normsPxx','freqHz','parFreqBandDetectionEMG','order','clusters','crossFreq');
        end
        

        %% getFreqBandDetection
        function [data]=getFreqBandDetection(obj,varargin)
            obj.checkFileRecording;

            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.par.DVRLFPCh{obj.currentPRec},@isnumeric);
            addParameter(parseObj,'fMax',30,@isnumeric); %max freq. to examine
            addParameter(parseObj,'dftPoints',2^10,@isnumeric);
            addParameter(parseObj,'tStart',0,@isnumeric);
            addParameter(parseObj,'win',1000*60*60,@isnumeric);
            addParameter(parseObj,'maxDendroClusters',2,@isnumeric);
            addParameter(parseObj,'saveFile',[]);
            addParameter(parseObj,'overwrite',0,@isnumeric);
            addParameter(parseObj,'segmentLength',1000);
            addParameter(parseObj,'WelchOL',0.5);
            addParameter(parseObj,'binDuration',10000);
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
            
            %make parameter structure
            parFreqBandDetection=parseObj.Results;
            
            %check if analysis was already done done
            if isempty(saveFile)
                obj.files.spectralClustering=[obj.currentAnalysisFolder filesep 'spectalClustering_ch' num2str(ch) '.mat'];
            else
                obj.files.spectralClustering=[saveFile '.mat'];
            end
            
            if exist(obj.files.spectralClustering,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.spectralClustering);
                else
                    disp('Spectral clustering analysis already exists for this recording');
                end
                return;
            end
            
            MLong=obj.currentDataObj.getData(ch,tStart,win);
            FMLong=obj.filt.F.getFilteredData(MLong);
            
            %calculate initial parameters
            segmentSamples = round(segmentLength/1000*obj.filt.FFs);
            samplesOL = round(segmentSamples*WelchOL);
            samplesBin = binDuration/1000*obj.filt.FFs;
            
            nBins=numel(FMLong)/samplesBin;

            FMLongB=reshape(FMLong,[samplesBin,nBins]);
            
            if (numel(FMLong)/samplesBin)~=round(numel(FMLong)/samplesBin)
                nBins=nBins-1;
                FMLong=FMLong(1:(samplesBin*nBins));
                disp('Last bin in recording not included due to a missmatch between recording duration and binDuration');
            end
                
            [pxx,f] = pwelch(FMLongB,segmentSamples,samplesOL,dftPoints,obj.filt.FFs);
            %plot(10*log10(pxx))
            p=find(f<fMax);
            pp=find(sum(pxx(p,:))<0.4e6);
            
            sPxx=pxx(p,pp);
            freqHz=f(p);
            normsPxx=bsxfun(@rdivide,sPxx,mean(sPxx,2));
            corrMat=corrcoef(normsPxx);
            
            if maxDendroClusters==2
                
                [DC,order,clusters]=DendrogramMatrix(corrMat,'linkMetric','euclidean','linkMethod','ward','maxClusters',maxDendroClusters);
                
                S1=mean(normsPxx(:,clusters==1),2);
                S2=mean(normsPxx(:,clusters==2),2);
                if mean(S1(1:3))>mean(S2(1:3))
                    crossFreq=freqHz(find(S2-S1>=0,1,'first'));
                else
                    crossFreq=freqHz(find(S1-S2>=0,1,'first'));
                end
            else
                crossFreq=[];order=[];clusters=[];
            end
            
            save(obj.files.spectralClustering,'corrMat','sPxx','normsPxx','freqHz','parFreqBandDetection','order','clusters','crossFreq');
        end
        
        %% getFileNames
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
        
        %% batchProcessData
        function [outArgAll]=batchProcessData(obj,method,recNumbers,varargin)
            
            nOut=nargout;
            nRec=numel(recNumbers);
            
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
                    tmpObj=obj.setCurrentRecording(recNumbers(i));
                    fprintf('%d ',recNumbers(i));
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
                    tmpObj=obj.setCurrentRecording(recNumbers(i));
                    fprintf('%d ',recNumbers(i));
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
        
        %% getDigitalTriggers
        function data=getDigitalTriggers(obj,varargin)
            
            parseObj = inputParser;
            addParameter(parseObj,'saveFile',[]);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            addParameter(parseObj,'overwrite',0,@isnumeric);
            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end
            
            %evaluate all input parameters in workspace
            for i=1:numel(parseObj.Parameters)
                eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
            end
            
            obj.checkFileRecording;
            
            %check if analysis was already done done
            if isempty(saveFile)
                obj.files.digiTrig=[obj.currentAnalysisFolder filesep 'digiTrig.mat'];
            else
                obj.files.digiTrig=[saveFile '.mat'];
            end
            
            if exist(obj.files.digiTrig,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.digiTrig);
                else
                    disp('Triggers already exists for this recording');
                end
                return;
            end
            
            [tTrig,~,tTrig_string]=obj.currentDataObj.getTrigger;
            save(obj.files.digiTrig,'tTrig','tTrig_string');
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
        
        %% setCurrentRecording
        function [obj]=setCurrentRecording(obj,recName) %if recNumber is negative loads all recording related files but not recording object
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
                    elseif strcmp(allFullFiles{1}(end-3:end),'.kwd') %Intan recording
                        obj.currentDataObj=KwikRecording(obj.currentMEAFiles);
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
        
        %% getFilters
        function [obj]=getFilters(obj,Fs)
            if nargin==1
                if isempty(obj.currentDataObj)
                    error('Sampling frequency is required as an input');
                else
                    Fs=obj.currentDataObj.samplingFrequency(1);
                    disp(['sampling frequency set to that of current recording:' num2str(Fs) '[Hz]']);
                end
            end
            
            obj.filt.F=filterData(Fs);
            obj.filt.F.downSamplingFactor=Fs/250;
            obj.filt.F=obj.filt.F.designDownSample;
            obj.filt.F.padding=true;
            obj.filt.FFs=obj.filt.F.filteredSamplingFrequency;
            
            obj.filt.DS4Hz=filterData(Fs);
            obj.filt.DS4Hz.downSamplingFactor=Fs/250;
            obj.filt.DS4Hz.lowPassCutoff=4;
            obj.filt.DS4Hz.padding=true;
            obj.filt.DS4Hz=obj.filt.DS4Hz.designDownSample;
            
            obj.filt.FH=filterData(Fs);
            obj.filt.FH.highPassPassCutoff=100;
            obj.filt.FH.highPassStopCutoff=80;
            obj.filt.FH.lowPassPassCutoff=1800;
            obj.filt.FH.lowPassStopCutoff=2000;
            obj.filt.FH.attenuationInLowpass=20;
            obj.filt.FH.attenuationInHighpass=20;
            obj.filt.FH=obj.filt.FH.designBandPass;
            obj.filt.FH.padding=true;

            obj.filt.FHR=filterData(Fs);
            obj.filt.FHR.highPassPassCutoff=60;
            obj.filt.FHR.highPassStopCutoff=50;
            obj.filt.FHR.lowPassPassCutoff=900;
            obj.filt.FHR.lowPassStopCutoff=1000;
            obj.filt.FHR.attenuationInLowpass=20;
            obj.filt.FHR.attenuationInHighpass=40;
            obj.filt.FHR=obj.filt.FHR.designBandPass;
            obj.filt.FHR.padding=true;
            
            obj.filt.FL=filterData(Fs);
            obj.filt.FL.lowPassPassCutoff=4.5;
            obj.filt.FL.lowPassStopCutoff=6;
            obj.filt.FL.attenuationInLowpass=20;
            obj.filt.FL=obj.filt.FL.designLowPass;
            obj.filt.FL.padding=true;
            
            obj.filt.FH2=filterData(Fs);
            obj.filt.FH2.highPassCutoff=100;
            obj.filt.FH2.lowPassCutoff=2000;
            obj.filt.FH2.filterDesign='butter';
            obj.filt.FH2=obj.filt.FH2.designBandPass;
            obj.filt.FH2.padding=true;
        end
         
        %% checkFile - check the existance of a data file and a recording object
        function checkFileRecording(obj,fileName,message)
            if isempty(obj.currentDataObj)
                error('No data recording object selected!!!!');
            end
            if nargin>1
                if ~exist(fileName,'file')
                    if nargin==2
                        error('Relevant analysis file missing, please first run relevant function');
                    else
                        error(message);
                    end
                end
            end
        end

    end
    
end