classdef MEAPAnalysis < MEAAnalysis
    
    properties

    end
    
    properties (Constant)

    end
    
    methods
               
        %% MEAPAnalysis - class constructor
        function [obj]=MEAPAnalysis(xlsFile)
            if nargin==0
                xlsFile=[];
            end
            obj=obj@MEAAnalysis(xlsFile);
        end
        
        function [hand,data]=plotPatchTriggerSpikeOnMEA(obj,varargin)
            %default parameters
            parseObj = inputParser;
            addParameter(parseObj,'DrawElectrodeNumbers',0,@isnumeric);
            addParameter(parseObj,'scaling','std',@isstr);
            addParameter(parseObj,'scaleFac',[],@isnumeric);
            addParameter(parseObj,'DrawGrid',true,@isnumeric);
            addParameter(parseObj,'traceColor',[0 0 0.4]);
            addParameter(parseObj,'gridColor',[0.8 0.8 0.8],@isnumeric);
            addParameter(parseObj,'LineWidth',1,@isnumeric);
            addParameter(parseObj,'transparentScale',true,@isnumeric);
            addParameter(parseObj,'avgSubstruction',false,@isnumeric);
            addParameter(parseObj,'showScaleBar',true,@isnumeric);
            addParameter(parseObj,'substractMean',true,@isnumeric);
            addParameter(parseObj,'lockXYRatio',true,@isnumeric);
            addParameter(parseObj,'medFilterMs',0,@isnumeric);
            addParameter(parseObj,'averagingMethod','mean',@isnumeric);%'mean','median'
            addParameter(parseObj,'spikeOrder2Include',[],@isnumeric); %Plot only a subset of spikes in each sweep e.g. 3:5 includes the 3rd to 5th spikes if existing
            addParameter(parseObj,'ser',[],@isnumeric); %if empty plots the average of all series
            addParameter(parseObj,'channels2Remove',[],@isnumeric);
            addParameter(parseObj,'timeStartEnd',[],@isnumeric); %the start and end times to show relative to spike time
            
            addParameter(parseObj,'eventRejection',false,@isnumeric);
            addParameter(parseObj,'precentile4EventRejection',80,@isnumeric);
            addParameter(parseObj,'substractRandSIF',false,@isnumeric);
            addParameter(parseObj,'minPreviousISI',0,@isnumeric);
            
            addParameter(parseObj,'plotRankAverage',false,@isnumeric); % Add a plot with separate averages for different spike ranks
            addParameter(parseObj,'plotTimeAverage',false,@isnumeric); % Add a plot with separate averages for different spike times
            addParameter(parseObj,'showColorBar',false,@isnumeric); % show color bar in time and rank plots 
            addParameter(parseObj,'maxRank4RankAvgPlot',7,@isnumeric); %The maximal rank of the neuron to plot, if vecot takes intervals (eg. [1 3 6 10])1-2,3-5,6-9
            addParameter(parseObj,'intervals4TimeAvgPlot',250,@isnumeric); %The maximal rank of the neuron to plot (if vector takes the range of values
            addParameter(parseObj,'figureFileNameTimeAvg',[]);
            addParameter(parseObj,'figureFileNameRankAvg',[]);
            
            addParameter(parseObj,'saveFigures',1,@isnumeric);
            addParameter(parseObj,'figureFileName',[]);
            addParameter(parseObj,'h',0,@ishandle);
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            
            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams, disp(parseObj.Results), return, end
            P=parseObj.Results;
            
            [funName] = dbstack;funName=funName.name;funName=strsplit(funName,'.');funName=funName{end};
            %saveFileName=obj.files.(funName);
            
            obj.checkFileRecording(obj.files.getPatchTriggerSpikeOnMEA,'getPatchTriggerSpikeOnMEA file missing, please run getPatchTriggerSpikeOnMEA');
            
            if isempty(P.ser)
                PM=load(obj.files.getPatchTriggerSpikeOnMEA);
                pSeries=find(cellfun(@(x) ~isempty(x),PM.MAll));
                pars=PM.par;
                if P.substractRandSIF
                    PMR=load(obj.files.getPatchTriggerSpikeOnMEARandomized);
                end
                pSpks=1:numel(PM.spkRankAll);
            else
                PM = matfile(obj.files.getPatchTriggerSpikeOnMEA);
                [~,pSeries]=intersect(obj.recTable.Series(obj.currentPRec),P.ser);
                pars=PM.par;
                if P.substractRandSIF
                    PMR=load(obj.files.getPatchTriggerSpikeOnMEARandomized);
                end
                tmpNSpikes=PM.nSpikes;
                pSpks=(sum(cell2mat(tmpNSpikes(1:pSeries-1)))+1):sum(cell2mat(tmpNSpikes(1:pSeries)));
            end
            spkRankAll=PM.spkRankAll; %this is required since indexing cant be used with matfile
            spkTimeAll=PM.spkTimeAll; %this is required since indexing cant be used with matfile
            
            MAll=[];MAllR=[];spkOnMEAAll=[];
            for i=1:numel(pSeries)
                tmp=PM.MAll(1,pSeries(i));
                MAll=cat(2,MAll,tmp{1});
                tmp=cellfun(@(x) cell2mat(x),PM.spikeTimesOnMcd(1,pSeries(i)),'UniformOutput',0);
                spkOnMEAAll=[spkOnMEAAll tmp{1}];
                if P.substractRandSIF
                    tmp=PMR.MAll(1,pSeries(i));
                    MAllR=cat(2,MAllR,tmp{1});
                end
                %MAll=cat(2,MAll,matObj.MAll(pSeries(i)));
            end
            nTotSpks=size(MAll,2);
            
            if isempty(MAll)
                disp('No spikes identified intra-cellular data');
                return;
            end
            
            if ~isempty(P.spikeOrder2Include)
                pRel=false(1,nTotSpks);
                pRel(commonElements(spkRankAll(pSpks),P.spikeOrder2Include,0))=true;
            else
                pRel=true(1,nTotSpks);
            end
            pRel=pRel & [Inf diff(spkOnMEAAll)]>P.minPreviousISI;
            
            if P.substractMean
                MAll=bsxfun(@minus,MAll,mean(mean(MAll(:,pRel,:),2),3));
                %MAllR=bsxfun(@minus,MAllR,mean(mean(MAllR,2),3));
            end
            
            if P.eventRejection
                traceMAD=mean(mad(MAll(1:end-1,:,:),1,3),1);
                eventThresh=prctile(traceMAD,P.precentile4EventRejection);
                pRel=pRel & (traceMAD<eventThresh);
            end
            
            if P.substractRandSIF
                avgMEASIF=squeeze(mean(MAll(1:end-1,pRel,:),2))-squeeze(mean(MAllR(1:end-1,pRel,:),2));
            end
            
            if strcmp(P.averagingMethod,'mean')
                avgMEASIF=squeeze(mean(MAll(1:end-1,pRel,:),2));
            elseif strcmp(P.averagingMethod,'median')
                avgMEASIF=squeeze(median(MAll(1:end-1,pRel,:),2));
            end
            
            Fs_ms=obj.currentDataObj.samplingFrequency/1000;
            medFilter=round(P.medFilterMs*Fs_ms);
            if P.medFilterMs~=0
                avgMEASIF=medfilt1(avgMEASIF,medFilter,[],2);
            end
            
            pCh=1:numel(obj.currentDataObj.channelNumbers);
            if ~isempty(P.channels2Remove)
                [~,pCh]=setdiff(obj.currentDataObj.channelNumbers,P.channels2Remove);
            end
            
            pT=1:size(avgMEASIF,2);
            winPlot=pars.prePatchSpikeOnMEA+pars.postPatchSpikeOnMEA;
            if ~isempty(P.timeStartEnd)
                pStart=(pars.prePatchSpikeOnMEA+P.timeStartEnd(1))*Fs_ms;
                pEnd=(pars.postPatchSpikeOnMEA-P.timeStartEnd(2))*Fs_ms;
                pT=pT(pStart+1:end-pEnd);
                winPlot=diff(P.timeStartEnd);
            end
            
            %regular plot
            f1=[];f2=[];
            f=figure;h=axes;
            [hPlot, scale, En]=activityTracePhysicalSpacePlot(h,obj.currentDataObj.channelNumbers(pCh),avgMEASIF(pCh,pT),obj.currentDataObj.chLayoutNumbers,...
                'DrawElectrodeNumbers',P.DrawElectrodeNumbers,'DrawGrid',P.DrawGrid,'traceColor',P.traceColor,'gridColor',P.gridColor,'LineWidth',P.LineWidth,...
                'averageSubstructionSamples',1:(pars.prePatchSpikeOnMEA*Fs_ms),'scaling',P.scaling,'averageSubstruction',P.avgSubstruction,'lockXYRatio',P.lockXYRatio,'scaleFac',P.scaleFac);
            [yE,xE]=size(En);

            if P.showScaleBar
                [hScaleBar]=addScaleBar(h,'xLim_real',xE*[0 winPlot],'yLim_real',round([0 yE*scale(2)]),'scaleFac',2,'transparentScale',P.transparentScale);
            end
            
            if P.plotRankAverage
                if numel(P.maxRank4RankAvgPlot)==1
                    rankValues=mat2cell((1:P.maxRank4RankAvgPlot)',ones(1,P.maxRank4RankAvgPlot));
                else
                    for i=1:numel(P.maxRank4RankAvgPlot)-1
                        rankValues{i}=P.maxRank4RankAvgPlot(i):(P.maxRank4RankAvgPlot(i+1)-1);
                    end
                end
                nRanks=numel(rankValues);

                f1=figure;h=axes;
                cMap=jet(nRanks);
                tmpAvgSIF=squeeze(mean(MAll(1:end-1,:,:),2));
                [hPlot, scaleFac]=activityTracePhysicalSpacePlot(h,obj.currentDataObj.channelNumbers(pCh),tmpAvgSIF(pCh,pT),obj.currentDataObj.chLayoutNumbers,...
                        'DrawElectrodeNumbers',P.DrawElectrodeNumbers,'DrawGrid',P.DrawGrid,'traceColor','w','gridColor',P.gridColor,'LineWidth',P.LineWidth,...
                        'averageSubstructionSamples',1:(pars.prePatchSpikeOnMEA*Fs_ms),'scaling',P.scaling,'lockXYRatio',P.lockXYRatio,'scaleFac',P.scaleFac);
                delete(hPlot);
                
                for i=1:nRanks
                    tmpAvgSIF=squeeze(mean(MAll(1:end-1,commonElements(spkRankAll(pSpks),rankValues{i},0),:),2));
                    if medFilter~=0
                        tmpAvgSIF=medfilt1(tmpAvgSIF,medFilter,[],2);
                    end
                    [hPlot, scale, En]=activityTracePhysicalSpacePlot(h,obj.currentDataObj.channelNumbers(pCh),tmpAvgSIF(pCh,pT),obj.currentDataObj.chLayoutNumbers,...
                        'DrawElectrodeNumbers',P.DrawElectrodeNumbers,'DrawGrid',P.DrawGrid,'traceColor',cMap(i,:),'gridColor',P.gridColor,'LineWidth',P.LineWidth,...
                        'averageSubstructionSamples',1:(pars.prePatchSpikeOnMEA*Fs_ms),'scaling',P.scaling,'scaleFac',P.scaleFac,'lockXYRatio',P.lockXYRatio);
                end
                
                if P.showScaleBar
                    [hScaleBar]=addScaleBar(h,'xLim_real',xE*[0 winPlot],'yLim_real',round([0 yE*scale(2)]),'scaleFac',2,'transparentScale',P.transparentScale);
                end
                
                if P.showColorBar
                    cb=colorbar;
                    cb.Position=[ 0.833630952380952         0.707142857142857        0.0163690476190481         0.219047619047621];
                    colormap(cb,jet);
                    cb.XTick=(0:(nRanks-1))/(nRanks-1)*2;
                    cb.XTickLabel=cellfun(@(x) num2str(x),rankValues,'UniformOutput',0);
                    ylabel(cb,'Rank');
                end

            end
            
            if P.plotTimeAverage
                if numel(P.intervals4TimeAvgPlot)==1
                    tmp=round(spkTimeAll(pSpks)/P.intervals4TimeAvgPlot);
                    edegs = prctile(tmp,[5 95]) * P.intervals4TimeAvgPlot;
                    intervals=edegs(1):P.intervals4TimeAvgPlot:edegs(end);
                else
                    intervals=P.intervals4TimeAvgPlot;
                    %intervals=2000+20*(exp(0:5)-1);
                end
                nIntervals=numel(intervals)-1;
                
                f2=figure;h=axes;
                cMap=jet(nIntervals);
                tmpAvgSIF=squeeze(mean(MAll(1:end-1,:,:),2));
                [hPlot, scaleFac]=activityTracePhysicalSpacePlot(h,obj.currentDataObj.channelNumbers(pCh),tmpAvgSIF(pCh,pT),obj.currentDataObj.chLayoutNumbers,...
                        'DrawElectrodeNumbers',P.DrawElectrodeNumbers,'DrawGrid',P.DrawGrid,'traceColor','w','gridColor',P.gridColor,'LineWidth',P.LineWidth,...
                        'averageSubstructionSamples',1:(pars.prePatchSpikeOnMEA*Fs_ms),'scaling',P.scaling,'lockXYRatio',P.lockXYRatio,'scaleFac',P.scaleFac);
                    delete(hPlot);
                for i=1:nIntervals
                    tmpAvgSIF=squeeze(mean(MAll(1:end-1,spkTimeAll(pSpks)>=intervals(i) & spkTimeAll(pSpks)<intervals(i+1),:),2));
                    if medFilter~=0
                        tmpAvgSIF=medfilt1(tmpAvgSIF,medFilter,[],2);
                    end
                    [hPlot, scale, En]=activityTracePhysicalSpacePlot(h,obj.currentDataObj.channelNumbers(pCh),tmpAvgSIF(pCh,pT),obj.currentDataObj.chLayoutNumbers,...
                        'DrawElectrodeNumbers',P.DrawElectrodeNumbers,'DrawGrid',P.DrawGrid,'traceColor',cMap(i,:),'gridColor',P.gridColor,'LineWidth',P.LineWidth,...
                        'averageSubstructionSamples',1:(pars.prePatchSpikeOnMEA*Fs_ms),'scaling',P.scaling,'scaleFac',P.scaleFac,'lockXYRatio',P.lockXYRatio);
                end

                if P.showScaleBar
                    [hScaleBar]=addScaleBar(h,'xLim_real',xE*[0 winPlot],'yLim_real',round([0 yE*scale(2)]),'scaleFac',2,'transparentScale',P.transparentScale);
                end
                
                if P.showColorBar
                    cb=colorbar;
                    cb.Position=[ 0.833630952380952         0.707142857142857        0.0163690476190481         0.219047619047621];
                    colormap(cb,jet);
                    cb.XTick=(0:(numel(intervals)-1))/(numel(intervals)-1)*2;
                    cb.XTickLabel=num2str(round(intervals'));
                    ylabel(cb,'Time [ms]');
                end
            end

            if P.saveFigures
                figureFileName=obj.saveFigure(f,P.figureFileName);
                obj.saveFigure(f1,['rankAvg_' P.figureFileName]);
                obj.saveFigure(f2,['timeAvg_' P.figureFileName]);
            end
            
            data=[];hand=[];
            if nargout>1
                data.avgSIF=avgMEASIF;
                data.pCh=pCh;
                data.pT=pT;
            end
            if nargout>0
                hand=[f1 f2];
            end
            
        end
        %%
        function [data]=getSortingInSession(obj,varargin)
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
            
            %make parameter structure
            par=parseObj.Results;
            
            %check if analysis was already done done
            [funName] = dbstack;funName=funName.name;funName=strsplit(funName,'.');funName=funName{end};
            saveFileName=obj.files.(funName);
            if exist(saveFileName,'file') & ~overwrite
                if nargout==1
                    data=load(saveFileName);
                else
                    disp('Analysis results already exist for this method, use overwrite if needed');
                end
                return;
            end
            
            obj.populateGridSorterObj;
            S=obj.gridSorterObj.getSortedData;
            

            obj.checkFileRecording(obj.files.getSpikeSorting,'Patch resampling file missing, please run getPatchData');
            Pc=load(obj.files.getPatchData);
        end
        %% getPatchTriggerSpikeOnMEA
          function [data]=getPatchTriggerSpikeOnMEA(obj,varargin)
            
            %default parameters
            parseObj = inputParser;
            addParameter(parseObj,'minSpike2NoiseAmpDiff',10,@isnumeric);
            addParameter(parseObj,'spikePeakDetectionWinMs',5,@isnumeric);
            addParameter(parseObj,'minFilterCrossMs',0.2,@isnumeric);
            addParameter(parseObj,'minSpikePeakValue',0,@isnumeric);
            addParameter(parseObj,'prePatchSpikeOnMEA',20,@isnumeric);
            addParameter(parseObj,'postPatchSpikeOnMEA',100,@isnumeric);
            addParameter(parseObj,'maxMissingTriggers',5,@isnumeric);
            addParameter(parseObj,'medFiltDev',10,@isnumeric);
            addParameter(parseObj,'medFiltWin',30,@isnumeric);
            addParameter(parseObj,'triggerChannel',1,@isnumeric);
            addParameter(parseObj,'substractBaseline',1,@isnumeric);
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
            par=parseObj.Results;
            
            %check if analysis was already done done
            [funName] = dbstack;funName=funName.name;funName=strsplit(funName,'.');funName=funName{end};
            saveFileName=obj.files.(funName);
            if exist(saveFileName,'file') & ~overwrite
                if nargout==1
                    data=load(saveFileName);
                else
                    disp('Analysis results already exist for this method, use overwrite if needed');
                end
                return;
            end
            
            obj.checkFileRecording(obj.files.getPatchData,'Patch resampling file missing, please run getPatchData');
            Pc=load(obj.files.getPatchData);
            
            obj.checkFileRecording(obj.files.getDigitalTriggers,'Triggers file missing, please run getDigitalTriggers');
            T=load(obj.files.getDigitalTriggers);
            
            %get required parameters
            if ~ismember('SamplingCorrection', obj.recTable.Properties.VariableNames)
                disp('No resampling factor found, sampling is assumed to be correct!');
                resamplingFactor=ones(numel(obj.currentPRec),1);
            else
                resamplingFactor=obj.recTable.SamplingCorrection(obj.currentPRec);
            end
            
            % get parameters
            Fs_ms=obj.currentDataObj.samplingFrequency/1000;
            spikePeakDetectionWinSamples=spikePeakDetectionWinMs*Fs_ms;
            minFilterCrossSamples=minFilterCrossMs*Fs_ms;
            nSamples=(prePatchSpikeOnMEA+postPatchSpikeOnMEA)*Fs_ms;
            medFiltWinSamples=medFiltWin*Fs_ms;
            
            patchSeries=obj.recTable.Series(obj.currentPRec);
            clamp=obj.recTable.Clamp(obj.currentPRec);
            delay2Trig=obj.recTable.delay2Trig_ms(obj.currentPRec);
            protocol=obj.recTable.protocol(obj.currentPRec);
            
            %check if recordings are in current clamp mode
            pNonCurrentClamp=~strcmp(clamp,'CC');
            if any(pNonCurrentClamp)
                fprintf('The following series are not in current clamp and were not analyzed: %d\n',patchSeries(pNonCurrentClamp));
            end
            pValidPatch=find(~isnan(patchSeries) & ~pNonCurrentClamp )';
            
            %check for trigger validity
            nTotalSweeps=sum(Pc.nPatchSweeps);
            sweepNStart=[1;cumsum(Pc.nPatchSweeps(1:end-1))+1];
            if ~isempty(T.tTrig)
                nTotalTrigs=numel(T.tTrig{triggerChannel});
            else
                nTotalTrigs=0;
                disp('No triggers exists in this recording!!! Could not align patch to MEA');
                return;
            end
            
            if nTotalTrigs==0
                triggerCount=cellfun(@(x) numel(x),T.tTrig);
                newTriggerChannel=find(triggerCount>0,1,'first');
                if ~isempty(newTriggerChannel)
                    fprintf('Selected trigger channel %d is empty, trying with trigger channels %d\n',triggerChannel,newTriggerChannel);
                    triggerChannel=newTriggerChannel;
                    nTotalTrigs=numel(T.tTrig{triggerChannel});
                end
            end
            
            nMissingTrigs=nTotalSweeps-nTotalTrigs;
            if nTotalSweeps==nTotalTrigs
                badRecording=false;
                disp([num2str(nTotalSweeps) ' sweeps detected, ' num2str(nTotalTrigs) ' triggers detect - no triggers missing']);
                trig=T.tTrig{triggerChannel};
            else
                disp([num2str(nTotalSweeps) ' sweeps detected, ' num2str(nTotalTrigs) ' - Trigger mismatch!!!!!!!!!!']);
                if nMissingTrigs<=maxMissingTriggers
                    disp(['Assuming that the first ' num2str(nMissingTrigs) ' sweeps have no triggers - attempting to continue']);
                else
                    error('Trigger mismatch is larger than maxMissingTriggers! Please check recording!!!!!!');
                end
                badRecording=true;
                trig=[nan(1,nMissingTrigs) T.tTrig{triggerChannel}];
            end
            
            nSeries=numel(patchSeries);
            nSpikes=cell(1,nSeries);
            pSpikes=cell(1,nSeries);
            startTimesOnMcd=cell(1,nSeries);
            MAll=cell(1,nSeries);
            for k=pValidPatch %go over series in patch data
                disp(['Calculating patch triggered MEA data - ' obj.currentRecName ' , Series ' num2str(patchSeries(k))]);
                    
                    cumsumNSpikes=0;

                    for m=1:Pc.nPatchSweeps(k)
                        if badRecording && (m<=nMissingTrigs)
                            pSpikes{k}{m}=[];
                            startTimesOnMcd{k}{m}=[];
                            nSpikes{k}(m)=0;
                        else
                            %calculate the peak value of a spike
                            peakValue=max(Pc.Mp{k}(m,:));
                            
                            %identify potential threshold crossings
                            med = fastmedfilt1d(Pc.Mp{k}(m,:), medFiltWinSamples,-fliplr(Pc.Mp{k}(m,1:(medFiltWinSamples/2))),-fliplr(Pc.Mp{k}(m,end+1-(medFiltWinSamples/2):end)))';
                            %absMed=abs(Pc.Mp{k}(m,:)-med);medDev = fastmedfilt1d(absMed,medFiltWinSamples,-fliplr(absMed(1:(medFiltWinSamples/2))),-fliplr(absMed(end+1-(medFiltWinSamples/2):end)))'*1.4826;
                            patchSpikesLogical=Pc.Mp{k}(m,:)>(med+(peakValue-med)*0.5) & (peakValue-med)>minSpike2NoiseAmpDiff;
                            pUpCross=find((patchSpikesLogical(2:end)-patchSpikesLogical(1:end-1))==1);
                            pDnCross=find((patchSpikesLogical(2:end)-patchSpikesLogical(1:end-1))==-1);

                            %check that the duration of med filter crossing is not too short
                            if numel(pDnCross)==(numel(pUpCross))
                                pUpCross((pDnCross-pUpCross)<minFilterCrossSamples)=[]; %down crossing is not used
                            elseif numel(pDnCross)==(numel(pUpCross)-1)
                                pDnCross=[pDnCross numel(patchSpikesLogical)];
                                pUpCross((pDnCross-pUpCross)<minFilterCrossSamples)=[];
                            else
                                %does not remove samples from pUpCross
                            end
                            
                            %find peak of intra cellular spike
                            if numel(pUpCross)>0
                                spikeLocal=zeros(spikePeakDetectionWinSamples,numel(pUpCross));
                                pPeakMp=bsxfun(@plus,pUpCross,(1:spikePeakDetectionWinSamples)');
                                pPeakMp(pPeakMp>numel(Pc.tMp{k}))=numel(Pc.tMp{k}); %remove samples outside of recording
                                spikeLocal(:)=Pc.Mp{k}(m,pPeakMp);
                                spikeLocalSmooth=csaps(1:spikePeakDetectionWinSamples,spikeLocal',0.02,1:spikePeakDetectionWinSamples,ones(1,spikePeakDetectionWinSamples))';
                                [vMax,pMax]=max(spikeLocal); %find spike peak
                                pSpikes{k}{m}=pUpCross+pMax; %update spike position
                                pSpikes{k}{m}(pSpikes{k}{m}<minSpikePeakValue)=[]; %remove spike with low potential
                                nSpikes{k}(m)=numel(pSpikes{k}{m});
                                % A criterion with spike durations of at least 1ms should be added to false positive spikes
                                %{
                                figure('Position', [50 50 1200 800]);plot(Pc.Mp{k}(m,:));hold on;plot(patchSpikesLogical*3,'r');plot(med+(peakValue-med)*0.5,'g');plot(pSpikes{m},Mp(m,pSpikes{m}),'xk');xlim([0.4 0.8]*10^5)
                                %}
                                nCh=numel(obj.currentDataObj.channelNumbers);
                                M=zeros(nCh+1,numel(pSpikes{k}{m}),nSamples);
                                patchTmp=zeros(nSpikes{k}(m),nSamples);
                                p=bsxfun(@plus,pSpikes{k}{m}',-(prePatchSpikeOnMEA*Fs_ms):(postPatchSpikeOnMEA*Fs_ms-1));
                                p(p<=0)=1;p(p>Pc.nPatchSamples(k))=Pc.nPatchSamples(k);
                                patchTmp(:)=Pc.Mp{k}(m,p);
                                M(nCh+1,:,:)=patchTmp;
                                
                                startTimes=round(Fs_ms*(trig(sweepNStart(k)+m-1)-delay2Trig(k)*resamplingFactor(k)+pSpikes{k}{m}/Fs_ms-prePatchSpikeOnMEA))/Fs_ms;
                                %startTimes=round(Fs_ms*(trig(m)-delay2Trig(k)*resamplingFactor(k)+pSpikes{m}/Fs_ms-prePatchSpikeOnMEA))/Fs_ms;
                                M(1:nCh,:,:)=obj.currentDataObj.getData(obj.currentDataObj.channelNumbers,startTimes,prePatchSpikeOnMEA+postPatchSpikeOnMEA);
                                MAll{k}(:,cumsumNSpikes+1:cumsumNSpikes+nSpikes{k}(m),:)=M;
                                cumsumNSpikes=cumsumNSpikes+nSpikes{k}(m);
                            else
                                pSpikes{k}{m}=[];
                                nSpikes{k}(m)=0;
                                startTimes=[];
                            end
                            spikeTimesOnMcd{k}{m}=startTimes+prePatchSpikeOnMEA;
                            
                        end
                    end
                    
                    if ~isempty(MAll{k}) && substractBaseline
                        MAll{k}(1:nCh,:,:)=bsxfun(@minus,MAll{k}(1:nCh,:,:),mean(mean(MAll{k}(1:nCh,:,1:(prePatchSpikeOnMEA*Fs_ms)),2),3)); %remove average
                    end

            end %for loop over series
            
            spkRankAll=[];spkTimeAll=[];
            nSpk=nSpikes;% has to be done in 2 steps since PM is a matlab io class
            for i=pValidPatch
                pStart=[1 cumsum( nSpk{i}(1:end-1) )+1];
                pEnd=cumsum( nSpk{i} );
                spkRank=cell(1,numel(pStart));
                for j=1:numel(pStart)
                    spkRank{j}=1:(pEnd(j)-pStart(j)+1);
                end
                spkRankAll=[spkRankAll cell2mat(spkRank)];
                spkTimeAll=[spkTimeAll cell2mat(pSpikes{i})];
            end
            spkTimeAll=spkTimeAll/Fs_ms;
            
            STPatch=cellfun(@(x) x(end,:,:),MAll,'UniformOutput',1);
            save(saveFileName,'par','MAll','STPatch','nSpikes','pSpikes','spikeTimesOnMcd','badRecording','spkTimeAll','spkRankAll','-v7.3','-mat');
            
          end %runMEASTAs
          
          %% getPatchTriggerSpikeOnMEARandomized
          function [data]=getPatchTriggerSpikeOnMEARandomized(obj,varargin)
            
            %default parameters
            parseObj = inputParser;
            addParameter(parseObj,'randomization','uniform',@isstr);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            addParameter(parseObj,'overwrite',false,@isnumeric);

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
            
            %check if analysis was already done done
            [funName] = dbstack;funName=funName.name;funName=strsplit(funName,'.');funName=funName{end};
            saveFileName=obj.files.(funName);
            if exist(saveFileName,'file') & ~overwrite
                if nargout==1
                    disp('Loading previouse analyzed results');
                    data=load(saveFileName);
                else
                    disp('Analysis results already exist for this method, use overwrite if needed');
                end
                return;
            end
            
            obj.checkFileRecording(obj.files.getPatchTriggerSpikeOnMEA,'getPatchTriggerSpikeOnMEA file missing, please run getPatchTriggerSpikeOnMEA');
            load(obj.files.getPatchTriggerSpikeOnMEA,'par','nSpikes','pSpikes','spikeTimesOnMcd','badRecording','spkTimeAll','spkRankAll');
            
            obj.checkFileRecording(obj.files.getDigitalTriggers,'Triggers file missing, please run getDigitalTriggers');
            T=load(obj.files.getDigitalTriggers);
            
            delay2Trig=obj.recTable.delay2Trig_ms(obj.currentPRec);
            Fs_ms=obj.currentDataObj.samplingFrequency/1000;
            
            c=0;
            spkTimeAll=[];
            for i=1:numel(spikeTimesOnMcd)
                Ttmp=T.tTrig{1}(c+1:c+numel(spikeTimesOnMcd{i}));
                c=numel(Ttmp);
                spkTimeAllTmp=cell(1,numel(Ttmp));
                if strcmp(randomization,'uniform')
                    for j=1:numel(Ttmp)
                        intervals=diff([Ttmp(j) spikeTimesOnMcd{i}{j}]);
                        prm=randperm(numel(intervals));
                        
                        spkTimeAllTmp{j}=cumsum(intervals(prm));
                        permInt=Ttmp(j)+spkTimeAllTmp{j};
                        pSpikes{i}{j}=round((delay2Trig(i)+permInt-Ttmp(j))*Fs_ms);
                        spikeTimesOnMcd{i}{j}=permInt;
                    end
                    spkTimeAll=[spkTimeAll cell2mat(spkTimeAllTmp)];
                elseif strcmp(randomization,'shuffleTrials')
                    
                    prm=randperm(numel(Ttmp));
                    
                    pSpikes{i}=pSpikes{i}(prm);
                    nSpikes{i}=nSpikes{i}(prm);
                    
                    for j=1:numel(Ttmp)
                        spikeTimesOnMcd{i}{j}=Ttmp(j)+pSpikes{i}{j}/Fs_ms-delay2Trig(i);
                    end
                    
                    spkTimeAll=[];
                    spkRankAll=[];
                
                end
                MAll{i}=obj.currentDataObj.getData([],cell2mat(spikeTimesOnMcd{i})-par.prePatchSpikeOnMEA,par.postPatchSpikeOnMEA+par.prePatchSpikeOnMEA);
                MAll{i}(end+1,:,:)=0;
            end
            
            save(saveFileName,'par','MAll','nSpikes','pSpikes','spikeTimesOnMcd','badRecording','spkTimeAll','spkRankAll','-v7.3','-mat');
            
          end
          
          
          function plotICSpikeIdentification(obj)
              
              obj.checkFileRecording(obj.files.getPatchData,'Patch resampling file missing, please run getPatchData');
              Pc=load(obj.files.getPatchData,'Mp','tMp');
              
              obj.checkFileRecording(obj.files.getPatchTriggerSpikeOnMEA,'Patch resampling file missing, please run getPatchTriggerSpikeOnMEA');
              PM=load(obj.files.getPatchTriggerSpikeOnMEA,'pSpikes');
              
              patchSeries=obj.recTable.Series(obj.currentPRec);
              pValidPatch=find(~isnan(patchSeries))';
              
              for j=pValidPatch
                  vMax=max(Pc.Mp{j}(:));
                  vMin=min(Pc.Mp{j}(:));
                  f=figure;
                  imagesc(Pc.tMp{j},1:size(Pc.Mp{j},1),Pc.Mp{j},[max(vMax-50,vMin) vMax]);hold on;
                  for i=1:numel(PM.pSpikes{j})
                      if ~isempty(PM.pSpikes{j}{i})
                          plot(Pc.tMp{j}(PM.pSpikes{j}{i}),i,'.r')
                      end
                  end
                  xlabel('Time [ms]');
                  ylabel('Trial #');
              end
          end
          
          %% getPatchData
          function dataOut=getPatchData(obj,varargin)
              
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
            
            %make parameter structure
            par=parseObj.Results;
            
            %check if analysis was already done done
            [funName] = dbstack;funName=funName.name;funName=strsplit(funName,'.');funName=funName{end};
            saveFileName=obj.files.(funName);
            if exist(saveFileName,'file') & ~overwrite
                if nargout==1
                    dataOut=load(saveFileName);
                else
                    disp('Analysis results already exist for this method, use overwrite if needed');
                end
                return;
            end
            
            %get required parameters
            if ~ismember('SamplingCorrection', obj.recTable.Properties.VariableNames)
                disp('No resampling factor found, sampling is assumed to be correct!');
                resamplingFactor=ones(numel(obj.currentPRec),1);
            else
                resamplingFactor=obj.recTable.SamplingCorrection(obj.currentPRec);
            end

            % Extract patch data (including upsampling
            patchSeries=obj.recTable.Series(obj.currentPRec);
            pValidPatch=find(~isnan(patchSeries))';
            clamp=obj.recTable.Clamp(obj.currentPRec);
            delay2Trig=obj.recTable.delay2Trig_ms(obj.currentPRec);
            fullPatchDataFilename=cellfun(@(x,y) [x filesep num2str(y)],obj.recTable.folder(obj.currentPRec),obj.recTable.PatchFile(obj.currentPRec),'UniformOutput',0);
            hekaCh=obj.recTable.HekaCh(obj.currentPRec);
            if ~isempty(strmatch('hekaCurrCh',obj.recTable.Properties.VariableNames))
                hekaCurrCh=obj.recTable.hekaCurrCh(obj.currentPRec);
            else
                hekaCurrCh=nan(size(clamp));
            end
            %To add - if more than one series exists, increase the number of series by concatenating
            
            %check if to only analyze specific series, else analyzes all series with exluded series where patch data is not available
            nSeries=numel(patchSeries);

            Mp=cell(nSeries,1);tMp=cell(nSeries,1);units=cell(nSeries,1);
            MpIn=cell(nSeries,1);unitsIn=cell(nSeries,1);
            nPatchSamples=zeros(nSeries,1);nPatchSweeps=zeros(nSeries,1);
            for k=pValidPatch %go over series in patch data
                
                disp(['Resampling patch data ' obj.currentRecName ' , Series ' num2str(patchSeries(k)) ,' , Channel ' num2str(hekaCh(k))]);                
                
                %load and parse patch data
                if k==1 || ~strcmp(fullPatchDataFilename{k},fullPatchDataFilename{k-1})
                    load(fullPatchDataFilename{k},'data'); %it is assumned that the same Neuron will always belong to the same patch recording
                end
                [ch , t]= parseHeka(data,patchSeries(k),hekaCh(k));
                if isnan(hekaCurrCh(k))
                    [chI , tI]= parseHeka(data,patchSeries(k),hekaCh(k)+2);
                else
                    [chI , tI]= parseHeka(data,patchSeries(k),hekaCurrCh(k));
                end

                %get parameters
                units{k}=ch.units;
                unitsIn{k}=chI.units;
                nPatchSweeps(k)=size(ch.vm,2);
                Fs=obj.currentDataObj.samplingFrequency;
                
                %clear time vector for interpolation
                tMp{k}=(1000/Fs/resamplingFactor(k)):(1000/Fs/resamplingFactor(k)):max(t.ms);
                nPatchSamples(k)=numel(tMp{k});
                
                Mp{k}=zeros(nPatchSweeps(k),nPatchSamples(k)); %initialize data array
                for m=1:nPatchSweeps(k) %upsample data
                    Mp{k}(m,:) = interp1(t.ms,ch.vm(:,m),tMp{k},'spline');
                    MpIn{k}(m,:) = interp1(t.ms,chI.vm(:,m),tMp{k},'spline');
                end
                %assume that the input current is equal in all sweeps
                if all(clamp{k}=='CC') %check that recording is in current clamp
                    stimBinary=MpIn{k}(1,:)>max(MpIn{k}(1,:))/2;
                    stimOnTimes{k}=(ones(nPatchSweeps(k),1)*tMp{k}(stimBinary(2:end)==1 & stimBinary(1:end-1)==0))';
                    stimOffTimes{k}=(ones(nPatchSweeps(k),1)*tMp{k}(stimBinary(2:end)==0 & stimBinary(1:end-1)==1))';
                else
                    stimOnTimes{k}=[];
                    stimOffTimes{k}=[];
                end
                
            end %for loop over series
            save(saveFileName,'par','Mp','MpIn','tMp','units','unitsIn','stimOnTimes','stimOffTimes','nPatchSamples','nPatchSweeps','-v7.3','-mat');            
        end %getPatchData

    end
end