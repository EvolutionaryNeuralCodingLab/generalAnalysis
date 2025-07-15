classdef rectGridAnalysis < VStimAnalysis

    properties

    end

    properties (Constant)
        trialType = 'imageTrials'
    end

    methods (Hidden)
        %class constructor - name of class should be identical to the visual stimulation with the addition of Analysis
        function [obj] = rectGridAnalysis(dataObj)

            obj = obj.initialize(dataObj);
        end
    end

    methods

        function obj = plotSpatialTuningLFP(obj,params)
            arguments (Input)
                obj
                params.overwrite logical = false
                params.analysisTime = datetime('now')
                params.inputParams = false
                params.normTrace = true
                params.updatePlots = true
            end
            if params.inputParams,disp(params),return,end
            
            LFP=obj.getStimLFP;
            vStimGrid=nan(obj.VST.rectGridSize,obj.VST.rectGridSize);
            vStimGrid(sub2ind([obj.VST.rectGridSize',obj.VST.rectGridSize],obj.VST.pos2Y,obj.VST.pos2X))=1:(obj.VST.rectGridSize^2);
            vStimGrid=flipud(vStimGrid);

            [~,pOrdered]=sort(obj.VST.pos);
            [nElecs,nTrials,nBins]=size(LFP.LFP_on);
            mFLP_on=squeeze(mean(reshape(LFP.LFP_on(:,pOrdered,:),[nElecs,obj.VST.trialsPerCategory,nTrials/obj.VST.trialsPerCategory,nBins]),2));
            mFLP_off=squeeze(mean(reshape(LFP.LFP_off(:,pOrdered,:),[nElecs,obj.VST.trialsPerCategory,nTrials/obj.VST.trialsPerCategory,nBins]),2));
            
            if params.normTrace
                mFLP_on=bsxfun(@minus,mFLP_on,mean(mFLP_on,3));
                mFLP_off=bsxfun(@minus,mFLP_off,mean(mFLP_off,3));
            end

            stdOn=std(mFLP_on(:));
            stdOff=std(mFLP_off(:));

            f1=figure('Name','On responses vs position');hL=tiledlayout(obj.VST.rectGridSize,obj.VST.rectGridSize,'TileSpacing','tight','Padding','tight','Visible','off');
            for i=1:obj.VST.rectGridSize^2
                nexttile;
                imagesc(squeeze(mFLP_on(:,i,:)));
                set(gca,'XTick',[],'YTick',[],'CLim',[-stdOn*2 stdOn*2]);
            end
            hL.Visible='on';
            if params.updatePlots,obj.printFig(f1,'OnVsPosition'),end

            f2=figure('Name','Off responses vs position');hL=tiledlayout(obj.VST.rectGridSize,obj.VST.rectGridSize,'TileSpacing','tight','Padding','tight','Visible','off');
            for i=1:obj.VST.rectGridSize^2
                nexttile;
                imagesc(squeeze(mFLP_off(:,i,:)));
                set(gca,'XTick',[],'YTick',[],'CLim',[-stdOff*2 stdOff*2]);
            end
            hL.Visible='on';

            f3=figure('Name','On averaged across positions');
            imagesc((1:nBins)/LFP.samplingFreqLFP*1000,1:nElecs,squeeze(mean(mFLP_on,2)));
            ylabel('Ch #');
            xlabel('Time [ms]');

            f4=figure('Name','Off averaged across positions');
            imagesc((1:nBins)/LFP.samplingFreqLFP*1000,1:nElecs,squeeze(mean(mFLP_off,2)));
            ylabel('Ch #');
            xlabel('Time [ms]');

            f5=figure('Name','On averaged across electrodes');
            h=axes;
            activityTracePhysicalSpacePlot(h,1:numel(vStimGrid),squeeze(mean(mFLP_on,1)),vStimGrid,'traceColor',[0.8 0.2 0.2],'DrawElectrodeNumbers',0);

            f6=figure('Name','Off averaged across electrodes');
            h=axes;
            activityTracePhysicalSpacePlot(h,1:numel(vStimGrid),squeeze(mean(mFLP_off,1)),vStimGrid,'traceColor',[0.8 0.2 0.2],'DrawElectrodeNumbers',0);


            %[hPlot]=activityMultiTracePhysicalSpacePlot(h,(1:32)',squeeze(M_Int),vStimGrid,'DrawElectrodeNumbers',0,'traceColor',[0.8 0.2 0.2]);
            %[hScaleBar]=addScaleBar(h,'scaleFac',1.4);
            %activityTracePhysicalSpacePlot(h,IntanData.channelNumbers(:),squeeze(mean(M_Intan(:,pOrdered(1),:),2)),IntanData.chLayoutNumbers,'traceColor',[0.8 0.2 0.2],'DrawElectrodeNumbers',0);



        end
    end
end