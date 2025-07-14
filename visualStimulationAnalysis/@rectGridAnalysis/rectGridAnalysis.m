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
            end


            LFP=vs.getStimLFP;
            vStimGrid=nan(vs.VST.rectGridSize,vs.VST.rectGridSize);
            vStimGrid(sub2ind([vs.VST.rectGridSize',vs.VST.rectGridSize],vs.VST.pos2Y,vs.VST.pos2X))=1:(vs.VST.rectGridSize^2);
            vStimGrid=flipud(vStimGrid);

            [~,pOrdered]=sort(vs.VST.pos);
            [nElecs,nTrials,nBins]=size(LFP.LFP_on);
    
            
            tmp=reshape(LFP.LFP_on(:,pOrdered,:),[nElecs,vs.VST.trialsPerCategory,nTrials/vs.VST.trialsPerCategory,nBins]);

            %tmp=reshape(permute(LFP.LFP_on(:,pOrdered,:),[3,1,2]),[nBins,nElecs,vs.VST.trialsPerCategory,nTrials/vs.VST.trialsPerCategory]);

            subplot(1,2,1);imagesc(squeeze(mean(tmp(:,:,2,:),2)));
            subplot(1,2,2);imagesc(squeeze(mean(LFP.LFP_on(:,pOrdered(16:30),:),2)));

            

            %permute([2 1 3]);
            %mFLP_on=squeeze(mean(LFP.LFP_on(:,pOrdered,:),2));

            mean(LFP.LFP_on)


            f=figure;h=axes;
            [hPlot]=activityMultiTracePhysicalSpacePlot(h,(1:32)',squeeze(M_Int),vStimGrid,'DrawElectrodeNumbers',0,'traceColor',[0.8 0.2 0.2]);
            [hScaleBar]=addScaleBar(h,'scaleFac',1.4);

                f=figure;h=axes;
                activityTracePhysicalSpacePlot(h,IntanData.channelNumbers(:),squeeze(mean(M_Intan(:,pOrdered(1),:),2)),IntanData.chLayoutNumbers,'traceColor',[0.8 0.2 0.2],'DrawElectrodeNumbers',0);
            %}
        end
    end
end