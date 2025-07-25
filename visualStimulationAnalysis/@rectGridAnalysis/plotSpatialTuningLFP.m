function plotSpatialTuningLFP(obj,params)
%Generates a few plots with the position LFP tuning for the rectGrid stim

arguments (Input)
    obj
    params.normTrace = true %normalized the traces such the the mean of each trial is zero
    params.updatePlots = true %if true overwrites the existing plots. If false only generates the plots without saving
    params.overwrite logical = false %if true overwrites results
    params.analysisTime = datetime('now') %extract the time at which analysis was performed
    params.inputParams = false %if true - prints out the iput parameters so that it is clear what can be manipulated in the method
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
if params.updatePlots,obj.printFig(f1,'On_VsPosition'),end

f2=figure('Name','Off responses vs position');hL=tiledlayout(obj.VST.rectGridSize,obj.VST.rectGridSize,'TileSpacing','tight','Padding','tight','Visible','off');
for i=1:obj.VST.rectGridSize^2
    nexttile;
    imagesc(squeeze(mFLP_off(:,i,:)));
    set(gca,'XTick',[],'YTick',[],'CLim',[-stdOff*2 stdOff*2]);
end
hL.Visible='on';
if params.updatePlots,obj.printFig(f2,'Off_VsPosition'),end

f3=figure('Name','On averaged across positions');
imagesc((1:nBins)/LFP.samplingFreqLFP*1000,1:nElecs,squeeze(mean(mFLP_on,2)));
ylabel('Ch #');
xlabel('Time [ms]');
if params.updatePlots,obj.printFig(f3,'On_AvgOnPositions'),end

f4=figure('Name','Off averaged across positions');
imagesc((1:nBins)/LFP.samplingFreqLFP*1000,1:nElecs,squeeze(mean(mFLP_off,2)));
ylabel('Ch #');
xlabel('Time [ms]');
if params.updatePlots,obj.printFig(f4,'Off_AvgOnPositions'),end

f5=figure('Name','On averaged across electrodes');
h=axes;
activityTracePhysicalSpacePlot(h,1:numel(vStimGrid),squeeze(mean(mFLP_on,1)),vStimGrid,'traceColor',[0.8 0.2 0.2],'DrawElectrodeNumbers',0);
if params.updatePlots,obj.printFig(f5,'On_AvgOnElectrodes'),end

f6=figure('Name','Off averaged across electrodes');
h=axes;
vs.activityTracePhysicalSpacePlot(h,1:numel(vStimGrid),squeeze(mean(mFLP_off,1)),vStimGrid,'traceColor',[0.8 0.2 0.2],'DrawElectrodeNumbers',0);
if params.updatePlots,obj.printFig(f6,'Off_AvgOnElectrodes'),end

%[hPlot]=activityMultiTracePhysicalSpacePlot(h,(1:32)',squeeze(M_Int),vStimGrid,'DrawElectrodeNumbers',0,'traceColor',[0.8 0.2 0.2]);
%[hScaleBar]=addScaleBar(h,'scaleFac',1.4);
%activityTracePhysicalSpacePlot(h,IntanData.channelNumbers(:),squeeze(mean(M_Intan(:,pOrdered(1),:),2)),IntanData.chLayoutNumbers,'traceColor',[0.8 0.2 0.2],'DrawElectrodeNumbers',0);



end
