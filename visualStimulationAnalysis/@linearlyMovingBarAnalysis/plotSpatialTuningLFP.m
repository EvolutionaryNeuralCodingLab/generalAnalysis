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

%add code to create multiple plots in case that more than 1 ball size or velocity was used.
nSpeeds=numel(obj.VST.speed);
nOffsets=numel(obj.VST.parallelsOffset);
nDirections=obj.VST.numberOfDirections;

Ax1Trials=obj.VST.directions;
Ax2Trials=obj.VST.offsets;
Ax1Vals=unique(obj.VST.directions);
Ax2Vals=unique(obj.VST.offsets);
nAx1=numel(Ax1Vals);
nAx2=numel(Ax2Vals);

vGrid=nan(nAx2,nAx1);
vGrid(:)=1:(nAx1*nAx2);
vGrid=flipud(vGrid');
[Xgrid,Ygrid]=meshgrid(1:nAx2,1:nAx1);

[~,pOrdered2]=sort(Ax2Trials);
[originalOrder,pOrdered1]=sort(Ax1Trials(pOrdered2));
pOrdered=pOrdered2(pOrdered1);

%figure;plotyy(1:numel(pOrder),Ax2Trials(pOrder),1:numel(pOrder),Ax1Trials(pOrder));

[nElecs,nTrials,nBins_on]=size(LFP.LFP_on);
[~,~,nBins_off]=size(LFP.LFP_off);

mFLP_on=squeeze(mean(reshape(LFP.LFP_on(:,pOrdered,:),[nElecs,obj.VST.trialsPerCategory,nTrials/obj.VST.trialsPerCategory,nBins_on]),2));
mFLP_off=squeeze(mean(reshape(LFP.LFP_off(:,pOrdered,:),[nElecs,obj.VST.trialsPerCategory,nTrials/obj.VST.trialsPerCategory,nBins_off]),2));

if params.normTrace
    mFLP_on=bsxfun(@minus,mFLP_on,mean(mFLP_on,3));
    mFLP_off=bsxfun(@minus,mFLP_off,mean(mFLP_off,3));
end

stdOn=std(mFLP_on(:));
stdOff=std(mFLP_off(:));

f1=figure('Name','On responses vs position','Position',[30 200 1800 600]);hL=tiledlayout(nAx1,nAx2,'TileSpacing','tight','Padding','compact','Visible','off');
for i=1:(nAx1*nAx2)
    nexttile;
    imagesc(squeeze(mFLP_on(:,i,:)));
    set(gca,'XTick',[],'YTick',[],'CLim',[-stdOn*2 stdOn*2]);
    if i<=nAx2
        title(num2str(Ax2Vals(i)),'FontSize',16)
    end
    if mod(i,nAx2)==1
        ylabel(Ax1Vals((i-1)/nAx2+1),'FontSize',16,'FontWeight','bold');
    end
end
hL.Visible='on';
if params.updatePlots,obj.printFig(f1,'On_VsPosition'),end

f2=figure('Name','Off responses vs position','Position',[30 200 1800 600]);hL=tiledlayout(nAx1,nAx2,'TileSpacing','tight','Padding','compact','Visible','off');
for i=1:(nAx1*nAx2)
    nexttile;
    imagesc(squeeze(mFLP_off(:,i,:)));
    set(gca,'XTick',[],'YTick',[],'CLim',[-stdOff*2 stdOff*2]);
    if i<=nAx2
        title(num2str(Ax2Vals(i)),'FontSize',16)
    end
    if mod(i,nAx2)==1
        ylabel(Ax1Vals((i-1)/nAx2+1),'FontSize',16,'FontWeight','bold');
    end
end
hL.Visible='on';
if params.updatePlots,obj.printFig(f2,'Off_VsPosition'),end

f3=figure('Name','On averaged across positions','Position',[30 200 1000 300]);
imagesc((1:nBins_on)/LFP.samplingFreqLFP*1000,1:nElecs,squeeze(mean(mFLP_on,2)));
ylabel('Ch #');
xlabel('Time [ms]');
if params.updatePlots,obj.printFig(f3,'On_AvgOnPositions'),end

f4=figure('Name','Off averaged across positions','Position',[30 200 1000 300]);
imagesc((1:nBins_off)/LFP.samplingFreqLFP*1000,1:nElecs,squeeze(mean(mFLP_off,2)));
ylabel('Ch #');
xlabel('Time [ms]');
if params.updatePlots,obj.printFig(f4,'Off_AvgOnPositions'),end

f5=figure('Name','On averaged across electrodes','Position',[30 200 1800 600]);
h=axes;
activityTracePhysicalSpacePlot(h,1:numel(vGrid),squeeze(mean(mFLP_on,1)),vGrid,'traceColor',[0.8 0.2 0.2],'DrawElectrodeNumbers',0);
if params.updatePlots,obj.printFig(f5,'On_AvgOnElectrodes'),end

f6=figure('Name','Off averaged across electrodes','Position',[30 200 1800 600]);
h=axes;
activityTracePhysicalSpacePlot(h,1:numel(vGrid),squeeze(mean(mFLP_off,1)),vGrid,'traceColor',[0.8 0.2 0.2],'DrawElectrodeNumbers',0);
if params.updatePlots,obj.printFig(f6,'Off_AvgOnElectrodes'),end

%[hPlot]=activityMultiTracePhysicalSpacePlot(h,(1:32)',squeeze(M_Int),vGrid,'DrawElectrodeNumbers',0,'traceColor',[0.8 0.2 0.2]);
%[hScaleBar]=addScaleBar(h,'scaleFac',1.4);
%activityTracePhysicalSpacePlot(h,IntanData.channelNumbers(:),squeeze(mean(M_Intan(:,pOrdered(1),:),2)),IntanData.chLayoutNumbers,'traceColor',[0.8 0.2 0.2],'DrawElectrodeNumbers',0);

end
