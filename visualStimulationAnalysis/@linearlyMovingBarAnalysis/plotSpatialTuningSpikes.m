function plotSpatialTuningSpikes(obj,params)
%Generates a few plots with the position spike tuning for the rectGrid stim
arguments (Input)
    obj
    params.bin = 10 %[ms] - bins size for the generated rasters
    params.win = [500,500] % duration [1,2] [ms] (for on and off) for LFP analysis
    params.getWinFromStimDuration = true %if this option is used, only the on response is calculated
    params.overwrite = true %if true overwrites the existing plots. If false only generates the plots without saving
    params.analysisTime = datetime('now') %extract the time at which analysis was performed
    params.inputParams = false %if true - prints out the iput parameters so that it is clear what can be manipulated in the method
end
if params.inputParams,disp(params),return,end

%load diode
stimTimes=obj.getSyncedDiodeTriggers;

s=obj.getSpikeTIcData;

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

%To add: for cases of low memory use a loop for calculating over groups of 10 trials and merge
if params.getWinFromStimDuration
    params.win=[obj.VST.stimDuration*1000, 500];
end

Mspk_on=BuildBurstMatrix(s.ic,round(s.t/params.bin),round(stimTimes.stimOnFlipTimes/params.bin),round(params.win(1)/params.bin));
Mspk_off=BuildBurstMatrix(s.ic,round(s.t/params.bin),round(stimTimes.stimOffFlipTimes/params.bin),round(params.win(2)/params.bin));

[nTrials,nNeu,nBins_on]=size(Mspk_on);
[~,~,nBins_off]=size(Mspk_off);

Mspk_on=reshape(Mspk_on(pOrdered,:,:),[obj.VST.trialsPerCategory,nTrials/obj.VST.trialsPerCategory,nNeu,nBins_on]);
Mspk_off=reshape(Mspk_off(pOrdered,:,:),[obj.VST.trialsPerCategory,nTrials/obj.VST.trialsPerCategory,nNeu,nBins_off]);

mSpk_on=squeeze(mean(Mspk_on,1));
mSpk_off=squeeze(mean(Mspk_off,1));

f1=figure('Name','On spike responses vs position','Position',[30 200 1800 600]);hL=tiledlayout(nAx1,nAx2,'TileSpacing','tight','Padding','compact','Visible','off');
for i=1:(nAx1*nAx2)
    nexttile;
    imagesc(squeeze(mSpk_on(i,:,:)));colormap(flipud(gray(16)))
    set(gca,'XTick',[],'YTick',[],'CLim',[0 2/params.bin]);
    if i<=nAx2
        title(num2str(Ax2Vals(i)),'FontSize',16)
    end
    if mod(i,nAx2)==1
        ylabel(Ax1Vals((i-1)/nAx2+1),'FontSize',16,'FontWeight','bold');
    end
end
hL.Visible='on';
if params.overwrite,obj.printFig(f1,'On_SpkVsPosition'),end

f2=figure('Name','Off spike responses vs position','Position',[30 200 1800 600]);hL=tiledlayout(nAx1,nAx2,'TileSpacing','tight','Padding','compact','Visible','off');
for i=1:(nAx1*nAx2)
    nexttile;
    imagesc(squeeze(mSpk_off(i,:,:)));colormap(flipud(gray(16)))
    set(gca,'XTick',[],'YTick',[],'CLim',[0 2/params.bin]);
    if i<=nAx2
        title(num2str(Ax2Vals(i)),'FontSize',16)
    end
    if mod(i,nAx2)==1
        ylabel(Ax1Vals((i-1)/nAx2+1),'FontSize',16,'FontWeight','bold');
    end
end
hL.Visible='on';
if params.overwrite,obj.printFig(f2,'Off_SpkVsPosition'),end

f3=figure('Name','On spike averaged across positions');
imagesc((1:nBins_on)*params.bin,1:nNeu,squeeze(mean(mSpk_on,2)));colormap(flipud(gray(16)));
clim([0 0.5/params.bin]);
ylabel('Ch #');
xlabel('Time [ms]');
if params.overwrite,obj.printFig(f3,'On_SpkAvgOnPositions'),end

f4=figure('Name','Off spike averaged across positions');
imagesc((1:nBins_off)*params.bin,1:nNeu,squeeze(mean(mSpk_off,2)));colormap(flipud(gray(16)));
clim([0 0.5/params.bin]);
ylabel('Ch #');
xlabel('Time [ms]');
if params.overwrite,obj.printFig(f4,'Off_SpkAvgOnPositions'),end

f5=figure('Name','On spike averaged across electrodes');
h=axes;
activityTracePhysicalSpacePlot(h,1:numel(vGrid),squeeze(mean(mSpk_on,2)),vGrid,'traceColor',[0.8 0.2 0.2],'DrawElectrodeNumbers',0);
if params.overwrite,obj.printFig(f5,'On_SpkAvgOnElectrodes'),end

f6=figure('Name','Off spike averaged across electrodes');
h=axes;
activityTracePhysicalSpacePlot(h,1:numel(vGrid),squeeze(mean(mSpk_off,2)),vGrid,'traceColor',[0.8 0.2 0.2],'DrawElectrodeNumbers',0);
if params.overwrite,obj.printFig(f6,'Off_SpkAvgOnElectrodes'),end

mM_on=squeeze(mean(Mspk_on,3));
f7=figure('Name','On spike all trials averaged across electrodes','Position',[30 200 1800 600]);
hL=tiledlayout(nAx1,nAx2,'TileSpacing','tight','Padding','compact','Visible','off');
for i=1:(nAx1*nAx2)
    nexttile;
    imagesc(squeeze(mM_on(:,i,:)));colormap(flipud(gray(16)))
    set(gca,'XTick',[],'YTick',[],'CLim',[0 1/params.bin]);
    if i<=nAx2
        title(num2str(Ax2Vals(i)),'FontSize',16)
    end
    if mod(i,nAx2)==1
        ylabel(Ax1Vals((i-1)/nAx2+1),'FontSize',16,'FontWeight','bold');
    end
end
hL.Visible='on';
if params.overwrite,obj.printFig(f7,'On_SpkTrialsVsPosition'),end

mM_off=squeeze(mean(Mspk_off,3));
f8=figure('Name','Off spike all trials averaged across electrodes','Position',[30 200 1800 600]);
hL=tiledlayout(nAx1,nAx2,'TileSpacing','tight','Padding','compact','Visible','off');
for i=1:(nAx1*nAx2)
    nexttile;
    imagesc(squeeze(mM_off(:,i,:)));colormap(flipud(gray(16)))
    set(gca,'XTick',[],'YTick',[],'CLim',[0 1/params.bin]);
       if i<=nAx2
        title(num2str(Ax2Vals(i)),'FontSize',16)
    end
    if mod(i,nAx2)==1
        ylabel(Ax1Vals((i-1)/nAx2+1),'FontSize',16,'FontWeight','bold');
    end
end
hL.Visible='on';
if params.overwrite,obj.printFig(f8,'Off_SpkTrialsVsPosition'),end
%{
f9=figure('Name','On spike all trials averaged across electrodes');
h=axes;
[hPlot]=activityMultiTracePhysicalSpacePlot(h,1:81,permute(squeeze(),[2,1,3]),vGrid,'scaling','noOverlap','DrawElectrodeNumbers',0,'traceColor',[0.8 0.2 0.2]);
%}

%[hPlot]=activityMultiTracePhysicalSpacePlot(h,(1:32)',squeeze(M_Int),vGrid,'DrawElectrodeNumbers',0,'traceColor',[0.8 0.2 0.2]);
%[hScaleBar]=addScaleBar(h,'scaleFac',1.4);
%activityTracePhysicalSpacePlot(h,IntanData.channelNumbers(:),squeeze(mean(M_Intan(:,pOrdered(1),:),2)),IntanData.chLayoutNumbers,'traceColor',[0.8 0.2 0.2],'DrawElectrodeNumbers',0);

end
