function plotRaster(obj,params)

arguments (Input) 
    obj
    params.overwrite logical = false
    params.analysisTime = datetime('now')
    params.inputParams = false
    params.preBase = 500
    params.bin = 30
    params.exNeurons = 1
    params.AllSomaticNeurons = false
    params.AllResponsiveNeurons = false
    params.fixedWindow = true
    params.MergeNtrials =1
    params.GaussianLength = 3
    params.oneTrial = false
    params.imageDir = 'W:\Large_scale_mapping_NP\NormalAndRandImages'
    params.plotRawData = false
end


NeuronResp = obj.ResponseWindow;
Stats = obj.ShufflingAnalysis;
directimesSorted = NeuronResp.C(:,1)';

goodU = NeuronResp.goodU;
p = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
phy_IDg = p.phy_ID(string(p.label') == 'good');
pvals = Stats.pvalsResponse;

stimDur = NeuronResp.stimDur;
stimInter = NeuronResp.stimInter;

% trialDivision = numel(directimesSorted)/numel(unique(NeuronResp.C(:,2)))/numel(unique(NeuronResp.C(:,3)))/...
%     numel(unique(NeuronResp.C(:,4)));


preBase = params.preBase;

if params.AllSomaticNeurons
    eNeuron = 1:size(goodU,2);
    pvals = [eNeuron;pvals(eNeuron)];
elseif params.AllResponsiveNeurons
    eNeuron = find(pvals<0.05);
    pvals = [eNeuron;pvals(eNeuron)];% Select all good neurons if not specified
    if isempty(eNeuron)
        fprintf('No responsive neurons.\n')
        return
    end
else
    eNeuron = params.exNeurons;
    pvals = [eNeuron;pvals(eNeuron)];
end

bin=params.bin;
win=stimDur+preBase*2;
%preBase = round(stimInter/20)*10;

[Mr]=BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round(win/bin));

Mr = ConvBurstMatrix(Mr,fspecial('gaussian',[1 3],3),'same');

[nT,nN,nB] = size(Mr);
%indRG --> sorted infexes

trialsPerCath = length(directimesSorted);

ur =1;
for u = eNeuron

    if params.MergeNtrials >1
        j=1;
        mergeTrials = params.MergeNtrials;

        Mr2 = zeros(nT/mergeTrials,nB);

        for i = 1:mergeTrials:nT

            meanb = mean(squeeze(Mr(i:min(i+mergeTrials-1, end),u,:)),1);

            Mr2(j,:) = meanb;

            j = j+1;

        end
    else
        Mr2=Mr(:,u,:);
        mergeTrials =1;
    end

     if params.fixedWindow  %%Select highest window stim type
        j =1;
        meanMr = zeros(1,nT/trialsPerCath);
        for i = 1:trialsPerCath:nT
            meanMr(j) = mean(Mr2(i:i+trialsPerCath-1,:),'all');
            j = j+1;
        end

        [maxResp,maxRespIn]= max(meanMr);
        %Figure paper
        start = -50;
        window = stimDur+100;
    else
        [maxResp,maxRespIn]= max(NeuronResp.NeuronVals(u,:,1));
        start = NeuronResp.NeuronVals(u,maxRespIn,3)*NeuronResp.params.binRaster -20;
        window = 500;
     end

    [T,B] = size(Mr2);

  
    fig = figure;
    %Build raster
    M = Mr2.*(1000/bin);
    [nTrials,nTimes]=size(M);
    imagesc((1:nTimes),1:nTrials,squeeze(M));colormap(flipud(gray(64)));
    xline(preBase/bin, LineWidth=1.5, Color="#77AC30");
    xline((stimDur+preBase)/bin, LineWidth=1.5, Color="#0072BD");
    xticks([preBase/bin (round(stimDur/100)*100+preBase)/bin]);
    xticklabels(xticks*bin)

    yline([1:T]+0.5,LineWidth=1)

    set(fig, 'Color', 'w');
    % Set the color of the figure and axes to black
    colorbar;
    %caxis([0 1]);
    title(sprintf('NaturalImage-raster-U%d-PhyU-%dpval-%s',u,phy_IDg(u),num2str(pvals(2,ur),'%.4f')))
    ylabel(sprintf('%d trials',nTrials*mergeTrials))
    xlabel('Time (ms)')


    fig.Position =  [147   270   662   446];%[147    58   994   658];

    %%Plot raw data

    if params.plotRawData

        maxRespIn = maxRespIn-1;
        trials = maxRespIn*trialsPerCath+1:maxRespIn*trialsPerCath + trialsPerCath;

        chan = goodU(1,u);

        startTimes = directimesSorted(trials)+start;

        freq = "AP"; %or "LFP"

        typeData = "line"; %or heatmap

        fig2 = figure;

        spikes = squeeze(BuildBurstMatrix(goodU(:,u),round(p.t),round(startTimes),round((window))));

        if params.oneTrial
            [mx ind] = max(sum(spikes,2)); %select trial with most spikes
        else
            ind = 1:size(spikes,1);
        end

        [fig2, mx, mn] = PlotRawDataNP(obj,fig = fig2,chan = chan, startTimes = startTimes(ind),...
            window = window,spikeTimes = spikes(ind,:),multFactor =1.5, stdMult = 3);

        xline(-start/1000,'LineWidth',1.5,Color="#77AC30")
        xline((stimDur+abs(start))/1000,'LineWidth',1.5,Color="#0072BD")
        xticks([0,abs(start)/1000:abs(start*10)/1000:stimDur/1000+abs(start*2)/1000+1])
        xticklabels([start,0:abs(start*10):stimDur+abs(start*2)])
        xlabel('Miliseconds')
        yticks([])
        ylabel('uV')
        title(sprintf('U.%d-Unit-phy-%d-p-%d',u,phy_IDg(u),pvals(2,ur)));
        fig2.Position =  [147   270   662   446];

        if params.overwrite,obj.printFig(fig2,sprintf('%s-rect-GRid-rawData-raster-eNeuron-%d',obj.dataObj.recordingName,u)),close(fig2),end

    end

    if params.overwrite,obj.printFig(fig,sprintf('%s-rect-GRid-raster-eNeuron-%d',obj.dataObj.recordingName,u)),close(fig),end
    %prettify_plot
    
    ur = ur+1;

end %end eNeuron for loop

end %end plotRaster