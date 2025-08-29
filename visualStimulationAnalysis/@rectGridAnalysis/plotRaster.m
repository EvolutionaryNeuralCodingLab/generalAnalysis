function plotRaster(obj,params)

arguments (Input) 
    obj
    params.overwrite logical = false
    params.analysisTime = datetime('now')
    params.inputParams = false
    params.preBase = 200
    params.bin = 30
    params.exNeurons = 1
    params.AllSomaticNeurons = false
    params.AllResponsiveNeurons = false
    params.fixedWindow = false
    params.MergeNtrials =1
    params.GaussianLength = 3
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

positionsMatrix = [obj.VST.pos2X,obj.VST.pos2Y];%NewExp

seqMatrix = obj.VST.pos;

% trialDivision = numel(directimesSorted)/numel(unique(NeuronResp.C(:,2)))/numel(unique(NeuronResp.C(:,3)))/...
%     numel(unique(NeuronResp.C(:,4)));
nSize = numel(unique(NeuronResp.C(:,3)));
nLum = numel(unique(NeuronResp.C(:,4)));

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

trialsPerCath = length(seqMatrix)/(length(unique(seqMatrix)));

ur =1;
for u = eNeuron

    t = tiledlayout(sqrt(max(seqMatrix)), sqrt(max(seqMatrix)),'TileSpacing','compact');

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

    [T,B] = size(Mr2);
    j=1;
    for i = 1:trialsPerCath/mergeTrials:T
        %Build raster
        M = Mr2(i:min(i+trialsPerCath/mergeTrials-1, end),:,:).*(1000/bin);
        [nTrials,nTimes]=size(M);
        nexttile
        imagesc((1:nTimes),1:nTrials,squeeze(M));colormap(flipud(gray(64)));
        xline(preBase/bin, LineWidth=1.5, Color="#77AC30");
        xline((stimDur+preBase)/bin, LineWidth=1.5, Color="#0072BD");
        xticks([preBase/bin (round(stimDur/100)*100+preBase)/bin]);
        xticklabels(xticks*bin)
        if nSize >1
            yline([trialsPerCath/mergeTrials/nSize:trialsPerCath/mergeTrials/nSize:trialsPerCath/mergeTrials-1]+0.5,LineWidth=1)
        end

        if nLum >1
            yline([trialsPerCath/mergeTrials/nLum:trialsPerCath/mergeTrials/nLum:trialsPerCath/mergeTrials-1]+0.5,LineWidth=1)
        end

        set(gca,'YTickLabel',[]);

        if i < T - (trialsPerCath/mergeTrials)*max(positionsMatrix(:))-1
            set(gca,'XTickLabel',[]);

        end

        j = j+1;
    end
    fig = gcf;
    set(fig, 'Color', 'w');

    % Set the color of the figure and axes to black
    colorbar;
    title(t,sprintf('Rect-GRid-raster-U%d-pval-%s',u,num2str(pvals(2,ur),'%.4f')))
    ylabel(t,sprintf('%d trials',nTrials*mergeTrials))
    xlabel(t,'Time (ms)')
    fig.Position =  [147   270   662   446];%[147    58   994   658];

    if params.overwrite,obj.printFig(fig,sprintf('%s-rect-GRid-raster-eNeuron-%d',obj.dataObj.recordingName,u)),end
    %prettify_plot
    close

    ur = ur+1;



end %end eNeuron for loop

end %end plotRaster