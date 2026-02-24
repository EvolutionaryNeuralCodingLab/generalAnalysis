function plotRaster(obj,params)

arguments (Input)
    obj
    params.overwrite logical = false
    params.analysisTime = datetime('now')
    params.inputParams = false
    params.preBase = 200
    params.bin = 40
    params.exNeurons = []
    params.AllSomaticNeurons = false
    params.AllResponsiveNeurons = true
    params.fixedWindow = true
    params.MergeNtrials =1
    params.GaussianLength = 50
    params.oneTrial = false
    params.selectedLum = []
    params.plotPatch logical = true
    params.PaperFig logical = false
    params.stim2show = 300
    
end


NeuronResp = obj.ResponseWindow;
Stats = obj.ShufflingAnalysis;
directimesSorted = NeuronResp.C(:,1)';

nSize = numel(unique(NeuronResp.C(:,3)));
nLum = numel(unique(NeuronResp.C(:,4)));

seqMatrix = obj.VST.pos;

if ~isempty(params.selectedLum)
    directimesSorted = directimesSorted(NeuronResp.C(:,4)==params.selectedLum);
    nLum = 1;
end

proportionTrials = 1/(numel(NeuronResp.C(:,1))/numel(directimesSorted));

goodU = NeuronResp.goodU;
p = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
phy_IDg = p.phy_ID(string(p.label') == 'good');
pvals = Stats.pvalsResponse;


stimDur = NeuronResp.stimDur;
stimInter = NeuronResp.stimInter;

positionsMatrix = [obj.VST.pos2X,obj.VST.pos2Y];%NewExp



% trialDivision = numel(directimesSorted)/numel(unique(NeuronResp.C(:,2)))/numel(unique(NeuronResp.C(:,3)))/...
%     numel(unique(NeuronResp.C(:,4)));


preBase = params.preBase;

if params.AllSomaticNeurons && isempty(params.exNeurons)
    eNeuron = 1:size(goodU,2);
    pvals = [eNeuron;pvals(eNeuron)];
elseif params.AllResponsiveNeurons && isempty(params.exNeurons)
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
win= preBase+params.stim2show; %stimDur+preBase*2;
%preBase = round(stimInter/20)*10;

[Mr]=BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round(win/bin));

%Mr = ConvBurstMatrix(Mr,fspecial('gaussian',[1 params.GaussianLength],3),'same');



ur =1;
for u = eNeuron

    [nT,nN,nB] = size(Mr);

    trialsPerCath = length(directimesSorted)/(length(unique(seqMatrix)));
    
    figure('Units', 'centimeters', 'Position', [10 10 6 4])

    t = tiledlayout(sqrt(max(seqMatrix)), sqrt(max(seqMatrix)),'TileSpacing','tight');

    if params.MergeNtrials >1
        j=1;
        mergeTrials = params.MergeNtrials;

        Mr2 = zeros(nT/mergeTrials,nB);

        for i = 1:mergeTrials:nT

            meanb = mean(squeeze(Mr(i:min(i+mergeTrials-1, end),u,:)),1);

            Mr2(j,:) = meanb;

            j = j+1;

        end
        trialsPerCath = trialsPerCath/mergeTrials;
        nT = nT/mergeTrials;
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
        window = 300;
    else
        [maxResp,maxRespIn]= max(NeuronResp.NeuronVals(u,:,1));
        start = NeuronResp.NeuronVals(u,maxRespIn,3)*NeuronResp.params.binRaster -20;
        window = 300;
    end

    [T,B] = size(Mr2);
    j=1;


    for i = 1:trialsPerCath:T
        %Build raster
        M = Mr2(i:min(i+trialsPerCath-1, end),:,:).*(1000/bin);
        [nTrials,nTimes]=size(M);

        ax = nexttile;

        imagesc((1:nTimes),1:nTrials,squeeze(M));colormap(flipud(gray(64)));

        %axis tight
        ax.Box = 'on';
        ax.XColor = [0.6 0.6 0.6];
        ax.YColor = [0.6 0.6 0.6];

        caxis([0 1]);   % choose meaningful limits
        set(gca, 'Color', 'w')

        %xline(preBase/bin, LineWidth=1.5, Color="#77AC30");
        %xline(preBase/bin, LineWidth=1.5,Color=[0.7 0.7 0.7])
        %xline((stimDur+preBase)/bin, LineWidth=1.5,Color=[0.7 0.7 0.7]);
        %xticks([preBase/bin (round(stimDur/100)*100+preBase)/bin]);
        %xticklabels(xticks*bin/1000)

        yl = ylim;
        
        patch([(preBase)/params.bin (preBase+params.stim2show)/params.bin (preBase+params.stim2show)/params.bin (preBase)/params.bin],...
                [yl(1) yl(1) yl(2) yl(2)],...
                'k','FaceAlpha',0.1,'EdgeColor','none')
        
        xticks([])
        ax = gca;
        ax.XAxis.FontSize = 8;
        ax.XAxis.FontName = 'helvetica';

        if nSize >1
            yline([trialsPerCath/nSize:trialsPerCath/nSize:trialsPerCath-1]+0.5,LineWidth=1)
        end

        if nLum >1
            yline([trialsPerCath/nLum:trialsPerCath/nLum:trialsPerCath-1]+0.5,LineWidth=1)
        end

        set(gca,'YTickLabel',[]);

        if i < T - (trialsPerCath)*max(positionsMatrix(:))-1
            set(gca,'XTickLabel',[]);

        end

        if j ==  maxRespIn && params.plotPatch && ~params.oneTrial

            ylims = ylim;

            patch([(preBase+start)/params.bin (preBase+start+window)/params.bin (preBase+start+window)/params.bin (preBase+start)/params.bin]...
                ,[ylims(1) ylims(1) ylims(2) ylims(2)],'r','FaceAlpha',0.3,'EdgeColor','none')


        elseif j ==  maxRespIn && params.plotPatch && params.MergeNtrials == 1

            [mx ind] = max(sum(squeeze(M(:,:,round((preBase+start)/params.bin):round((preBase+start+window)/params.bin))),2));

            patch([(preBase+start)/params.bin (preBase+start+window)/params.bin (preBase+start+window)/params.bin (preBase+start)/params.bin],...
                [ind-0.5 ind-0.5 ind+0.5 ind+0.5],...
                'r','FaceAlpha',0.3,'EdgeColor','none')
        end

        j = j+1;
    end

    fig = gcf;
    set(fig, 'Color', 'w');
    % % Set the color of the figure and axes to black
    % cb = colorbar;
    % cb.FontName = 'helvetica';
    % cb.FontSize = 8;
    % title(cb, '[spk/s]', 'FontSize', 10, 'FontName','helvetica');
    
    
    %%%% scale bar (outside, below last tile – same geometry)

    hold(ax, 'on')

    % ---- data-space limits ----
    xl = xlim(ax);
    yl = ylim(ax);

    % ---- parameters (UNCHANGED from your logic) ----
    xLen = params.stim2show / bin;          % ms
    yLen = 0.12 * range(yl);
    marginY = 0.1;

    % ---- anchor in DATA coordinates (same as before) ----
    x0 = preBase / bin;
    y0 = max(yl) - marginY * range(yl); %#ok<NASGU> (kept for clarity)

    % ---- convert DATA-x to AXES-normalized coordinates ----
    xNorm0   = (x0   - xl(1)) / range(xl);
    xNormLen =  xLen / range(xl);

    % ---- get axes position in FIGURE coordinates ----
    axPos = ax.Position;    % [x y width height]

    % ---- place bar BELOW the tile ----
    gap = 0.02;             % gap below tile (figure units)
    yFig = axPos(2) - gap;

    % ---- horizontal bar (above text) ----
    annotation('line', ...
        axPos(1) + [xNorm0 xNorm0 + xNormLen] * axPos(3), ...
        [yFig yFig], ...
        'Color','k', 'LineWidth', 2)

    % ---- vertical bar (short, downward, right end) ----
    yLenFig = 0.12 * axPos(4);
    dx = 0.02;  % fraction of tile width to shift left
    annotation('line', ...
        axPos(1) + ((xNorm0 + xNormLen - dx) * axPos(3)) * [1 1], ...
        [yFig yFig - yLenFig-yLenFig*0.3], ...
        'Color','k', 'LineWidth', 2)

    % ---- label (below the horizontal bar) ----
    textOffset = 0.050;  % fraction of figure height to move text below the horizontal bar

    annotation('textbox', ...
        [axPos(1) + xNorm0 * axPos(3), ...
        yFig - textOffset, ...        % moved down from the horizontal line
        xNormLen * axPos(3), ...
        0.05], ...
        'String', sprintf('%d ms', (round(params.stim2show/10)*10)), ...
        'EdgeColor','none', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize', 10, 'FontName','helvetica', ...
        'Interpreter','none','FitBoxToText','on')


    hold(ax, 'off')

    title(t,sprintf('Rect-GRid-raster-U%d-PhyU-%dpval-%s',u,phy_IDg(u),num2str(pvals(2,ur),'%.4f')))
    ylabel(t,sprintf('%d Trials/Location',nTrials*mergeTrials),'FontSize',10,'FontName','helvetica')
    xlabel(t,'Time','FontSize',10,'FontName','helvetica')
    
    set(fig, 'Units', 'centimeters');
    set(fig, 'Position', [20 20 6 6]);
    
    
   % fig.Position =  [1120         651         466         307];%[147    58   994   658];

    % pos1 = cb.Position(1);
    % cb.Position(1) = pos1 + 0.03;


    if params.PaperFig
        obj.printFig(fig,sprintf('%s-rect-GRid-raster-eNeuron-%d',obj.dataObj.recordingName,u),PaperFig = params.PaperFig)
    elseif params.overwrite
        obj.printFig(fig,sprintf('%s-rect-GRid-raster-eNeuron-%d',obj.dataObj.recordingName,u))
    end


    %%Plot raw data

    maxRespIn = maxRespIn-1;

    trialsPerCath = length(directimesSorted)/(length(unique(seqMatrix)));
    trials = maxRespIn*trialsPerCath+1:maxRespIn*trialsPerCath + trialsPerCath;

    chan = goodU(1,u);

    startTimes = directimesSorted(trials)+start;

    freq = "AP"; %or "LFP"

    typeData = "line"; %or heatmap

    spikes = squeeze(BuildBurstMatrix(goodU(:,u),round(p.t),round(startTimes),round((window))));

    if params.oneTrial
        [mx ind] = max(sum(spikes,2)); %select trial with most spikes
    else
        ind = 1:size(spikes,1);
    end
    fig2 = figure;

    [fig2, mx, mn] = PlotRawDataNP(obj,fig = fig2,c = chan, startTimes = startTimes(ind),...
        window = window,spikeTimes = spikes(ind,:));
    %
    % xline(-start/1000,'r','LineWidth',1.5)
    % xline((stimDur+abs(start))/1000,'r','LineWidth',1.5)
    xticks([0,1/10:1/10:window/1000])
    xticklabels([0:abs(100):obj.VST.stimDuration*1000+abs(start*2)])
    xlabel('Time [ms]','FontSize',10,'FontName','helvetica')
    ax = gca;
    ax.XAxis.FontSize = 8;
    ax.XAxis.FontName = 'helvetica';
    ax.YAxis.FontSize = 8;
    ax.YAxis.FontName = 'helvetica';

    ylabel('[\muV]','FontSize',10,'FontName','helvetica')
    title(sprintf('U.%d-Unit-phy-%d-p-%d',u,phy_IDg(u),pvals(2,ur)));
    
    if params.oneTrial
        tr = 1;
        set(fig2, 'Units', 'centimeters');
        set(fig2, 'Position', [20 20 10 3]);
    else
        tr = numel(ind);
    end

    if params.PaperFig
        obj.printFig(fig2,sprintf('%s-rect-GRid-rawData-%d-Trials-raster-eNeuron-%d',obj.dataObj.recordingName,tr,u),PaperFig = params.PaperFig)
    elseif params.overwrite
        obj.printFig(fig2,sprintf('%s-rect-GRid-rawData-%d-Trials-raster-eNeuron-%d',obj.dataObj.recordingName,u))
    end

    %prettify_plot

    if u ~= eNeuron(end)
        close all
    end

    ur = ur+1;


end %end eNeuron for loop

end %end plotRaster