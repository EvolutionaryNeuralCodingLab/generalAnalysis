function plotRaster(obj, params)
% plotRaster  Natural-image raster, PSTH, and raw trace figure.
%
% Layout (manually positioned axes — no subplot):
%
%    +-----------------------------+  top = 1 - topMargin
%    |                             |
%    |     RASTER     |  THUMB     |
%    |                             |
%    +-----------------------------+  psthTop + midGap
%    |        midGap (empty)       |
%    +-----------------------------+  psthTop
%    |           PSTH              |
%    +-----------------------------+  psthBottom = bottomMargin
%    |   bottomMargin (xlabel)     |
%    +-----------------------------+  y = 0
%
% ResponseWindow fields used:
%   C columns = [stimOnTime, imageIndex, state]
%   stimDur, stimInter, NeuronVals

% --------------------------------------------------------------------------
% BUG FIXES vs original:
%  1. Mr2 = Mr(:,u,:) left 3-D when MergeNtrials==1 -> squeeze().
%  2. Row indices into Mr2 used raw trialsPerCath when merged
%     -> rowsPerCath = trialsPerCath/mergeTrials.
%  3. maxRespIn 0-based shift applied before trial arithmetic -> clarified.
%  4. 30x10 subplot grid on a 5 cm-tall figure collapsed raster and PSTH
%     onto each other -> replaced with explicit ax.Position overrides
%     in normalised figure units.
% --------------------------------------------------------------------------

arguments (Input)
    obj
    params.overwrite            logical = false
    params.analysisTime                 = datetime('now')
    params.inputParams          logical = false
    params.preBase                      = 500             % Pre/post-stimulus baseline (ms)
    params.bin                          = 20              % Raster bin size (ms/bin)
    params.exNeurons                    = 1               % Neuron index into good-unit list
    params.AllSomaticNeurons    logical = false
    params.AllResponsiveNeurons logical = false
    params.fixedWindow          logical = true            % true: fixed window; false: NeuronVals
    params.MergeNtrials                 = 1               % Trials averaged per raster row
    params.GaussianLength               = 3               % Gaussian kernel length (bins)
    params.oneTrial             logical = false           % Raw: only trial with most spikes
    params.imageDir                     = 'W:\Large_scale_mapping_NP\NormalAndRandImages'
    params.windowDur                    = 500             % Sliding-window width (ms)
    params.psthBinWidth                 = 100              % PSTH bin width (ms)
    params.MaxVal_1                     = true
    params.plotRawData                  = false
    params.selectCats                   = []
    params.PaperFig                     = false
    params.bottomMargin                 = 0.18            % xlabel space (norm. units)
    params.midGap                       = 0.04            % Gap raster<->PSTH (norm. units)
    params.psthHeightFrac               = 0.18            % PSTH height (norm. units)
end

if params.inputParams, disp(params); return; end

% ==========================================================================
% 1. LOAD PRE-COMPUTED RESULTS
% ==========================================================================

NeuronResp = obj.ResponseWindow;
Stats      = obj.StatisticsPerNeuron;

C = NeuronResp.C;
nImages = numel(unique(C(:,2)));
nState  = numel(unique(C(:,3)));
directimesSorted = C(:,1)';

trialsPerCath = size(C,1) / nImages;   % Unmerged trials per image category

% ==========================================================================
% 2. IMAGE METADATA AND LOADING
% ==========================================================================

imagesNames = cell(1, numel(obj.VST.imgNames));
imagesNames(1:numel(imagesNames)/2)     = obj.VST.imgNames(1:2:end);
imagesNames(numel(imagesNames)/2+1:end) = obj.VST.imgNames(2:2:end);

cd(params.imageDir);
imgs        = cellfun(@imread, imagesNames, 'UniformOutput', false);
combinedImg = cat(1, imgs{:});

if ~isempty(params.selectCats)
    indexes = [];
    for i = 1:numel(params.selectCats)
        indexes = [indexes, ...
            params.selectCats(i)*trialsPerCath - (trialsPerCath-1) : params.selectCats(i)*trialsPerCath];
    end
    C                = C(indexes, :);
    nImages          = numel(unique(C(:,2)));
    nState           = numel(unique(C(:,3)));
    imagesNames      = imagesNames(params.selectCats);
    imgs             = cellfun(@imread, imagesNames, 'UniformOutput', false);
    combinedImg      = cat(1, imgs{:});
    directimesSorted = C(:,1)';
end

% ==========================================================================
% 3. SPIKE SORTING METADATA
% ==========================================================================

p       = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
phy_IDg = p.phy_ID(string(p.label') == 'good');
pvals   = Stats.pvalsResponse;
label   = string(p.label');
goodU   = p.ic(:, label == 'good');

% ==========================================================================
% 4. STIMULUS TIMING
% ==========================================================================

stimDur = NeuronResp.stimDur;
preBase = params.preBase;
bin     = params.bin;
win     = stimDur + preBase*2;

% ==========================================================================
% 5. NEURON SELECTION
% ==========================================================================

if params.AllSomaticNeurons
    eNeuron = 1:size(goodU, 2);
    pvals   = [eNeuron; pvals(eNeuron)];
elseif params.AllResponsiveNeurons
    eNeuron = find(pvals < 0.05);
    pvals   = [eNeuron; pvals(eNeuron)];
    if isempty(eNeuron)
        fprintf('No responsive neurons.\n'); return
    end
else
    eNeuron = params.exNeurons;
    pvals   = [eNeuron; pvals(eNeuron)];
end

% ==========================================================================
% 6. BUILD RASTER MATRIX
% ==========================================================================

Mr = BuildBurstMatrix(goodU, round(p.t/bin), ...
    round((directimesSorted - preBase)/bin), round(win/bin));
Mr = ConvBurstMatrix(Mr, fspecial('gaussian', [1 params.GaussianLength], 3), 'same');

[nT, ~, nB] = size(Mr);

% ==========================================================================
% 7. PER-NEURON FIGURE LOOP
% ==========================================================================

ur = 1;

for u = eNeuron

    % ------------------------------------------------------------------
    % 7a. 2-D raster Mr2 (BUG FIX 1: squeeze)
    % ------------------------------------------------------------------

    mergeTrials = params.MergeNtrials;

    if mergeTrials > 1
        nRows = floor(nT / mergeTrials);
        Mr2   = zeros(nRows, nB);
        for i = 1:nRows
            src      = (i-1)*mergeTrials + 1 : min(i*mergeTrials, nT);
            Mr2(i,:) = mean(squeeze(Mr(src, u, :)), 1);
        end
    else
        Mr2   = squeeze(Mr(:, u, :));
        nRows = nT;
    end

    rowsPerCath = trialsPerCath / mergeTrials;   % BUG FIX 2

    % ------------------------------------------------------------------
    % 7b. Best image category
    % ------------------------------------------------------------------

    meanMr = zeros(1, nImages);
    for i = 1:nImages
        r1 = (i-1)*rowsPerCath + 1;
        r2 =  i   *rowsPerCath;
        meanMr(i) = mean(Mr2(r1:r2, :), 'all');
    end

    [~, bestCat]    = max(meanMr);
    bestCatRowStart = (bestCat-1)*rowsPerCath + 1;
    bestCatRowEnd   =  bestCat   *rowsPerCath;
    trialsAbsolute  = (bestCat-1)*trialsPerCath + 1 : bestCat*trialsPerCath;

    % ------------------------------------------------------------------
    % 7c. Sliding-window search
    % ------------------------------------------------------------------

    window   = params.windowDur;
    nWinBins = round(window / bin);
    X        = min(Mr2(bestCatRowStart:bestCatRowEnd, :), 1);
    nWinPos  = nB - nWinBins + 1;

    window_means = zeros(rowsPerCath, nWinPos);
    for col = 1:nWinPos
        window_means(:, col) = mean(X(:, col:col+nWinBins-1), 2);
    end

    [~, linear_idx]      = max(window_means(:));
    [best_row, best_col] = ind2sub(size(window_means), linear_idx);

    bestDisplayRow = bestCatRowStart + best_row - 1;
    bestTrialIdx   = trialsAbsolute((best_row-1)*mergeTrials + 1);

    % ==========================================================================
    % 8. FIGURE — explicit axes positions in normalised figure units
    %
    % BUG FIX 4: the 30x10 subplot grid put the PSTH at tile rows 21-30 and
    % the raster at rows 1-18. On a 5 cm-tall figure each tile is <2 mm high
    % and MATLAB's automatic inset margins around each subplot are larger
    % than the tile itself — so the raster and PSTH visually overlap.
    %
    % Fix: compute three panel rectangles directly in normalised figure units
    % and create axes with axes('Position', [left bottom width height]).
    % This guarantees the PSTH sits below the raster with a clean gap
    % (midGap) and leaves an explicit bottom strip (bottomMargin) for
    % xlabel and xticklabels to render without being clipped.
    % ==========================================================================

    fig = figure;
    set(fig, 'Units', 'centimeters');
    set(fig, 'Position', [20 20 9 4]);

    % Horizontal layout
    leftMargin  = 0.12;
    rightMargin = 0.04;
    thumbWidth  = 0.18;
    gapColumn   = 0.02;

    % Vertical layout
    topMargin    = 0.04;
    bottomMargin = params.bottomMargin;
    midGap       = params.midGap;
    psthHeight   = params.psthHeightFrac;

    psthBottom   = bottomMargin;
    psthTop      = bottomMargin + psthHeight;
    rasterBottom = psthTop + midGap;
    rasterTop    = 1 - topMargin;
    rasterHeight = rasterTop - rasterBottom;

    rasterWidth  = 1 - leftMargin - rightMargin - thumbWidth - gapColumn;
    rasterLeft   = leftMargin;
    thumbLeft    = rasterLeft + rasterWidth + gapColumn;

    ax_raster = axes('Position', [rasterLeft, rasterBottom, rasterWidth, rasterHeight]);
    ax_thumb  = axes('Position', [thumbLeft,  rasterBottom, thumbWidth,  rasterHeight]);
    ax_psth   = axes('Position', [rasterLeft, psthBottom,   rasterWidth, psthHeight]);

    % ==========================================================================
    % 9. RASTER PANEL
    % ==========================================================================

    axes(ax_raster);

    M = Mr2 .* (1000 / bin);
    imagesc(1:nB, 1:nRows, M);
    colormap(ax_raster, flipud(gray(64)));
    hold on;

    xline(preBase/bin,           'k', 'LineWidth', 1.5);
    xline((stimDur+preBase)/bin, 'k', 'LineWidth', 1.5);
    ticks = trialsPerCath:trialsPerCath:nT;
    yticks(ticks(1:end-1))

    xticks([]);   % Time axis is labeled on the PSTH below

    yline((rowsPerCath:rowsPerCath:nRows-1) + 0.5, 'LineWidth', 1);
    yline((nRows/nState:nRows/nState:nRows-1) + 0.5, 'LineWidth', 3, 'Color', 'k');

    % Grey patch: best image category
    patch([1, nB, nB, 1], ...
          [bestCatRowStart-0.5, bestCatRowStart-0.5, ...
           bestCatRowEnd+0.5,   bestCatRowEnd+0.5], ...
          'k', 'FaceAlpha', 0.12, 'EdgeColor', 'none');

    % Red patch: best trial x best window
    patch([best_col, best_col+nWinBins, best_col+nWinBins, best_col], ...
          [bestDisplayRow-0.5, bestDisplayRow-0.5, ...
           bestDisplayRow+0.5, bestDisplayRow+0.5], ...
          'r', 'FaceAlpha', 0.35, 'EdgeColor', 'none');

    if params.MaxVal_1
        clim([0 1]);
    else
        colorbar;
    end

    ylabel(sprintf('%d trials', nRows * mergeTrials), ...
        'FontSize', 10, 'FontName', 'helvetica');
    ax_raster.FontSize = 8;
    ax_raster.FontName = 'helvetica';

    % ==========================================================================
    % 10. THUMBNAIL PANEL
    % ==========================================================================

    axes(ax_thumb);
    imagesc(combinedImg);
    axis image off;

    % ==========================================================================
    % 11. PSTH PANEL
    % ==========================================================================

    axes(ax_psth);

    MRhist = BuildBurstMatrix(goodU(:, u), round(p.t), ...
        round(directimesSorted(trialsAbsolute) - preBase), round(win));
    MRhist = squeeze(MRhist);

    [nT2, nB2] = size(MRhist);
    spikeTimes = repmat(1:nB2, nT2, 1);
    spikeTimes = spikeTimes(logical(MRhist));

    psthBin    = params.psthBinWidth;
    edges      = 1:psthBin:round(win);
    psthCounts = histcounts(spikeTimes, edges);
    psthRate   = (psthCounts / (psthBin * nT2)) * 1000;

    b = bar(edges(1:end-1), psthRate, 'histc');
    b.FaceColor = 'k';  b.FaceAlpha = 0.3;  b.MarkerEdgeColor = 'none';
    hold on;

    xline(preBase,           'k', 'LineWidth', 1.5);
    xline(stimDur + preBase, 'k', 'LineWidth', 1.5);

    xlim([0, win]);

    try
        ylim([0, max(psthRate) + std(psthRate)]);
    catch
    end

    ylims = ylim;
    yticks([round(ylims(2)/2), round(ylims(2))]);

    xticks([0, preBase:preBase:preBase+ceil(stimDur/1000)*1000, win]);

    xticklabels(arrayfun(@(x) sprintf('%.1f', x), ...
        [-preBase, 0:preBase:ceil(stimDur/1000)*1000, stimDur+preBase] / 1000, ...
        'UniformOutput', false));

    xlabel('Time [s]', 'FontSize', 10, 'FontName', 'helvetica');
    ylabel('[spk/s]',  'FontSize', 10, 'FontName', 'helvetica');
    ax_psth.FontSize = 8;
    ax_psth.FontName = 'helvetica';

    % ==========================================================================
    % 12. RAW DATA FIGURE
    % ==========================================================================

    if params.plotRawData
        if params.fixedWindow
            rawStart  = -50;
            rawWindow = stimDur + 100;
        else
            [~, maxRespIn] = max(NeuronResp.NeuronVals(u, :, 1));
            rawStart       = NeuronResp.NeuronVals(u, maxRespIn, 3) * NeuronResp.params.binRaster - 20;
            rawWindow      = 500;
            maxRespIn_0    = maxRespIn - 1;
            trialsAbsolute = maxRespIn_0*trialsPerCath + 1 : maxRespIn_0*trialsPerCath + trialsPerCath;
            bestTrialIdx   = trialsAbsolute(1);
        end

        startTimes = directimesSorted(trialsAbsolute) + rawStart;
        spikes = squeeze(BuildBurstMatrix(goodU(:,u), round(p.t), ...
            round(startTimes), round(rawWindow)));

        if params.oneTrial
            [~, ind] = max(sum(spikes, 2));
        else
            ind = 1:size(spikes, 1);
        end

        fig2 = figure;
        fig2.Position = [147 270 662 446];

        [fig2, ~, ~] = PlotRawDataNP(obj, fig=fig2, chan=goodU(1,u), ...
            startTimes=startTimes(ind), window=rawWindow, ...
            spikeTimes=spikes(ind,:), multFactor=1.5, stdMult=3);

        xline(-rawStart/1000,               'LineWidth', 1.5, 'Color', '#77AC30');
        xline((stimDur+abs(rawStart))/1000, 'LineWidth', 1.5, 'Color', '#0072BD');

        xticks([0, abs(rawStart)/1000 : abs(rawStart)/1000 : ...
            obj.VST.stimDuration + abs(rawStart*2)/1000 + 1]);
        xticklabels([rawStart, 0:abs(rawStart):obj.VST.stimDuration*1000+abs(rawStart*2)]);
        xlabel('Milliseconds');
        yticks([]);
        ylabel('uV');
        title(sprintf('U.%d  Phy-%d  p=%.4f', u, phy_IDg(u), pvals(2,ur)));

        if params.PaperFig
            obj.printFig(fig2, sprintf('%s-NatImg-rawData-eNeuron-%d', ...
                obj.dataObj.recordingName, u), "PaperFig", true);
        elseif params.overwrite
            obj.printFig(fig2, sprintf('%s-NatImg-rawData-eNeuron-%d', ...
                obj.dataObj.recordingName, u));
        end
    end

    % ==========================================================================
    % 13. EXPORT RASTER FIGURE
    % ==========================================================================

    if params.PaperFig
        obj.printFig(fig, sprintf('%s-NatImg-raster-eNeuron-%d', ...
            obj.dataObj.recordingName, u), "PaperFig", true);
    elseif params.overwrite
        obj.printFig(fig, sprintf('%s-NatImg-raster-eNeuron-%d', ...
            obj.dataObj.recordingName, u));
    end

    ur = ur + 1;

end  % end neuron loop

end  % end plotRaster