function plotRaster_MultiExp(exList, params)

arguments
    exList double
    params.stimTypes  (1,:) string  = ["rectGrid", "linearlyMovingBall"]
    params.binWidth   double        = 10
    params.smooth     double        = 0
    params.statType   string        = "BootstrapPerNeuron"
    params.speed      string        = "max"
    params.alpha      double        = 0.05
    params.postStim   double        = 500
    params.preBase    double        = 200
    params.overwrite  logical       = false
    params.TakeTopPercentTrials double = 0.3
    params.zScore     logical       = true   % default true — more meaningful for raster
    params.sortBy     string        = "peak" % "peak" = sort by peak response time, "depth" = sort by depth
    params.PaperFig   logical       = false
    params.climPrctile double = 90   % percentile for color limit — lower = more contrast
    params.climNeg double = 0   % fixed negative z-score limit (absolute value)
     params.colormap string = "gray"
end

% -------------------------------------------------------------------------
% Load depth info if sorting by depth
% -------------------------------------------------------------------------
if params.sortBy == "depth"
    depthFile = 'W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\NeuronDepths.mat';
    if ~exist(depthFile, 'file')
        error('NeuronDepths.mat not found. Run getNeuronDepths() first.');
    end
    D = load(depthFile);
    depthTable = D.depthTable;
end

% -------------------------------------------------------------------------
% Build save path
% -------------------------------------------------------------------------
NP_first = loadNPclassFromTable(exList(1));
vs_first  = linearlyMovingBallAnalysis(NP_first);

p = extractBefore(vs_first.getAnalysisFileName, 'lizards');
p = [p 'lizards'];
if ~exist([p '\Combined_lizard_analysis'], 'dir')
    cd(p)
    mkdir Combined_lizard_analysis
end
saveDir = [p '\Combined_lizard_analysis'];

stimLabel  = strjoin(params.stimTypes, '-');
nameOfFile = sprintf('\\Ex_%d-%d_Raster_%s.mat', exList(1), exList(end), stimLabel);

% -------------------------------------------------------------------------
% Decide whether to recompute or load
% -------------------------------------------------------------------------
if exist([saveDir nameOfFile], 'file') == 2 && ~params.overwrite
    S = load([saveDir nameOfFile]);
    if isequal(S.expList, exList)
        fprintf('Loading saved raster data from:\n  %s\n', [saveDir nameOfFile]);
        forloop = false;
    else
        fprintf('Experiment list mismatch — recomputing.\n');
        forloop = true;
    end
else
    forloop = true;
end

% =========================================================================
% EXPERIMENT LOOP
% =========================================================================
if forloop

    nStim = numel(params.stimTypes);
    nExp  = numel(exList);

    % rasterAll{s} grows one row per responsive neuron across all experiments
    % each row = mean PSTH of one neuron in spk/s (or z-score)
    rasterAll  = cell(1, nStim);   % nNeurons x nBins
    depthAll   = cell(1, nStim);   % nNeurons x 1 — depth of each neuron
    expAll     = cell(1, nStim);   % nNeurons x 1 — which experiment each neuron came from

    for s = 1:nStim
        rasterAll{s} = [];
        depthAll{s}  = [];
        expAll{s}    = [];
    end

    lockedPreBase = [];
    lockedNBins   = [];
    lockedEdges   = [];
    tAxis         = [];

    for ei = 1:nExp

        ex = exList(ei);
        fprintf('\n=== Experiment %d ===\n', ex);

        try
            NP = loadNPclassFromTable(ex);
        catch ME
            warning('Could not load experiment %d: %s', ex, ME.message);
            continue
        end

        for s = 1:nStim

            stimType = params.stimTypes(s);

            try
                switch stimType
                    case "rectGrid"
                        obj = rectGridAnalysis(NP);
                    case "linearlyMovingBall"
                        obj = linearlyMovingBallAnalysis(NP);
                    case "StaticGrating"
                        obj = StaticDriftingGratingAnalysis(NP);
                    case "MovingGrating"
                        obj = StaticDriftingGratingAnalysis(NP);
                    otherwise
                        error('Unknown stimType: %s', stimType);
                end
            catch ME
                warning('Could not build %s for exp %d: %s', stimType, ex, ME.message);
                continue
            end

            NeuronResp = obj.ResponseWindow;

            if params.statType == "BootstrapPerNeuron"
                Stats = obj.BootstrapPerNeuron;
            else
                Stats = obj.ShufflingAnalysis;
            end

            % Resolve field name and stim start
            fieldName = '';
            startStim = 0;
            if params.speed ~= "max" && isequal(obj.stimName, 'linearlyMovingBall')
                fieldName = 'Speed2';
            elseif isequal(obj.stimName, 'linearlyMovingBall')
                fieldName = 'Speed1';
            elseif isequal(stimType, 'StaticGrating')
                fieldName = 'Static';
            elseif isequal(stimType, 'MovingGrating')
                fieldName = 'Moving';
                startStim = obj.VST.static_time * 1000;
            end

            p_sort = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
            label  = string(p_sort.label');
            goodU  = p_sort.ic(:, label == 'good');

            try
                pvals = Stats.(fieldName).pvalsResponse;
            catch
                pvals = Stats.pvalsResponse;
            end

            try
                C = NeuronResp.(fieldName).C;
            catch
                C = NeuronResp.C;
            end
            directimesSorted = C(:, 1)' + startStim;

            preBase     = params.preBase;
            windowTotal = preBase + params.postStim;

            if isempty(lockedPreBase)
                lockedPreBase = preBase;
                lockedEdges   = 0 : params.binWidth : windowTotal;
                lockedNBins   = numel(lockedEdges) - 1;
                tAxis         = lockedEdges(1:end-1);
                fprintf('  Locked window: preBase=%d ms, postStim=%d ms, nBins=%d\n', ...
                    lockedPreBase, params.postStim, lockedNBins);
            end

            eNeurons = find(pvals < params.alpha);

            if isempty(eNeurons)
                fprintf('  [%s] No responsive neurons in exp %d.\n', stimType, ex);
                continue
            end

            fprintf('  [%s] %d responsive neuron(s) in exp %d.\n', stimType, numel(eNeurons), ex);

            % ----------------------------------------------------------
            % Build per-neuron PSTH
            % ----------------------------------------------------------
            for ni = 1:numel(eNeurons)
                u = eNeurons(ni);

                MRhist = BuildBurstMatrix( ...
                    goodU(:, u), ...
                    round(p_sort.t), ...
                    round(directimesSorted - lockedPreBase), ...
                    round(windowTotal));
                MRhist = squeeze(MRhist);

                if ~isempty(params.TakeTopPercentTrials)
                    MeanTrial  = mean(MRhist, 2);
                    [~, ind]   = sort(MeanTrial, 'descend');
                    takeTrials = ind(1:round(numel(MeanTrial)*params.TakeTopPercentTrials));
                    MRhist     = MRhist(takeTrials, :);
                end

                nTrials    = size(MRhist, 1);
                spikeTimes = repmat((1:size(MRhist,2)), nTrials, 1);
                spikeTimes = spikeTimes(logical(MRhist));
                counts     = histcounts(spikeTimes, lockedEdges);
                neuronPSTH = (counts / (params.binWidth * nTrials)) * 1000;  % spk/s

                % Z-score using baseline
                if params.zScore
                    baselineBins = tAxis < lockedPreBase;
                    bMean = mean(neuronPSTH(baselineBins));
                    bStd  = std(neuronPSTH(baselineBins));
                    if bStd > 0
                        neuronPSTH = (neuronPSTH - bMean) / bStd;
                    else
                        continue  % skip neuron if baseline std is zero
                    end
                end

                % Smooth if requested
                if params.smooth > 0
                    smoothBins = round(params.smooth / params.binWidth);
                    neuronPSTH = smoothdata(neuronPSTH, 'gaussian', smoothBins);
                end

                % Append neuron row
                rasterAll{s} = [rasterAll{s}; neuronPSTH];

                % Get depth for this neuron
                if params.sortBy == "depth"
                    depthRow = depthTable.Experiment == ex & depthTable.Unit == u;
                    if any(depthRow)
                        depthAll{s}(end+1) = depthTable.Depth_um(depthRow);
                    else
                        depthAll{s}(end+1) = NaN;
                    end
                else
                    depthAll{s}(end+1) = NaN;
                end

                expAll{s}(end+1) = ex;
            end

        end % stim loop
    end % experiment loop

    % ------------------------------------------------------------------
    % Save
    % ------------------------------------------------------------------
    S.expList       = exList;
    S.lockedEdges   = lockedEdges;
    S.lockedPreBase = lockedPreBase;
    S.params        = params;

    for s = 1:numel(params.stimTypes)
        stimField = matlab.lang.makeValidName(params.stimTypes(s));
        S.(sprintf('%s_raster',  stimField)) = rasterAll{s};
        S.(sprintf('%s_depth',   stimField)) = depthAll{s};
        S.(sprintf('%s_exp',     stimField)) = expAll{s};
    end

    save([saveDir nameOfFile], '-struct', 'S');
    fprintf('\nSaved raster data to:\n  %s\n', [saveDir nameOfFile]);

else
    % Load from disk
    lockedEdges   = S.lockedEdges;
    lockedPreBase = S.lockedPreBase;

    rasterAll = cell(1, numel(params.stimTypes));
    depthAll  = cell(1, numel(params.stimTypes));
    expAll    = cell(1, numel(params.stimTypes));

    for s = 1:numel(params.stimTypes)
        stimField = matlab.lang.makeValidName(params.stimTypes(s));
        rasterAll{s} = S.(sprintf('%s_raster', stimField));
        depthAll{s}  = S.(sprintf('%s_depth',  stimField));
        expAll{s}    = S.(sprintf('%s_exp',    stimField));
    end
end

tAxis     = lockedEdges(1:end-1);
tAxisPlot = tAxis - lockedPreBase;

% =========================================================================
% SORT NEURONS
% =========================================================================
for s = 1:numel(params.stimTypes)
    data = rasterAll{s};
    if isempty(data); continue; end

    if params.sortBy == "peak"
        % Sort by time of peak response in the post-stimulus window
        postStimBins = tAxis >= lockedPreBase;
        [~, peakBin] = max(data(:, postStimBins), [], 2);
        [~, sortIdx] = sort(peakBin);
    elseif params.sortBy == "depth"
        % Sort by depth (shallow to deep)
        [~, sortIdx] = sort(depthAll{s}, 'ascend');
    else
        sortIdx = 1:size(data, 1);  % no sorting
    end

    rasterAll{s} = data(sortIdx, :);
    depthAll{s}  = depthAll{s}(sortIdx);
    expAll{s}    = expAll{s}(sortIdx);
end

% =========================================================================
% PLOT
% =========================================================================


stimLegendMap = containers.Map(...
    {'linearlyMovingBall', 'rectGrid', 'MovingGrating', 'StaticGrating'}, ...
    {'MB',                 'SB',       'MG',            'SG'});

nStim = numel(params.stimTypes);

% ------------------------------------------------------------------
% ------------------------------------------------------------------
% Global color limits — use lower percentile for better contrast
allValues = [];
for s = 1:nStim
    if ~isempty(rasterAll{s})
        allValues = [allValues, rasterAll{s}(:)']; %#ok<AGROW>
    end
end

if params.zScore
    cLimPos = prctile(allValues, params.climPrctile);   % data-driven positive limit
    cLims   = [-params.climNeg, cLimPos];               % asymmetric: fixed neg, data-driven pos
else
    cLims = [prctile(allValues, 2), prctile(allValues, params.climPrctile)];
end

% ------------------------------------------------------------------
% Figure and tiled layout
% ------------------------------------------------------------------
fig = figure;
set(fig, 'Units', 'centimeters', 'Position', [5 5 5*nStim + 2, 10]);

tl = tiledlayout(fig, 1, nStim, 'TileSpacing', 'compact', 'Padding', 'compact');

axAll = gobjects(1, nStim);

for s = 1:nStim

    data = rasterAll{s};
    stimKey = char(params.stimTypes(s));
    if isKey(stimLegendMap, stimKey)
        shortName = stimLegendMap(stimKey);
    else
        shortName = stimKey;
    end

    axAll(s) = nexttile(tl);
    ax = axAll(s);

    if isempty(data)
        title(ax, shortName, 'FontName', 'helvetica', 'FontSize', 8);
        axis(ax, 'off');
        continue
    end

    % imagesc: x = time, y = neuron index
    imagesc(ax, tAxisPlot, 1:size(data,1), data);
    clim(ax, cLims);
    colormap(ax, flipud(gray));  % white = low, black = high
    if params.zScore && params.colormap ~= "gray"
        cLimPos = prctile(allValues, params.climPrctile);
        cLims   = [-params.climNeg, cLimPos];

        % Proportion of colors for each side — white stays at zero
        nColors  = 256;
        nNeg     = round(nColors * params.climNeg / (params.climNeg + cLimPos));
        nPos     = nColors - nNeg;

        blueHalf = [linspace(0.1, 1, nNeg)', linspace(0.2, 1, nNeg)', linspace(0.8, 1, nNeg)'];
        redHalf  = [linspace(1, 0.9, nPos)', linspace(1, 0.2, nPos)', linspace(1, 0.05, nPos)'];
        colormap(ax, [blueHalf; redHalf]);
    else
        cLims = [prctile(allValues, 2), prctile(allValues, params.climPrctile)];
        colormap(ax, flipud(gray));
    end

    % ------------------------------------------------------------------
    % Depth bin boundary lines (only when sorted by depth)
    % ------------------------------------------------------------------
    if params.sortBy == "depth" && ~isempty(depthAll{s})

        % Load bin edges
        depthFile = 'W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\NeuronDepths.mat';
        D = load(depthFile);
        depthBinEdges = D.depthBinEdges;

        binLabelsDepth = {sprintf('%.0f-%.0f um', depthBinEdges(1), depthBinEdges(2)), ...
            sprintf('%.0f-%.0f um', depthBinEdges(2), depthBinEdges(3)), ...
            sprintf('%.0f-%.0f um', depthBinEdges(3), depthBinEdges(4))};

        % Find the last neuron index belonging to each bin boundary
        for edge = 2:3  % edges 2 and 3 are the internal boundaries
            %lastInBin = find(depthAll{s} <= depthBinEdges(edge), 1, 'last');
            %lastInBin = find(~isnan(depthAll{s}) & depthAll{s} <= depthBinEdges(edge), 1, 'last');
            depthCombined = depthAll{s};
            depthCombined = depthCombined(~isnan(depthCombined));
            lastInBin = find(depthCombined <= depthBinEdges(edge), 1, 'last');
            if ~isempty(lastInBin) && lastInBin < size(data,1)
                yline(ax, lastInBin + 0.5, 'k-', 'LineWidth', 1.2);
                % Label on the right side showing the bin range
                text(ax, tAxisPlot(5), lastInBin - size(data,1)*0.02, ...
                    binLabelsDepth{edge-1}, ...
                    'Color', 'w', 'FontSize', 6, 'FontName', 'helvetica', ...
                    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
            end
        end
        % Label for the deepest bin
        text(ax, tAxisPlot(5), size(data,1), ...
            binLabelsDepth{3}, ...
            'Color', 'w', 'FontSize', 6, 'FontName', 'helvetica', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    end

    % Stim onset and offset lines
    xline(ax, 0,               'k--', 'LineWidth', 1.0);
    xline(ax, params.postStim, 'k--', 'LineWidth', 1.0);

    xlim(ax, [tAxisPlot(1), tAxisPlot(end)]);
    ylim(ax, [0.5, size(data,1)+0.5]);
    xticks(ax, -params.preBase : 100 : params.postStim);

    xlabel(ax, 'Time re. stim onset [ms]', 'FontName', 'helvetica', 'FontSize', 8);
    if s == 1
        ylabel(ax, 'Neuron #', 'FontName', 'helvetica', 'FontSize', 8);
    end
    title(ax, sprintf('%s  (n=%d)', shortName, size(data,1)), ...
        'FontName', 'helvetica', 'FontSize', 8);

    ax.FontName = 'helvetica';
    ax.FontSize  = 8;
    ax.YDir      = 'normal';  % neuron 1 at bottom
 

end

% ------------------------------------------------------------------
% Single colorbar for the whole layout
% ------------------------------------------------------------------
cb = colorbar(axAll(end));
if params.zScore
    cb.Label.String = 'Z-score';
else
    cb.Label.String = 'Firing rate [spk/s]';
end
cb.Label.FontName = 'helvetica';
cb.Label.FontSize = 8;
cb.FontName       = 'helvetica';
cb.FontSize       = 8;

sgtitle(sprintf('N = %d experiments', numel(exList)), ...
    'FontName', 'helvetica', 'FontSize', 10);

if params.PaperFig
    vs_first.printFig(fig, sprintf('Raster-%s', stimLabel), PaperFig=params.PaperFig);
end

end