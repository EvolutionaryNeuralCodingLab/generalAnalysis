function plotRaster_MultiExp(exList, params)
% plotRaster_MultiExp  Build and display a multi-experiment raster plot.
%
%   plotRaster_MultiExp(exList, params)
%
%   For each experiment in exList the function:
%     1. Loads spike-sorted data for each stimulus type.
%     2. Identifies statistically responsive neurons.
%     3. Builds a per-neuron PSTH (optionally z-scored and smoothed).
%     4. Sorts neurons by peak-response time or recording depth.
%     5. Displays an imagesc raster for each stimulus type side-by-side,
%        with a shared x-label and a single colorbar.
%
% -------------------------------------------------------------------------
%  BUGS FIXED
% -------------------------------------------------------------------------
%   BUG 1 — Sort-by-peak double-smoothing
%     ConvBurstMatrix was applied in-place and the smoothed version stored
%     back into rasterAll, causing neurons to be smoothed twice when
%     params.smooth > 0. Fixed: peak-finding uses a local copy (dataForSort);
%     rasterAll always stores the original data.
%
%   BUG 2 — Wrong colour limits for zScore=true + colormap="gray"
%     The else-branch in the tile loop overwrote cLims with raw firing-rate
%     percentiles, ignoring climNeg. Fixed: cLims and colormap are computed
%     once before the tile loop and reused.
%
%   BUG 3 — Stim-offset xline in wrong units
%     xline used the ms value on a seconds axis. Fixed: xline now uses
%     stimDurAll(s) / 1000.

arguments
    exList double                                                   % vector of experiment IDs to include
    params.stimTypes  (1,:) string = ["rectGrid","linearlyMovingBall"] % stimulus types — one tile each
    params.binWidth   double       = 10                             % PSTH bin width in ms
    params.smooth     double       = 0                             % Gaussian smoothing SD in ms (0 = off)
    params.statType   string       = "MaxPermuteTest"              % which statistics field to use
    params.speed      string       = "max"                         % ball-speed selector: "max" or other
    params.alpha      double       = 0.05                          % significance threshold
    params.postStim   double       = 0                           % post-stimulus window in ms;
                                                                   %   when useCompleteWindow=true this is
                                                                   %   added as a post-offset buffer on top
                                                                   %   of the actual stimulus duration
    params.preBase    double       = 200                           % pre-stimulus baseline in ms
    params.overwrite  logical      = false                         % if true, recompute even if cache exists
    params.TakeTopPercentTrials double = 0.3                       % top fraction of trials to keep by mean
                                                                   %   firing rate; set empty to keep all
    params.zScore     logical      = true                          % z-score each neuron using its baseline
    params.sortBy     string       = "peak"                        % "peak" | "depth" | "none"
    params.PaperFig   logical      = false                         % if true, export figure via printFig
    params.climPrctile double      = 90                            % upper percentile for colour scale
    params.climNeg    double       = 0                             % fixed negative z-score colour limit
    params.colormap   string       = "gray"                        % "gray" -> flipud(gray); else -> diverging
    params.GaussianLength          = 5                             % Gaussian kernel half-width (bins) for sorting
    params.useCompleteWindow logical = true                       % if true, read stimulus duration from
                                                                   %   NeuronResp.(fieldName).stimDur and use
                                                                   %   params.postStim as a post-offset buffer
end

% -------------------------------------------------------------------------
% Load depth table when sorting by cortical depth
% -------------------------------------------------------------------------
if params.sortBy == "depth"
    depthFile = 'W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\NeuronDepths.mat';
    if ~exist(depthFile, 'file')                                    % abort early if depth file is missing
        error('NeuronDepths.mat not found. Run getNeuronDepths() first.');
    end
    D          = load(depthFile);                                   % load MAT file containing depth info
    depthTable = D.depthTable;                                      % table with columns: Experiment, Unit, Depth_um
end

% -------------------------------------------------------------------------
% Derive save/load path from the first experiment
% -------------------------------------------------------------------------
NP_first  = loadNPclassFromTable(exList(1));                        % load NP object for first experiment (path only)
vs_first  = linearlyMovingBallAnalysis(NP_first);                  % build analysis object to access file-system path

p = extractBefore(vs_first.getAnalysisFileName, 'lizards');        % root directory up to 'lizards' token
p = [p 'lizards'];                                                  % reconstruct path including 'lizards'
if ~exist([p '\Combined_lizard_analysis'], 'dir')                   % create output folder if absent
    cd(p)
    mkdir Combined_lizard_analysis
end
saveDir = [p '\Combined_lizard_analysis'];                          % full path to output directory

stimLabel  = strjoin(params.stimTypes, '-');                        % e.g. "rectGrid-linearlyMovingBall"
nameOfFile = sprintf('\\Ex_%d-%d_Raster_%s.mat', ...               % unique filename from experiment range + stim types
    exList(1), exList(end), stimLabel);

% -------------------------------------------------------------------------
% Decide whether to recompute or reload from cache
% -------------------------------------------------------------------------
if exist([saveDir nameOfFile], 'file') == 2 && ~params.overwrite   % cache exists and overwrite not forced
    S = load([saveDir nameOfFile]);                                 % load cached struct
    if isequal(S.expList, exList)                                   % verify cached experiment list matches
        fprintf('Loading saved raster data from:\n  %s\n', [saveDir nameOfFile]);
        forloop = false;                                            % skip computation
    else
        fprintf('Experiment list mismatch — recomputing.\n');
        forloop = true;                                             % list differs: recompute
    end
else
    forloop = true;                                                 % no cache or overwrite requested
end

% =========================================================================
% EXPERIMENT LOOP — collect responsive neurons across all experiments
% =========================================================================
if forloop

    nStim = numel(params.stimTypes);                                % number of stimulus conditions
    nExp  = numel(exList);                                          % number of experiments

    % Accumulators: one entry per stimulus type, growing one row per neuron
    rasterAll = cell(1, nStim);                                     % nNeurons x nBins PSTH matrix per stim
    depthAll  = cell(1, nStim);                                     % recording depth (um) per neuron
    expAll    = cell(1, nStim);                                     % experiment ID per neuron row

    for s = 1:nStim
        rasterAll{s} = [];
        depthAll{s}  = [];
        expAll{s}    = [];
    end

    % lockedPreBase is shared (same baseline for all stimuli).
    % Time-axis variables are per-stimulus (cell/array) because each
    % stimulus can have a different total window when useCompleteWindow=true.
    lockedPreBase = [];                                             % baseline duration (ms) — locked on first exp
    lockedEdges   = cell(1, nStim);                                 % bin edges (ms) — one set per stimulus
    lockedNBins   = zeros(1, nStim);                                % number of bins per stimulus
    tAxis         = cell(1, nStim);                                 % left bin edge (ms) per stimulus
    stimDurAll    = zeros(1, nStim);                                % stim duration (ms) per stimulus for xline

    % Pre-scan: find the shortest stimulus duration per stimulus type.
    % All experiments are truncated to this minimum so that no trial
    % window exceeds what every session can provide.
    minStimDur = inf(1, nStim);                                     % one minimum per stimulus type

    for ei = 1:nExp
        for s = 1:nStim
            try
                NPtmp = loadNPclassFromTable(exList(ei));
                switch params.stimTypes(s)
                    case "rectGrid";            objTmp = rectGridAnalysis(NPtmp);
                    case "linearlyMovingBall";  objTmp = linearlyMovingBallAnalysis(NPtmp);
                    case "StaticGrating";       objTmp = StaticDriftingGratingAnalysis(NPtmp);
                    case "MovingGrating";       objTmp = StaticDriftingGratingAnalysis(NPtmp);
                end
                NRtmp = objTmp.ResponseWindow;

                % Resolve fieldName using the same logic as the main loop
                if params.speed ~= "max" && isequal(objTmp.stimName, 'linearlyMovingBall')
                    fn = 'Speed2';
                elseif isequal(objTmp.stimName, 'linearlyMovingBall')
                    fn = 'Speed1';
                elseif isequal(params.stimTypes(s), 'StaticGrating')
                    fn = 'Static';
                elseif isequal(params.stimTypes(s), 'MovingGrating')
                    fn = 'Moving';
                else
                    fn = '';                                         % rectGrid: flat struct, no sub-field
                end

                try
                    dur = NRtmp.(fn).stimDur;                       % named sub-field (e.g. Speed1, Moving)
                catch
                    dur = NRtmp.stimDur;                            % fallback: flat struct (e.g. rectGrid)
                end

                minStimDur(s) = min(minStimDur(s), dur);            % keep shortest duration for this stimulus
            catch
                % skip quietly if experiment/stim cannot be loaded
            end
        end
    end

    fprintf('Minimum stimulus durations per type (ms):');
    for s = 1:nStim
        fprintf('  %s = %.0f ms', params.stimTypes(s), minStimDur(s));
    end
    fprintf('\n');

    % -----------------------------------------------------------------
    % Main loop: experiments x stimulus types
    % -----------------------------------------------------------------
    for ei = 1:nExp

        ex = exList(ei);
        fprintf('\n=== Experiment %d ===\n', ex);

        try
            NP = loadNPclassFromTable(ex);                          % load Neuropixels data object
        catch ME
            warning('Could not load experiment %d: %s', ex, ME.message);
            continue                                                 % skip experiment on load failure
        end

        for s = 1:nStim

            stimType = params.stimTypes(s);                         % current stimulus type string

            % Build stimulus-specific analysis object
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
                continue                                             % skip this stim/exp on failure
            end

            NeuronResp = obj.ResponseWindow;                        % full response-window struct for this stimulus

            % Select the correct statistics struct
            if params.statType == "BootstrapPerNeuron"
                Stats = obj.BootstrapPerNeuron;                     % bootstrap-based p-values
            else
                Stats = obj.StatisticsPerNeuron;                    % default: permutation-test p-values
            end

            % Resolve sub-field name and any intra-window stim onset offset
            % SUGGESTION: a switch/case or a method on each analysis class
            % would be cleaner than chained if-elseif here.
            fieldName = '';                                         % sub-field key into Stats / NeuronResp
            startStim = 0;                                          % ms offset to align stim onset to t = 0
            if params.speed ~= "max" && isequal(obj.stimName, 'linearlyMovingBall')
                fieldName = 'Speed2';                               % slower ball-speed condition
            elseif isequal(obj.stimName, 'linearlyMovingBall')
                fieldName = 'Speed1';                               % fastest (max) ball-speed condition
            elseif isequal(stimType, 'StaticGrating')
                fieldName = 'Static';                               % static grating sub-field
            elseif isequal(stimType, 'MovingGrating')
                fieldName = 'Moving';                               % moving grating sub-field
                startStim = obj.VST.static_time * 1000;            % moving phase starts after static (s -> ms)
            end

            % Convert Phy sorting output to time x unit matrix
            p_sort = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder); % Phy -> tIc format
            label  = string(p_sort.label');                         % quality label per unit
            goodU  = p_sort.ic(:, label == 'good');                 % keep only manually curated 'good' units

            % Extract response p-values; try named sub-field first
            try
                pvals = Stats.(fieldName).pvalsResponse;            % stim-specific p-values
            catch
                pvals = Stats.pvalsResponse;                        % fallback: flat struct
            end

            % Extract stimulus onset times from the condition matrix
            try
                C = NeuronResp.(fieldName).C;                       % condition matrix for this sub-field
            catch
                C = NeuronResp.C;                                   % fallback
            end
            directimesSorted = C(:, 1)' + startStim;               % stim onset times in ms

            % ----------------------------------------------------------
            % Determine the total window for this stimulus
            % ----------------------------------------------------------
            preBase = params.preBase;                               % baseline duration (ms)

            if params.useCompleteWindow
                rawStimDur_ms = minStimDur(s);                              % shortest duration for this stimulus type
                windowTotal   = preBase + rawStimDur_ms + params.postStim;  % baseline + truncated stim + post-offset buffer
            else
                rawStimDur_ms = params.postStim;                            % fixed window as before
                windowTotal   = preBase + params.postStim;
            end

            % Lock the baseline duration on the very first experiment
            if isempty(lockedPreBase)
                lockedPreBase = preBase;                            % shared across all stimuli
            end

            % Lock bin edges per stimulus on the first experiment that
            % provides them — all subsequent experiments use the same edges
            % so that every row in rasterAll{s} has identical length.
            if isempty(lockedEdges{s})
                lockedEdges{s} = 0 : params.binWidth : windowTotal; % bin edges from 0 to windowTotal (ms)
                lockedNBins(s) = numel(lockedEdges{s}) - 1;         % number of bins for this stimulus
                tAxis{s}       = lockedEdges{s}(1:end-1);           % left edge of each bin (ms)
                stimDurAll(s)  = rawStimDur_ms;                     % store stim duration for xline in plot
                fprintf('  [%s] Locked window: preBase=%d ms, stimDur=%.0f ms, nBins=%d\n', ...
                    stimType, lockedPreBase, rawStimDur_ms, lockedNBins(s));
            end

            eNeurons = find(pvals < params.alpha);                  % indices of significantly responsive neurons

            if isempty(eNeurons)
                fprintf('  [%s] No responsive neurons in exp %d.\n', stimType, ex);
                continue
            end

            fprintf('  [%s] %d responsive neuron(s) in exp %d.\n', stimType, numel(eNeurons), ex);

            % ----------------------------------------------------------
            % Build per-neuron PSTH and append to rasterAll
            % ----------------------------------------------------------
            for ni = 1:numel(eNeurons)

                u = eNeurons(ni);                                   % index of this neuron in the 'good' list

                % BuildBurstMatrix returns a trials x time-samples binary matrix (1 ms resolution)
                MRhist = BuildBurstMatrix( ...
                    goodU(:, u), ...                                % binary spike train for this unit
                    round(p_sort.t), ...                            % sample timestamps (rounded to avoid float errors)
                    round(directimesSorted - lockedPreBase), ...    % trial start = stim onset minus baseline
                    round(windowTotal));                            % window duration (ms, integer)
                MRhist = squeeze(MRhist);                           % collapse singleton dimensions

                % Optionally restrict to the highest-firing trials
                if ~isempty(params.TakeTopPercentTrials)
                    MeanTrial  = mean(MRhist, 2);                   % mean spike count per trial (full window)
                    % SUGGESTION: computing the mean over the post-stim
                    % window only (columns where tAxis{s} >= lockedPreBase)
                    % would avoid selecting high-baseline trials over
                    % high-response trials.
                    [~, ind]   = sort(MeanTrial, 'descend');        % rank trials by mean activity
                    takeTrials = ind(1 : round(numel(MeanTrial) * params.TakeTopPercentTrials));
                    MRhist     = MRhist(takeTrials, :);             % keep top fraction of trials
                end

                nTrials    = size(MRhist, 1);                       % number of trials used
                spikeTimes = repmat((1:size(MRhist,2)), nTrials, 1); % bin index for every position in MRhist
                spikeTimes = spikeTimes(logical(MRhist));           % keep only bins where spikes occurred
                counts     = histcounts(spikeTimes, lockedEdges{s}); % spike count per bin (per-stimulus edges)
                neuronPSTH = (counts / (params.binWidth * nTrials)) * 1000; % convert counts -> spk/s

                % Z-score using the pre-stimulus baseline
                if params.zScore
                    baselineBins = tAxis{s} < lockedPreBase;        % logical mask for baseline bins (per-stim tAxis)
                    bMean = mean(neuronPSTH(baselineBins));         % baseline mean firing rate
                    bStd  = std(neuronPSTH(baselineBins));          % baseline standard deviation
                    if bStd > 0
                        neuronPSTH = (neuronPSTH - bMean) / bStd;  % z-score
                    else
                        % SUGGESTION: silently dropping zero-SD neurons
                        % biases the population toward cells with measurable
                        % spontaneous activity. Consider logging dropped
                        % units or using a minimum-SD floor.
                        continue                                     % skip: silent baseline -> undefined z-score
                    end
                end

                % Smooth PSTH if requested
                if params.smooth > 0
                    smoothBins = round(params.smooth / params.binWidth); % smoothing SD in bins
                    neuronPSTH = smoothdata(neuronPSTH, 'gaussian', smoothBins); % Gaussian smooth
                end

                rasterAll{s} = [rasterAll{s}; neuronPSTH];         % append neuron as a new row

                % Store recording depth (needed only for depth-sorted plots)
                if params.sortBy == "depth"
                    depthRow = depthTable.Experiment == ex & depthTable.Unit == u;
                    if any(depthRow)
                        depthAll{s}(end+1) = depthTable.Depth_um(depthRow); % depth in um
                    else
                        depthAll{s}(end+1) = NaN;                  % not found in depth table
                    end
                else
                    depthAll{s}(end+1) = NaN;                      % depth unused; keeps vectors aligned
                end

                expAll{s}(end+1) = ex;                             % record source experiment

            end % neuron loop

        end % stimulus loop
    end % experiment loop

    % ------------------------------------------------------------------
    % Save processed data to disk
    % ------------------------------------------------------------------
    S.expList       = exList;                                       % saved for validation on reload
    S.lockedEdges   = lockedEdges;                                  % cell array of per-stimulus bin edges
    S.lockedPreBase = lockedPreBase;                                % shared baseline duration
    S.stimDurAll    = stimDurAll;                                   % per-stimulus stim duration for xline
    S.params        = params;                                       % full parameter set alongside data

    for s = 1:numel(params.stimTypes)
        stimField = matlab.lang.makeValidName(params.stimTypes(s)); % valid MATLAB struct field name
        S.(sprintf('%s_raster', stimField)) = rasterAll{s};
        S.(sprintf('%s_depth',  stimField)) = depthAll{s};
        S.(sprintf('%s_exp',    stimField)) = expAll{s};
    end

    save([saveDir nameOfFile], '-struct', 'S');
    fprintf('\nSaved raster data to:\n  %s\n', [saveDir nameOfFile]);

else
    % ------------------------------------------------------------------
    % Reload cached data from disk
    % ------------------------------------------------------------------
    lockedEdges   = S.lockedEdges;                                  % restore per-stimulus bin edges
    lockedPreBase = S.lockedPreBase;                                % restore baseline duration
    stimDurAll    = S.stimDurAll;                                   % restore per-stimulus stim durations

    rasterAll = cell(1, numel(params.stimTypes));
    depthAll  = cell(1, numel(params.stimTypes));
    expAll    = cell(1, numel(params.stimTypes));

    for s = 1:numel(params.stimTypes)
        stimField    = matlab.lang.makeValidName(params.stimTypes(s));
        rasterAll{s} = S.(sprintf('%s_raster', stimField));
        depthAll{s}  = S.(sprintf('%s_depth',  stimField));
        expAll{s}    = S.(sprintf('%s_exp',    stimField));
    end

    % Reconstruct per-stimulus tAxis from the stored edges
    tAxis = cell(1, numel(params.stimTypes));
    for s = 1:numel(params.stimTypes)
        tAxis{s} = lockedEdges{s}(1:end-1);                        % left edge of each bin (ms)
    end

end

% =========================================================================
% SORT NEURONS
% =========================================================================
for s = 1:numel(params.stimTypes)

    data = rasterAll{s};                                            % nNeurons x nBins matrix
    if isempty(data); continue; end

    if params.sortBy == "peak"
        postStimBins = tAxis{s} >= lockedPreBase;                   % post-stim mask using per-stimulus tAxis

        % BUG FIX: previously ConvBurstMatrix overwrote 'data' and the
        % smoothed version was stored back into rasterAll (double-smoothing).
        % dataForSort is a local copy used only for peak detection.

        if size(data,2) > 100
            dataForSort  = ConvBurstMatrix( ...
            data, fspecial('gaussian', [1 params.GaussianLength+20], 2), 'same'); % smooth copy for peak detection only

        else
            dataForSort  = ConvBurstMatrix( ...
            data, fspecial('gaussian', [1 params.GaussianLength], 2), 'same'); % smooth copy for peak detection only
        end

        [~, peakBin] = max(dataForSort(:, postStimBins), [], 2);   % column of peak per neuron
        [~, sortIdx] = sort(peakBin);                              % ascending: early-peaking neurons first

    elseif params.sortBy == "depth"
        [~, sortIdx] = sort(depthAll{s}, 'ascend');                % shallowest recording sites first

    else
        sortIdx = 1:size(data, 1);                                 % no reordering
    end

    rasterAll{s} = data(sortIdx, :);                               % reorder rows of the ORIGINAL data
    depthAll{s}  = depthAll{s}(sortIdx);                           % reorder depth vector to match
    expAll{s}    = expAll{s}(sortIdx);                             % reorder experiment-ID vector to match

end

% =========================================================================
% PLOT
% =========================================================================

stimLegendMap = containers.Map( ...
    {'linearlyMovingBall', 'rectGrid', 'MovingGrating', 'StaticGrating'}, ...
    {'MB',                 'SB',       'MG',            'SG'});     % short display labels

nStim = numel(params.stimTypes);

% ------------------------------------------------------------------
% Global colour limits — computed once and shared across all tiles
% so the colour scale is directly comparable between stimuli.
% BUG FIX: previously cLims was recomputed incorrectly inside the loop.
% ------------------------------------------------------------------
allValues = [];
for s = 1:nStim
    if ~isempty(rasterAll{s})
        allValues = [allValues, rasterAll{s}(:)']; %#ok<AGROW>     % pool all data values for percentile calculation
    end
end

if params.zScore
    cLimPos = prctile(allValues, params.climPrctile);               % data-driven upper z-score limit
    cLims   = [-params.climNeg, cLimPos];                           % asymmetric: fixed lower, data-driven upper
else
    cLims = [prctile(allValues, 2), prctile(allValues, params.climPrctile)]; % symmetric percentile range
end

% ------------------------------------------------------------------
% Build diverging colormap once (used when colormap ~= "gray")
% ------------------------------------------------------------------
if params.zScore && params.colormap ~= "gray"
    nColors   = 256;                                                % total colour table entries
    nNeg      = round(nColors * params.climNeg / (params.climNeg + cLimPos)); % entries for negative half
    nPos      = nColors - nNeg;                                     % entries for positive half
    blueHalf  = [linspace(0.1,1,nNeg)', linspace(0.2,1,nNeg)', linspace(0.8,1,nNeg)']; % blue -> white
    redHalf   = [linspace(1,0.9,nPos)', linspace(1,0.2,nPos)', linspace(1,0.05,nPos)']; % white -> red
    cmapToUse = [blueHalf; redHalf];                                % full diverging colormap
else
    cmapToUse = flipud(gray);                                       % default: white = low, black = high
end

% ------------------------------------------------------------------
% Figure and tiled layout
% ------------------------------------------------------------------
fig = figure;
set(fig, 'Units', 'centimeters', 'Position', [5 5 5*nStim + 2, 10]); % width scales with number of tiles

tl    = tiledlayout(fig, 1, nStim, 'TileSpacing', 'compact', 'Padding', 'compact'); % 1-row tile grid
axAll = gobjects(1, nStim);                                         % pre-allocate axes handles

for s = 1:nStim

    data    = rasterAll{s};                                         % nNeurons x nBins for this stimulus
    stimKey = char(params.stimTypes(s));                            % char for containers.Map lookup
    if isKey(stimLegendMap, stimKey)
        shortName = stimLegendMap(stimKey);                         % abbreviated label
    else
        shortName = stimKey;                                        % fallback to full name
    end

    axAll(s) = nexttile(tl);                                       % create tile and capture axes handle
    ax = axAll(s);

    if isempty(data)
        title(ax, shortName, 'FontName', 'helvetica', 'FontSize', 8);
        axis(ax, 'off');                                            % nothing to plot: hide axes
        continue
    end

    % Per-stimulus time axis in seconds
    tAxisPlot = tAxis{s} - lockedPreBase;                          % shift so stim onset = 0 ms
    tAxisSec  = tAxisPlot / 1000;                                  % convert to seconds

    imagesc(ax, tAxisSec, 1:size(data,1), data);                   % raster: x = seconds, y = neuron index
    clim(ax, cLims);                                               % shared colour limits
    colormap(ax, cmapToUse);                                       % shared colormap

    % ------------------------------------------------------------------
    % Depth-bin boundary lines (only when sortBy = "depth")
    % ------------------------------------------------------------------
    if params.sortBy == "depth" && ~isempty(depthAll{s})

        D2 = load(depthFile);                                       % load bin edges (file already validated above)
        depthBinEdges = D2.depthBinEdges;                           % edges defining depth layers (um)

        binLabelsDepth = { ...
            sprintf('%.0f-%.0f um', depthBinEdges(1), depthBinEdges(2)), ...
            sprintf('%.0f-%.0f um', depthBinEdges(2), depthBinEdges(3)), ...
            sprintf('%.0f-%.0f um', depthBinEdges(3), depthBinEdges(4))};

        depthCombined = depthAll{s}(~isnan(depthAll{s}));           % exclude NaN before boundary search
        labelX = tAxisSec(1) + 0.05 * range(tAxisSec);             % x position for labels: 5% from left edge

        for edge = 2:3                                              % internal boundaries only
            lastInBin = find(depthCombined <= depthBinEdges(edge), 1, 'last');
            if ~isempty(lastInBin) && lastInBin < size(data,1)
                yline(ax, lastInBin + 0.5, 'k-', 'LineWidth', 1.2); % horizontal separator between depth bins
                text(ax, labelX, lastInBin - size(data,1)*0.02, ...
                    binLabelsDepth{edge-1}, ...
                    'Color', 'w', 'FontSize', 6, 'FontName', 'helvetica', ...
                    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
            end
        end
        text(ax, labelX, size(data,1), binLabelsDepth{3}, ...       % label for the deepest bin
            'Color', 'w', 'FontSize', 6, 'FontName', 'helvetica', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    end

    % Stim onset and offset lines in seconds
    xline(ax, 0,                    'k--', 'LineWidth', 1.0);      % stim onset at t = 0 s
    xline(ax, stimDurAll(s)/1000,   'k--', 'LineWidth', 1.0);      % stim offset — per-stimulus, in seconds

    % Per-tile x-axis range — reflects this stimulus's own duration
    xlim(ax, [tAxisSec(1), tAxisSec(end)]);
    ylim(ax, [0.5, size(data,1) + 0.5]);                           % half-row margin above and below

    % Equal-spaced ticks at interval =  seconds.
    ticksSec    = linspace(tAxisSec(1), tAxisSec(end), 5);  % 5 perfectly equally spaced ticks
    [~, iz]     = min(abs(ticksSec));                        % find the tick nearest to 0
    ticksSec(iz) = 0;                                        % snap it to exactly 0 for a clean label
    xticks(ax, ticksSec);
    xticklabels(ax, arrayfun(@(v) sprintf('%.2g', v), ticksSec, 'UniformOutput', false));

    if s == 1
        ylabel(ax, 'Neuron #', 'FontName', 'helvetica', 'FontSize', 8); % y-label on leftmost tile only
    end

    title(ax, sprintf('%s  (n=%d)', shortName, size(data,1)), ...
        'FontName', 'helvetica', 'FontSize', 8);

    ax.FontName = 'helvetica';
    ax.FontSize  = 8;
    ax.YDir      = 'normal';                                        % neuron 1 at bottom, increasing upward
    ax.TickDir   = 'out';                                           % outward ticks (publication convention)
    ax.Box       = 'off';                                           % remove top/right border

end

% ------------------------------------------------------------------
% Shared x-label via tiledlayout — one label centred below all tiles
% ------------------------------------------------------------------
xlabel(tl, 'Time relative to stimulus onset (s)', ...
    'FontName', 'helvetica', 'FontSize', 8);

% ------------------------------------------------------------------
% Single colorbar on the rightmost tile
% ------------------------------------------------------------------
cb = colorbar(axAll(end));                                         % attach colorbar to last axes
if params.zScore
    cb.Label.String = 'Z-score';
else
    cb.Label.String = 'Firing rate (spk/s)';
end
cb.Label.FontName = 'helvetica';
cb.Label.FontSize = 8;
cb.FontName       = 'helvetica';
cb.FontSize       = 8;
set(fig, 'Units', 'centimeters', 'Position', [20 20 9 12]);

sgtitle(sprintf('N = %d experiments', numel(exList)), ...
    'FontName', 'helvetica', 'FontSize', 10);                      % super-title above the entire figure

if params.PaperFig
    vs_first.printFig(fig, sprintf('Raster-%s', stimLabel), PaperFig=params.PaperFig); % export to file
end

end