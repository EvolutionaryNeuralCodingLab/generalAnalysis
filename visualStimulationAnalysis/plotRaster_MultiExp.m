function plotRaster_MultiExp(exList, params)
% plotRaster_MultiExp  Build and display a multi-experiment raster plot.
%
%   plotRaster_MultiExp(exList, params)
%
%   For each experiment in exList the function:
%     1. Loads spike-sorted data for each stimulus type.
%     2. Identifies statistically responsive neurons.
%     3. Builds a per-neuron PSTH (optionally z-scored and smoothed).
%     4. Sorts neurons by peak-response time, recording depth,
%        spatial-tuning index, or leaves them unsorted.
%     5. Displays an imagesc raster for each stimulus type side-by-side,
%        with a shared x-label and a single colorbar.
%
% -------------------------------------------------------------------------
%  CHANGE LOG
% -------------------------------------------------------------------------
%   BUG 1 (fixed earlier) — Sort-by-peak double-smoothing.
%   BUG 2 (fixed earlier) — Wrong colour limits for zScore + gray.
%   BUG 3 (fixed earlier) — Stim-offset xline in wrong units.
%   FIX 4 — Trial selection now uses post-onset window only, avoiding
%           bias toward high-baseline trials.
%   FIX 5 — Zero-SD neurons are counted and logged instead of silently
%           dropped.
%   FIX 6 — TakeTopPercentTrials = 0 is handled explicitly.
%   FIX 7 — Baseline/binWidth alignment assertion added.
%   FIX 8 — Diverging colormap warns when climNeg = 0.
%   FIX 9 — Accumulator alignment assertion after experiment loop.
%   NEW  — sortBy = "spatialTuning": sort neurons by a column from an
%          external spatial-tuning table, matched via Phy cluster ID.
%   NEW  — phyAll accumulator tracks Phy cluster IDs for each neuron row.

arguments
    exList double                                                        % vector of experiment IDs to include
    params.stimTypes       (1,:) string  = ["rectGrid","linearlyMovingBall"] % stimulus types — one tile each
    params.binWidth        double        = 10                            % PSTH bin width in ms
    params.smooth          double        = 0                             % Gaussian smoothing SD in ms (0 = off)
    params.statType        string        = "MaxPermuteTest"              % which statistics field to use
    params.speed           string        = "max"                         % ball-speed selector: "max" or other
    params.alpha           double        = 0.05                          % significance threshold
    params.postStim        double        = 0                             % post-stimulus window in ms
    params.preBase         double        = 200                           % pre-stimulus baseline in ms
    params.overwrite       logical       = false                         % if true, recompute even if cache exists
    params.TakeTopPercentTrials double   = 0.3                           % fraction (0,1] of trials to keep; [] or 0 = keep all
    params.zScore          logical       = true                          % z-score each neuron using its baseline
    params.sortBy          string        = "spatialTuning"                        % "peak" | "depth" | "spatialTuning" | "none"
    params.PaperFig        logical       = false                         % if true, export figure via printFig
    params.climPrctile     double        = 90                            % upper percentile for colour scale
    params.climNeg         double        = 0                             % fixed negative z-score colour limit
    params.colormap        string        = "gray"                        % "gray" -> flipud(gray); else -> diverging
    params.GaussianLength                = 5                             % Gaussian kernel half-width (bins) for sort smoothing
    params.useCompleteWindow logical     = true                          % if true, use actual stim duration from data
    % --- Spatial-tuning sort parameters ---
    params.tuningIndexCol  string        = "L_amplitude_diff"            % column in tuning table to sort by
    params.tuningSortOrder string        = "descend"                     % "descend" = most-tuned at top
    params.tuningFile      string        = ""                            % full path to tuning .mat; "" = auto-construct
end

% -------------------------------------------------------------------------
% Sanity check: baseline must be an integer multiple of bin width
% -------------------------------------------------------------------------
assert(mod(params.preBase, params.binWidth) == 0, ...
    'preBase (%g ms) must be a multiple of binWidth (%g ms).', ...
    params.preBase, params.binWidth);

% -------------------------------------------------------------------------
% Load depth table when sorting by cortical depth
% -------------------------------------------------------------------------
if params.sortBy == "depth"
    depthFile = 'W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\NeuronDepths.mat';
    if ~exist(depthFile, 'file')                                         % abort if file is missing
        error('NeuronDepths.mat not found. Run getNeuronDepths() first.');
    end
    D          = load(depthFile);                                        % load depth struct
    depthTable = D.depthTable;                                           % table: Experiment, Unit, Depth_um
end

% -------------------------------------------------------------------------
% Derive save/load path from the first experiment
% -------------------------------------------------------------------------
NP_first = loadNPclassFromTable(exList(1));                              % load NP object for path info
vs_first = linearlyMovingBallAnalysis(NP_first);                         % analysis object for filesystem path

p = extractBefore(vs_first.getAnalysisFileName, 'lizards');              % root up to 'lizards' token
p = [p 'lizards'];                                                       % include 'lizards' folder

if ~exist([p '\Combined_lizard_analysis'], 'dir')                        % create output dir if absent
    cd(p)
    mkdir Combined_lizard_analysis
end
saveDir = [p '\Combined_lizard_analysis'];                               % output directory

stimLabel  = strjoin(params.stimTypes, '-');                              % e.g. "rectGrid-linearlyMovingBall"
nameOfFile = sprintf('\\Ex_%d-%d_Raster_%s.mat', ...
    exList(1), exList(end), stimLabel);                                   % cache filename

% -------------------------------------------------------------------------
% Decide whether to recompute or reload from cache
% -------------------------------------------------------------------------
if exist([saveDir nameOfFile], 'file') == 2 && ~params.overwrite         % cache exists and overwrite not forced
    S = load([saveDir nameOfFile]);                                      % load cached struct
    if isequal(S.expList, exList)                                        % cached list matches current request
        fprintf('Loading saved raster data from:\n  %s\n', [saveDir nameOfFile]);
        forloop = false;                                                 % skip computation
    else
        fprintf('Experiment list mismatch — recomputing.\n');
        forloop = true;                                                  % mismatch: recompute
    end
else
    forloop = true;                                                      % no cache or overwrite requested
end

% =========================================================================
% EXPERIMENT LOOP — collect responsive neurons across all experiments
% =========================================================================
if forloop

    nStim = numel(params.stimTypes);                                     % number of stimulus conditions
    nExp  = numel(exList);                                               % number of experiments

    % --- Accumulators: one cell per stimulus type, one row per neuron ---
    rasterAll = cell(1, nStim);                                          % nNeurons x nBins PSTH per stim
    depthAll  = cell(1, nStim);                                          % recording depth (um) per neuron
    expAll    = cell(1, nStim);                                          % experiment ID per neuron row
    phyAll    = cell(1, nStim);                                          % Phy cluster ID per neuron row

    for s = 1:nStim
        rasterAll{s} = [];                                               % initialise empty
        depthAll{s}  = [];
        expAll{s}    = [];
        phyAll{s}    = [];
    end

    % Counter for neurons dropped due to zero-SD baseline
    nDroppedZeroSD = zeros(1, nStim);                                    % per-stim counter

    % --- Shared time-axis variables ---
    lockedPreBase = [];                                                  % baseline duration (ms) — locked on first exp
    lockedEdges   = cell(1, nStim);                                      % bin edges (ms) per stimulus
    lockedNBins   = zeros(1, nStim);                                     % bin count per stimulus
    tAxis         = cell(1, nStim);                                      % left bin edge (ms) per stimulus
    stimDurAll    = zeros(1, nStim);                                     % stim duration (ms) per stim for xline

    % -----------------------------------------------------------------
    % Pre-scan: find shortest stimulus duration per type across all exps
    % -----------------------------------------------------------------
    minStimDur = inf(1, nStim);                                          % initialise with inf

    for ei = 1:nExp
        for s = 1:nStim
            try
                NPtmp = loadNPclassFromTable(exList(ei));                % load NP object
                switch params.stimTypes(s)
                    case "rectGrid";            objTmp = rectGridAnalysis(NPtmp);
                    case "linearlyMovingBall";  objTmp = linearlyMovingBallAnalysis(NPtmp);
                    case "StaticGrating";       objTmp = StaticDriftingGratingAnalysis(NPtmp);
                    case "MovingGrating";       objTmp = StaticDriftingGratingAnalysis(NPtmp);
                end
                NRtmp = objTmp.ResponseWindow;                           % response-window struct

                % Resolve fieldName (same logic as main loop)
                if params.speed ~= "max" && isequal(objTmp.stimName, 'linearlyMovingBall')
                    fn = 'Speed2';
                elseif isequal(objTmp.stimName, 'linearlyMovingBall')
                    fn = 'Speed1';
                elseif isequal(params.stimTypes(s), 'StaticGrating')
                    fn = 'Static';
                elseif isequal(params.stimTypes(s), 'MovingGrating')
                    fn = 'Moving';
                else
                    fn = '';                                             % rectGrid: flat struct
                end

                try
                    dur = NRtmp.(fn).stimDur;                            % sub-field duration
                catch
                    dur = NRtmp.stimDur;                                 % flat struct fallback
                end

                minStimDur(s) = min(minStimDur(s), dur);                 % keep shortest for this stim
            catch
                % skip quietly if experiment/stim cannot be loaded
            end
        end
    end

    fprintf('Minimum stimulus durations per type (ms):');
    for s = 1:nStim
        fprintf('  %s = %.0f ms', params.stimTypes(s), minStimDur(s));   % display per-stim duration
    end
    fprintf('\n');

    % -----------------------------------------------------------------
    % Main loop: experiments x stimulus types
    % -----------------------------------------------------------------
    for ei = 1:nExp

        ex = exList(ei);                                                 % current experiment ID
        fprintf('\n=== Experiment %d ===\n', ex);

        try
            NP = loadNPclassFromTable(ex);                               % load Neuropixels data object
        catch ME
            warning('Could not load experiment %d: %s', ex, ME.message);
            continue                                                      % skip on load failure
        end

        for s = 1:nStim

            stimType = params.stimTypes(s);                              % current stimulus string

            % --- Build stimulus-specific analysis object ---
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
                continue                                                  % skip this stim/exp
            end

            NeuronResp = obj.ResponseWindow;                             % response-window struct

            % --- Select statistics struct ---
            if params.statType == "BootstrapPerNeuron"
                Stats = obj.BootstrapPerNeuron;                          % bootstrap p-values
            else
                Stats = obj.StatisticsPerNeuron;                         % default: permutation test
            end

            % --- Resolve sub-field name and stim-onset offset ---
            fieldName = '';                                              % sub-field key
            startStim = 0;                                               % ms offset for stim onset
            if params.speed ~= "max" && isequal(obj.stimName, 'linearlyMovingBall')
                fieldName = 'Speed2';                                    % slower speed condition
            elseif isequal(obj.stimName, 'linearlyMovingBall')
                fieldName = 'Speed1';                                    % fastest speed condition
            elseif isequal(stimType, 'StaticGrating')
                fieldName = 'Static';                                    % static grating sub-field
            elseif isequal(stimType, 'MovingGrating')
                fieldName = 'Moving';                                    % moving grating sub-field
                startStim = obj.VST.static_time * 1000;                 % moving phase onset (s -> ms)
            end

            % --- Convert Phy sorting to tIc format ---
            p_sort = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
            label  = string(p_sort.label');                              % quality label per unit
            goodU  = p_sort.ic(:, label == 'good');                      % keep only curated 'good' units
        
            % --- Extract Phy cluster IDs for good units ---
            goodPhyIDs = p_sort.phy_ID(label == 'good');                 % Phy cluster IDs matching goodU columns                                 % Phy cluster ID per good unit

            % --- Response p-values ---
            try
                pvals = Stats.(fieldName).pvalsResponse;                 % stim-specific
            catch
                pvals = Stats.pvalsResponse;                             % flat struct fallback
            end

            % --- Stimulus onset times ---
            try
                C = NeuronResp.(fieldName).C;                            % condition matrix
            catch
                C = NeuronResp.C;                                        % fallback
            end
            directimesSorted = C(:, 1)' + startStim;                    % onset times in ms

            % --- Determine total trial window ---
            preBase = params.preBase;                                    % baseline in ms

            if params.useCompleteWindow
                rawStimDur_ms = minStimDur(s);                           % truncate to shortest duration
                windowTotal   = preBase + rawStimDur_ms + params.postStim;
            else
                rawStimDur_ms = params.postStim;                         % fixed window
                windowTotal   = preBase + params.postStim;
            end

            % Lock baseline on first experiment
            if isempty(lockedPreBase)
                lockedPreBase = preBase;                                 % shared across all stimuli
            end

            % Lock bin edges per stim on first encounter
            if isempty(lockedEdges{s})
                lockedEdges{s} = 0 : params.binWidth : windowTotal;      % bin edges from 0 to windowTotal
                lockedNBins(s) = numel(lockedEdges{s}) - 1;              % number of bins
                tAxis{s}       = lockedEdges{s}(1:end-1);                % left edge per bin (ms)
                stimDurAll(s)  = rawStimDur_ms;                          % for xline in plot
                fprintf('  [%s] Locked window: preBase=%d ms, stimDur=%.0f ms, nBins=%d\n', ...
                    stimType, lockedPreBase, rawStimDur_ms, lockedNBins(s));
            end

            % --- Find responsive neurons ---
            eNeurons = find(pvals < params.alpha);                       % indices of significant neurons

            if isempty(eNeurons)
                fprintf('  [%s] No responsive neurons in exp %d.\n', stimType, ex);
                continue
            end

            fprintf('  [%s] %d responsive neuron(s) in exp %d.\n', ...
                stimType, numel(eNeurons), ex);

            % ----------------------------------------------------------
            % Per-neuron PSTH
            % ----------------------------------------------------------
            for ni = 1:numel(eNeurons)

                u = eNeurons(ni);                                        % index into the good-unit list

                % Binary spike matrix: trials x time at 1 ms resolution
                MRhist = BuildBurstMatrix( ...
                    goodU(:, u), ...                                     % spike identity for this unit
                    round(p_sort.t), ...                                 % rounded sample timestamps
                    round(directimesSorted - lockedPreBase), ...         % trial start = onset - baseline
                    round(windowTotal));                                 % window length (ms, integer)
                MRhist = squeeze(MRhist);                                % remove singleton dims

                % --- Optionally keep only the highest-firing trials ---
                if ~isempty(params.TakeTopPercentTrials) && ...
                        params.TakeTopPercentTrials > 0 && ...
                        params.TakeTopPercentTrials < 1                  % FIX 6: guard against 0 and >=1

                    % FIX 4: rank trials by post-onset activity only,
                    % avoiding bias toward high-baseline trials.
                    postOnsetCols = (lockedPreBase + 1) : size(MRhist, 2); % columns after stim onset (1 ms res)
                    MeanTrial     = mean(MRhist(:, postOnsetCols), 2);     % mean post-onset spike count per trial
                    [~, ind]      = sort(MeanTrial, 'descend');            % rank descending
                    nKeep         = max(1, round(numel(MeanTrial) * params.TakeTopPercentTrials)); % at least 1
                    takeTrials    = ind(1:nKeep);                          % indices of top trials
                    MRhist        = MRhist(takeTrials, :);                 % keep only top fraction
                end

                nTrials    = size(MRhist, 1);                            % number of trials used
                spikeTimes = repmat((1:size(MRhist,2)), nTrials, 1);     % column index for every position
                spikeTimes = spikeTimes(logical(MRhist));                % keep only spike positions
                counts     = histcounts(spikeTimes, lockedEdges{s});     % spike count per bin
                neuronPSTH = (counts / (params.binWidth * nTrials)) * 1000; % convert to spk/s

                % --- Z-score using pre-stimulus baseline ---
                if params.zScore
                    baselineBins = tAxis{s} < lockedPreBase;             % mask for baseline bins
                    bMean = mean(neuronPSTH(baselineBins));              % baseline mean
                    bStd  = std(neuronPSTH(baselineBins));               % baseline SD
                    if bStd > 0
                        neuronPSTH = (neuronPSTH - bMean) / bStd;       % z-score
                    else
                        % FIX 5: count and log rather than silently skip
                        nDroppedZeroSD(s) = nDroppedZeroSD(s) + 1;
                        continue                                          % skip: undefined z-score
                    end
                end

                % --- Smooth PSTH if requested ---
                if params.smooth > 0
                    smoothBins = round(params.smooth / params.binWidth); % smoothing SD in bin units
                    neuronPSTH = smoothdata(neuronPSTH, 'gaussian', smoothBins);
                end

                % --- Append neuron row to accumulators ---
                rasterAll{s} = [rasterAll{s}; neuronPSTH];              % PSTH row
                phyAll{s}(end+1)  = goodPhyIDs(u);                      % Phy cluster ID for this unit
                expAll{s}(end+1)  = ex;                                  % source experiment

                % Store recording depth
                if params.sortBy == "depth"
                    depthRow = depthTable.Experiment == ex & depthTable.Unit == u;
                    if any(depthRow)
                        depthAll{s}(end+1) = depthTable.Depth_um(depthRow);
                    else
                        depthAll{s}(end+1) = NaN;                       % not found
                    end
                else
                    depthAll{s}(end+1) = NaN;                            % unused; keeps vector aligned
                end

            end % neuron loop

        end % stimulus loop
    end % experiment loop

    % FIX 5: report zero-SD dropped neurons
    for s = 1:nStim
        if nDroppedZeroSD(s) > 0
            fprintf('  [%s] Dropped %d neuron(s) with zero baseline SD.\n', ...
                params.stimTypes(s), nDroppedZeroSD(s));
        end
    end

    % FIX 9: verify accumulator alignment after all experiments
    for s = 1:nStim
        nRows = size(rasterAll{s}, 1);
        assert(numel(expAll{s})   == nRows, 'expAll{%d} length mismatch.', s);
        assert(numel(phyAll{s})   == nRows, 'phyAll{%d} length mismatch.', s);
        assert(numel(depthAll{s}) == nRows, 'depthAll{%d} length mismatch.', s);
    end

    % ------------------------------------------------------------------
    % Save processed data to disk
    % ------------------------------------------------------------------
    S.expList       = exList;                                            % for validation on reload
    S.lockedEdges   = lockedEdges;                                       % per-stim bin edges
    S.lockedPreBase = lockedPreBase;                                     % shared baseline
    S.stimDurAll    = stimDurAll;                                        % per-stim duration for xline
    S.params        = params;                                            % full parameter set

    for s = 1:numel(params.stimTypes)
        stimField = matlab.lang.makeValidName(params.stimTypes(s));      % valid struct field name
        S.(sprintf('%s_raster', stimField)) = rasterAll{s};
        S.(sprintf('%s_depth',  stimField)) = depthAll{s};
        S.(sprintf('%s_exp',    stimField)) = expAll{s};
        S.(sprintf('%s_phy',    stimField)) = phyAll{s};                 % NEW: Phy cluster IDs
    end

    save([saveDir nameOfFile], '-struct', 'S');
    fprintf('\nSaved raster data to:\n  %s\n', [saveDir nameOfFile]);

else
    % ------------------------------------------------------------------
    % Reload cached data
    % ------------------------------------------------------------------
    lockedEdges   = S.lockedEdges;                                       % restore bin edges
    lockedPreBase = S.lockedPreBase;                                     % restore baseline
    stimDurAll    = S.stimDurAll;                                        % restore stim durations

    rasterAll = cell(1, numel(params.stimTypes));
    depthAll  = cell(1, numel(params.stimTypes));
    expAll    = cell(1, numel(params.stimTypes));
    phyAll    = cell(1, numel(params.stimTypes));

    for s = 1:numel(params.stimTypes)
        stimField    = matlab.lang.makeValidName(params.stimTypes(s));
        rasterAll{s} = S.(sprintf('%s_raster', stimField));
        depthAll{s}  = S.(sprintf('%s_depth',  stimField));
        expAll{s}    = S.(sprintf('%s_exp',    stimField));

        % phyAll may be absent in old caches — require recompute if needed
        phyField = sprintf('%s_phy', stimField);
        if isfield(S, phyField)
            phyAll{s} = S.(phyField);                                    % restore Phy IDs
        elseif params.sortBy == "spatialTuning"
            error(['Cache file lacks Phy IDs (old format). ' ...
                   'Re-run with params.overwrite = true.']);
        else
            phyAll{s} = nan(1, size(rasterAll{s}, 1));                   % fill NaN if unused
        end
    end

    % Reconstruct per-stimulus tAxis from stored edges
    tAxis = cell(1, numel(params.stimTypes));
    for s = 1:numel(params.stimTypes)
        tAxis{s} = lockedEdges{s}(1:end-1);                             % left edge per bin (ms)
    end

end

% =========================================================================
% LOAD SPATIAL-TUNING TABLE  (only when sorting by spatial tuning)
% =========================================================================
tuningAll = cell(1, numel(params.stimTypes));                            % one tuning vector per stim

if params.sortBy == "spatialTuning"

    % --- Resolve tuning file path ---
    % Auto-construction tries both stim orderings because the on-disk
    % filename may have been saved with stimTypes in a different order.
    % The data inside is the same (rows are filtered by the stimulus column).
    if strlength(params.tuningFile) == 0
        cand1 = sprintf('%s\\Ex_%d-%d_SpatialTuningIndex_%s_RF_prefDir_allResp.mat', ...
            saveDir, exList(1), exList(end), stimLabel);                 % primary: same stim order
        cand2 = sprintf('%s\\Ex_%d-%d_SpatialTuningIndex_%s_RF_prefDir_allResp.mat', ...
            saveDir, exList(1), exList(end), ...
            strjoin(flip(params.stimTypes), '-'));                        % fallback: reversed order
        if exist(cand1, 'file')
            tuningFile = cand1;                                          % prefer primary
        elseif exist(cand2, 'file')
            tuningFile = cand2;                                          % accept reversed order
        else
            error('Spatial-tuning file not found at:\n  %s\nor\n  %s', cand1, cand2);
        end
    else
        tuningFile = char(params.tuningFile);                            % explicit path overrides
    end

    fprintf('Loading spatial-tuning table from:\n  %s\n', tuningFile);

    % --- Load and validate the tuning table ---
    Ttmp  = load(tuningFile);                                            % load .mat
    tflds = fieldnames(Ttmp);                                            % variable names in file
    isTab = cellfun(@(f) istable(Ttmp.(f)), tflds);                      % find first table variable
    if ~any(isTab)
        error('No table variable found inside %s.', tuningFile);
    end
    tuningTable = Ttmp.(tflds{find(isTab, 1)});                          % grab the tuning table

    % Convert categorical columns to native types so == comparisons work.
    % double(categorical) returns category indices, not values — must go
    % through string first to recover the original numeric values.
    if iscategorical(tuningTable.experimentNum)
        tuningTable.experimentNum = str2double(string(tuningTable.experimentNum));
    end
    if iscategorical(tuningTable.phyID)
        tuningTable.phyID = str2double(string(tuningTable.phyID));
    end
    if iscategorical(tuningTable.stimulus)
        tuningTable.stimulus = string(tuningTable.stimulus);
    end

    % Check required columns exist
    varNames = string(tuningTable.Properties.VariableNames);
    assert(ismember("phyID",         varNames), 'Tuning table missing "phyID" column.');
    assert(ismember("experimentNum", varNames), 'Tuning table missing "experimentNum" column.');
    assert(ismember("stimulus",      varNames), 'Tuning table missing "stimulus" column.');
    assert(ismember(params.tuningIndexCol, varNames), ...
        'Tuning table missing requested column "%s".', params.tuningIndexCol);

    % --- Build per-stim tuning vectors aligned to rasterAll rows ---
    for s = 1:numel(params.stimTypes)
        nNeu         = size(rasterAll{s}, 1);                            % neurons for this stim
        tuningAll{s} = nan(1, nNeu);                                     % default NaN = sorted to end

        if nNeu == 0; continue; end                                      % nothing to do

        % Pre-filter to this stimulus for speed
        stimMask = string(tuningTable.stimulus) == params.stimTypes(s);
        subT     = tuningTable(stimMask, :);                             % subtable for this stim

        for k = 1:nNeu
            % Match by experiment AND Phy cluster ID
            row = subT.experimentNum == expAll{s}(k) & ...
                  subT.phyID         == phyAll{s}(k);
            if any(row)
                tuningAll{s}(k) = subT.(params.tuningIndexCol)(find(row, 1));
            end
        end

        nMissing = sum(isnan(tuningAll{s}));                             % unmatched rows
        if nMissing > 0
            warning('plotRaster:tuningMissing', ...
                '[%s] %d / %d neurons missing from tuning table — sorted to end.', ...
                params.stimTypes(s), nMissing, nNeu);
        end
    end
end

% =========================================================================
% SORT NEURONS
% =========================================================================
for s = 1:numel(params.stimTypes)

    data = rasterAll{s};                                                 % nNeurons x nBins
    if isempty(data); continue; end

    if params.sortBy == "peak"
        % --- Sort by peak-response latency ---
        postStimBins = tAxis{s} >= lockedPreBase;                        % post-onset mask

        % Local smoothed copy for peak detection only (avoids double-smoothing)
        if size(data, 2) > 100
            dataForSort = ConvBurstMatrix( ...
                data, fspecial('gaussian', [1 params.GaussianLength+20], 2), 'same');
        else
            dataForSort = ConvBurstMatrix( ...
                data, fspecial('gaussian', [1 params.GaussianLength], 2), 'same');
        end

        [~, peakBin] = max(dataForSort(:, postStimBins), [], 2);         % peak column per neuron
        [~, sortIdx] = sort(peakBin);                                    % early-peaking first

    elseif params.sortBy == "depth"
        % --- Sort by recording depth ---
        [~, sortIdx] = sort(depthAll{s}, 'ascend');                      % shallowest first

    elseif params.sortBy == "spatialTuning"
        % --- Sort by spatial-tuning index ---
        % MissingPlacement='last' sends NaN (unmatched) neurons to the bottom
        [~, sortIdx] = sort(tuningAll{s}, params.tuningSortOrder, ...
            'MissingPlacement', 'last');

    else
        % --- No reordering ---
        sortIdx = 1:size(data, 1);
    end

    % Apply sort to all parallel vectors
    rasterAll{s} = data(sortIdx, :);                                     % reorder PSTH rows
    depthAll{s}  = depthAll{s}(sortIdx);                                 % reorder depths
    expAll{s}    = expAll{s}(sortIdx);                                   % reorder experiment IDs
    phyAll{s}    = phyAll{s}(sortIdx);                                   % reorder Phy IDs

    if params.sortBy == "spatialTuning"
        tuningAll{s} = tuningAll{s}(sortIdx);                            % keep tuning vector aligned
    end

end

% =========================================================================
% PLOT
% =========================================================================

% Short display labels for each stimulus type
stimLegendMap = containers.Map( ...
    {'linearlyMovingBall', 'rectGrid', 'MovingGrating', 'StaticGrating'}, ...
    {'MB',                 'SB',       'MG',            'SG'});

nStim = numel(params.stimTypes);

% ------------------------------------------------------------------
% Global colour limits — computed once, shared across all tiles
% ------------------------------------------------------------------
allValues = [];
for s = 1:nStim
    if ~isempty(rasterAll{s})
        allValues = [allValues, rasterAll{s}(:)']; %#ok<AGROW>          % pool all values
    end
end

if params.zScore
    cLimPos = prctile(allValues, params.climPrctile);                    % data-driven upper limit
    cLims   = [-params.climNeg, cLimPos];                                % fixed lower, data-driven upper
else
    % Asymmetric percentile clipping: 2nd pctl as floor, climPrctile as ceiling
    cLims = [prctile(allValues, 2), prctile(allValues, params.climPrctile)];
end

% ------------------------------------------------------------------
% Build colormap once
% ------------------------------------------------------------------
if params.zScore && params.colormap ~= "gray"

    % FIX 8: warn if climNeg=0 collapses the negative half
    if params.climNeg == 0
        warning('plotRaster:noNegRange', ...
            'climNeg = 0 with a diverging colormap: negative half is empty. Set climNeg > 0 for a true diverging scale.');
    end

    nColors  = 256;                                                      % total colour entries
    nNeg     = round(nColors * params.climNeg / (params.climNeg + cLimPos)); % entries for negative half
    nPos     = nColors - nNeg;                                           % entries for positive half
    blueHalf = [linspace(0.1,1,nNeg)', linspace(0.2,1,nNeg)', linspace(0.8,1,nNeg)']; % blue -> white
    redHalf  = [linspace(1,0.9,nPos)', linspace(1,0.2,nPos)', linspace(1,0.05,nPos)']; % white -> red
    cmapToUse = [blueHalf; redHalf];                                     % diverging map
else
    cmapToUse = flipud(gray);                                            % white = low, black = high
end

% ------------------------------------------------------------------
% Figure and tiled layout
% ------------------------------------------------------------------
fig = figure;
set(fig, 'Units', 'centimeters', 'Position', [5 5 5*nStim + 2, 10]);    % width scales with tile count

tl    = tiledlayout(fig, 1, nStim, 'TileSpacing', 'compact', 'Padding', 'compact');
axAll = gobjects(1, nStim);                                              % pre-allocate axes handles

for s = 1:nStim

    data    = rasterAll{s};                                              % PSTH matrix for this stim
    stimKey = char(params.stimTypes(s));                                  % char key for Map lookup
    if isKey(stimLegendMap, stimKey)
        shortName = stimLegendMap(stimKey);                              % abbreviated label
    else
        shortName = stimKey;                                             % fallback
    end

    axAll(s) = nexttile(tl);                                             % create tile
    ax = axAll(s);

    if isempty(data)
        title(ax, shortName, 'FontName', 'helvetica', 'FontSize', 8);
        axis(ax, 'off');                                                 % nothing to plot
        continue
    end

    % --- Per-stimulus time axis in seconds ---
    tAxisPlot = tAxis{s} - lockedPreBase;                                % stim onset = 0 ms
    tAxisSec  = tAxisPlot / 1000;                                        % convert to seconds

    imagesc(ax, tAxisSec, 1:size(data,1), data);                         % raster image
    clim(ax, cLims);                                                     % shared colour limits
    colormap(ax, cmapToUse);                                             % shared colormap

    % ------------------------------------------------------------------
    % Depth-bin boundary lines (depth sort only)
    % ------------------------------------------------------------------
    if params.sortBy == "depth" && ~isempty(depthAll{s})

        D2 = load(depthFile);                                            % reload bin edges
        depthBinEdges = D2.depthBinEdges;                                % depth layer edges (um)

        binLabelsDepth = { ...
            sprintf('%.0f-%.0f um', depthBinEdges(1), depthBinEdges(2)), ...
            sprintf('%.0f-%.0f um', depthBinEdges(2), depthBinEdges(3)), ...
            sprintf('%.0f-%.0f um', depthBinEdges(3), depthBinEdges(4))};

        depthCombined = depthAll{s}(~isnan(depthAll{s}));                % exclude NaN
        labelX = tAxisSec(1) + 0.05 * range(tAxisSec);                  % label x pos: 5% from left

        for edge = 2:3                                                   % internal boundaries only
            lastInBin = find(depthCombined <= depthBinEdges(edge), 1, 'last');
            if ~isempty(lastInBin) && lastInBin < size(data,1)
                yline(ax, lastInBin + 0.5, 'k-', 'LineWidth', 1.2);     % horizontal separator
                text(ax, labelX, lastInBin - size(data,1)*0.02, ...
                    binLabelsDepth{edge-1}, ...
                    'Color', 'w', 'FontSize', 6, 'FontName', 'helvetica', ...
                    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
            end
        end
        text(ax, labelX, size(data,1), binLabelsDepth{3}, ...           % deepest bin label
            'Color', 'w', 'FontSize', 6, 'FontName', 'helvetica', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    end

    % --- Stim onset / offset lines (seconds) ---
    xline(ax, 0,                  'k--', 'LineWidth', 1.0);             % onset at t = 0
    xline(ax, stimDurAll(s)/1000, 'k--', 'LineWidth', 1.0);             % offset in seconds

    % --- Axis formatting ---
    xlim(ax, [tAxisSec(1), tAxisSec(end)]);                              % per-stim x range
    ylim(ax, [0.5, size(data,1) + 0.5]);                                 % half-row margin

    ticksSec     = linspace(tAxisSec(1), tAxisSec(end), 5);              % 5 equally spaced ticks
    [~, iz]      = min(abs(ticksSec));                                   % tick nearest to 0
    ticksSec(iz) = 0;                                                    % snap to exactly 0
    xticks(ax, ticksSec);
    xticklabels(ax, arrayfun(@(v) sprintf('%.2g', v), ticksSec, 'UniformOutput', false));

    if s == 1
        ylabel(ax, 'Neuron #', 'FontName', 'helvetica', 'FontSize', 8); % y-label on leftmost tile
    end

    title(ax, sprintf('%s  (n=%d)', shortName, size(data,1)), ...
        'FontName', 'helvetica', 'FontSize', 8);

    ax.FontName = 'helvetica';
    ax.FontSize = 8;
    ax.YDir     = 'normal';                                              % neuron 1 at bottom
    ax.TickDir  = 'out';                                                 % outward ticks
    ax.Box      = 'off';                                                 % no top/right border

end

% ------------------------------------------------------------------
% Shared x-label
% ------------------------------------------------------------------
xlabel(tl, 'Time relative to stimulus onset (s)', ...
    'FontName', 'helvetica', 'FontSize', 8);

% ------------------------------------------------------------------
% Single colorbar on the rightmost tile
% ------------------------------------------------------------------
cb = colorbar(axAll(end));                                               % attach to last axes
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
    'FontName', 'helvetica', 'FontSize', 10);                            % super-title

if params.PaperFig
    vs_first.printFig(fig, sprintf('Raster-%s', stimLabel), PaperFig=params.PaperFig);
end

end