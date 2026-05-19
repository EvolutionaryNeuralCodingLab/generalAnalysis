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
%        spatial-tuning index, preferred category level, or leaves
%        them unsorted.
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
%   NEW  — sortBy = "preferredCategory": group neurons by their preferred
%          category level (e.g. direction, size), show only
%          preferred-level trials in each PSTH row, secondary sort by
%          mean post-onset firing rate within each group.

arguments
    exList double                                                        % vector of experiment IDs to include
    params.stimTypes       (1,:) string  = ["RG","MB"]                    % stimulus abbreviations — one tile each (RG|MB|MBR|SDGs|SDGm|NI|NV|FFF)
    params.binWidth        double        = 10                            % PSTH bin width in ms
    params.smooth          double        = 0                             % Gaussian smoothing SD in ms (0 = off)
    params.statType        string        = "MaxPermuteTest"              % which statistics field to use
    params.speed           string        = "max"                         % ball-speed selector: "max" or other
    params.alpha           double        = 0.05                          % significance threshold
    params.postStim        double        = 0                             % post-stimulus window in ms
    params.preBase         double        = 200                           % pre-stimulus baseline in ms
    params.overwrite       logical       = false                         % if true, recompute even if cache exists
    params.TakeTopPercentTrials double   = 1                           % fraction (0,1] of trials to keep; [] or 0 = keep all
    params.zScore          logical       = true                          % z-score each neuron using its baseline
    params.sortBy          string        = "spatialTuning"               % "peak" | "depth" | "spatialTuning" | "preferredCategory" | "none"
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
    % --- Preferred-category sort parameters ---
    params.splitCategory   (1,:) string  = ""                            % per-stim category: e.g. ["size","direction"]; scalar is broadcast
    params.splitLevels     cell          = {}                            % per-stim levels: e.g. {[1 2],[0 90]}; {} = all levels for every stim
    % --- Per-category statistics parameters (used when splitCategory is active) ---
    params.nBootCategory   double        = 10000                         % bootstrap iterations for StatisticsPerNeuronPerCategory
    params.overwriteCatStats logical     = false                         % force recomputation of per-category statistics
    params.catBaseRespWindow double      = 100                           % base response window (ms) for per-category statistics
    params.catApplyFDR     logical       = false                         % apply FDR inside StatisticsPerNeuronPerCategory
    params.useGeneralFilter logical      = false                         % true = use general StatisticsPerNeuron p-values even in category mode
end

% -------------------------------------------------------------------------
% Sanity check: baseline must be an integer multiple of bin width
% -------------------------------------------------------------------------
assert(mod(params.preBase, params.binWidth) == 0, ...                    % prevents misaligned baseline bin edges
    'preBase (%g ms) must be a multiple of binWidth (%g ms).', ...
    params.preBase, params.binWidth);

% -------------------------------------------------------------------------
% Validate and broadcast preferredCategory parameters
% -------------------------------------------------------------------------
if params.sortBy == "preferredCategory"

    nStimTypes = numel(params.stimTypes);                                % number of stimulus types

    % --- Broadcast splitCategory ---
    % A single string is applied to all stim types; a vector must match nStim.
    % An empty string "" for a given slot means that stim uses all-trial mode.
    if isscalar(params.splitCategory)
        params.splitCategory = repmat(params.splitCategory, 1, nStimTypes); % broadcast scalar to all stim types
    end
    assert(numel(params.splitCategory) == nStimTypes, ...
        'splitCategory must be scalar or have one entry per stimType (%d).', nStimTypes);

    % --- Broadcast splitLevels ---
    % An empty cell {} means "all levels" for every stim.
    % A cell with one element is broadcast to all stim types.
    % A cell with nStim elements maps one-to-one.
    if isempty(params.splitLevels)
        params.splitLevels = repmat({[]}, 1, nStimTypes);                % all levels for every stim
    elseif numel(params.splitLevels) == 1
        params.splitLevels = repmat(params.splitLevels, 1, nStimTypes);  % broadcast single cell to all stim
    end
    assert(numel(params.splitLevels) == nStimTypes, ...
        'splitLevels must be empty, scalar cell, or have one entry per stimType (%d).', nStimTypes);

    % --- Ensure at least one stim has a non-empty category ---
    assert(any(strlength(params.splitCategory) > 0), ...
        'sortBy="preferredCategory" requires at least one non-empty splitCategory entry.');
end

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
    cd(p)                                                                % change to parent dir
    mkdir Combined_lizard_analysis                                       % create sub-directory
end
saveDir = [p '\Combined_lizard_analysis'];                               % output directory

stimLabel  = strjoin(params.stimTypes, '-');                              % e.g. "RG-MB"

% --- Include splitCategory in cache filename to avoid collisions ---
% Join all per-stim categories into a compact suffix, e.g. "_cat-size-direction"
if params.sortBy == "preferredCategory"
    nonEmpty = params.splitCategory(strlength(params.splitCategory) > 0); % only non-empty categories
    catSuffix = sprintf('_cat-%s', strjoin(nonEmpty, '-'));               % e.g. "_cat-size-direction"
else
    catSuffix = '';                                                       % no suffix for other sort modes
end
nameOfFile = sprintf('\\Ex_%d-%d_Raster_%s%s.mat', ...
    exList(1), exList(end), stimLabel, catSuffix);                       % cache filename

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
    rasterAll    = cell(1, nStim);                                       % nNeurons x nBins PSTH per stim
    depthAll     = cell(1, nStim);                                       % recording depth (um) per neuron
    expAll       = cell(1, nStim);                                       % experiment ID per neuron row
    phyAll       = cell(1, nStim);                                       % Phy cluster ID per neuron row
    prefLevelAll = cell(1, nStim);                                       % preferred category level per neuron (NEW)

    for s = 1:nStim
        rasterAll{s}    = [];                                            % initialise empty PSTH matrix
        depthAll{s}     = [];                                            % initialise empty depth vector
        expAll{s}       = [];                                            % initialise empty exp-ID vector
        phyAll{s}       = [];                                            % initialise empty Phy-ID vector
        prefLevelAll{s} = [];                                            % initialise empty preferred-level vector
    end

    % Counter for neurons dropped due to zero-SD baseline
    nDroppedZeroSD = zeros(1, nStim);                                    % per-stim counter

    % Counter for experiments skipped due to insufficient category levels
    nSkippedCat = zeros(1, nStim);                                       % per-stim counter (NEW)

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
                stimKey = params.stimTypes(s);                           % current stim abbreviation

                % --- Build analysis object from abbreviation ---
                gt = detectGratingType(stimKey);                         % '' for non-grating, 'moving'/'static' for SDG
                switch stimKey
                    case "RG";   objTmp = rectGridAnalysis(NPtmp);
                    case "MB";   objTmp = linearlyMovingBallAnalysis(NPtmp);
                    case "MBR";  objTmp = linearlyMovingBarAnalysis(NPtmp);
                    case {"SDGs"}
                        objTmp = StaticDriftingGratingAnalysis(NPtmp);
                    case {"SDGm"}
                        objTmp = StaticDriftingGratingAnalysis(NPtmp);
                    case "NI";   objTmp = imageAnalysis(NPtmp);
                    case "NV";   objTmp = movieAnalysis(NPtmp);
                    case "FFF";  objTmp = fullFieldFlashAnalysis(NPtmp);
                    otherwise;   error('Unknown stimType abbreviation: %s', stimKey);
                end
                NRtmp = objTmp.ResponseWindow;                           % response-window struct

                % Resolve fieldName (same logic as main loop)
                fn = resolveFieldName(stimKey, params.speed);            % sub-field key for this stim type

                try
                    dur = NRtmp.(fn).stimDur;                            % sub-field duration
                catch
                    dur = NRtmp.stimDur;                                  % flat struct fallback
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
            continue                                                     % skip on load failure
        end

        for s = 1:nStim

            stimType = params.stimTypes(s);                              % current stim abbreviation

            % --- Build stimulus-specific analysis object ---
            try
                gt = detectGratingType(stimType);                        % '' for non-grating, 'moving'/'static' for SDG
                switch stimType
                    case "RG"
                        obj = rectGridAnalysis(NP);                      % receptive-field grid stimulus
                    case "MB"
                        obj = linearlyMovingBallAnalysis(NP);            % moving ball stimulus
                    case "MBR"
                        obj = linearlyMovingBarAnalysis(NP);             % moving bar stimulus
                    case {"SDGs","SDGm"}
                        obj = StaticDriftingGratingAnalysis(NP); % grating stimulus
                    case "NI"
                        obj = imageAnalysis(NP);                         % natural image stimulus
                    case "NV"
                        obj = movieAnalysis(NP);                         % natural movie stimulus
                    case "FFF"
                        obj = fullFieldFlashAnalysis(NP);                % full-field flash stimulus
                    otherwise
                        error('Unknown stimType abbreviation: %s', stimType);
                end
                NeuronResp = obj.ResponseWindow;                             % response-window struct
            catch ME
                warning('Could not build %s for exp %d: %s', stimType, ex, ME.message);
                continue                                                 % skip this stim/exp
            end

          

            % --- Select statistics struct ---
            if params.statType == "BootstrapPerNeuron"
                Stats = obj.BootstrapPerNeuron;                          % bootstrap p-values
            else
                Stats = obj.StatisticsPerNeuron;                         % default: permutation test
            end

            % --- Resolve sub-field name and stim-onset offset ---
            fieldName = resolveFieldName(stimType, params.speed);        % sub-field key for this stim type
            startStim = 0;                                               % ms offset for stim onset
            if stimType == "SDGm"
                startStim = obj.VST.static_time * 1000;                  % moving phase onset (s -> ms)
            end

            % --- Convert Phy sorting to tIc format ---
            p_sort = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
            label  = string(p_sort.label');                              % quality label per unit
            goodU  = p_sort.ic(:, label == 'good');                      % keep only curated 'good' units

            % --- Extract Phy cluster IDs for good units ---
            goodPhyIDs = p_sort.phy_ID(label == 'good');                 % Phy cluster IDs matching goodU columns

            % --- General response p-values (StatisticsPerNeuron) ---
            % These are always extracted: used directly in standard mode,
            % or as a fallback when useGeneralFilter = true in category mode.
            try
                pvals = Stats.(fieldName).pvalsResponse;                 % stim-specific p-values
            catch
                pvals = Stats.pvalsResponse;                             % flat struct fallback
            end

            % --- Stimulus onset times and condition matrix ---
            try
                C = NeuronResp.(fieldName).C;                            % condition matrix (sub-field)
            catch
                C = NeuronResp.C;                                        % condition matrix (flat)
            end
            directimesSorted = C(:, 1)' + startStim;                     % onset times in ms

            % =============================================================
            % CATEGORY DETECTION (for preferredCategory mode)
            % =============================================================
            % isCatMode is per-stimulus: true only if sortBy is
            % preferredCategory AND this stim slot has a non-empty category.
            isCatMode = params.sortBy == "preferredCategory" && ...
                        strlength(params.splitCategory(s)) > 0;          % per-stim flag

            if isCatMode

                thisCatName   = params.splitCategory(s);                 % category for THIS stim (e.g. "direction")
                thisCatLevels = params.splitLevels{s};                   % requested levels for THIS stim ([] = all)

                % --- Extract column names from ResponseWindow ---
                try
                    allColNames = NeuronResp.(fieldName).colNames{1};    % sub-field column names
                catch
                    try
                        allColNames = NeuronResp.colNames{1};            % flat struct column names
                    catch
                        warning('[%s] exp %d: colNames not found — skipping.', stimType, ex);
                        nSkippedCat(s) = nSkippedCat(s) + 1;            % count skipped experiment
                        continue                                         % skip this stim/exp
                    end
                end

                % Column names: first 4 are metadata (onset, etc.),
                % entries 5+ are stimulus-parameter names.
                catColNames = string(allColNames(5:end));                 % stimulus-parameter names only

                % --- Find column index for the requested category ---
                catIdx = find(strcmpi(catColNames, thisCatName), 1);     % case-insensitive match

                if isempty(catIdx)
                    warning('[%s] exp %d: category "%s" not found in colNames [%s] — skipping.', ...
                        stimType, ex, thisCatName, strjoin(catColNames, ', '));
                    nSkippedCat(s) = nSkippedCat(s) + 1;                 % count skipped experiment
                    continue                                             % skip: category absent
                end

                % BUG 10 FIX: C column 1 = onset times; columns 2+ = stimulus
                % parameters matching catColNames.  Since catColNames is
                % already colNames{1}(5:end), catIdx is a 1-based index
                % into the parameter names.  Offset by 1 for the time column.
                catColInC = catIdx + 1;                                  % +1 for onset-time column 1

                % --- Extract category values per trial ---
                trialCatValues = C(:, catColInC);                        % category value for each trial

                % --- Determine available levels ---
                availableLevels = unique(trialCatValues);                % unique levels in this experiment

                % --- Filter to requested levels if specified ---
                if ~isempty(thisCatLevels)
                    availableLevels = intersect(availableLevels, thisCatLevels(:)); % keep only requested
                end

                % --- Skip if fewer than 2 levels remain ---
                if numel(availableLevels) < 2
                    warning('[%s] exp %d: only %d level(s) of "%s" after filtering — skipping.', ...
                        stimType, ex, numel(availableLevels), thisCatName);
                    nSkippedCat(s) = nSkippedCat(s) + 1;                 % count skipped experiment
                    continue                                             % skip: can't determine preference
                end

                fprintf('  [%s] exp %d: %d levels of "%s": [%s]\n', ...
                    stimType, ex, numel(availableLevels), thisCatName, ...
                    strjoin(string(availableLevels'), ', '));

                % ==========================================================
                % PER-CATEGORY STATISTICS (StatisticsPerNeuronPerCategory)
                % ==========================================================
                % When useGeneralFilter is false (default), we use per-level
                % p-values from StatisticsPerNeuronPerCategory to filter
                % neurons.  A neuron passes if it is significant for at
                % least one level (OR across levels).  This mirrors the
                % AllExpAnalysis Mode 2 convention and is methodologically
                % correct: neurons should be responsive to the specific
                % category being split, not just to the stimulus in general.
                if ~params.useGeneralFilter
                    gt = detectGratingType(stimType);                    % '' for non-grating, 'moving'/'static' for SDG

                    % Build name-value args for StatisticsPerNeuronPerCategory
                    catStatsArgs = { ...
                        'compareCategory', char(thisCatName), ...        % category to decompose
                        'nBoot',           params.nBootCategory, ...     % bootstrap iterations
                        'overwrite',       params.overwriteCatStats, ... % force recomputation flag
                        'BaseRespWindow',  params.catBaseRespWindow, ... % response window (ms)
                        'applyFDR',        params.catApplyFDR};          % FDR inside the stats function
                    if ~isempty(gt)
                        catStatsArgs = [catStatsArgs, {'GratingType', gt}]; % pass GratingType only for SDG stim
                    end

                    % Run per-category statistics (results are cached by the object)
                    catStats = obj.StatisticsPerNeuronPerCategory(catStatsArgs{:});

                    % Build OR mask: neuron passes if significant for ANY available level
                    nUnits = numel(pvals);                               % total good units (same as general stats)
                    orMask = false(nUnits, 1);                           % initialise all-false

                    for li = 1:numel(availableLevels)
                        fName = levelToFieldName(char(thisCatName), availableLevels(li)); % e.g. 'direction_0'
                        if isfield(catStats, fName)
                            levelP = catStats.(fName).pvalsResponse(:);  % per-level p-values
                            orMask = orMask | (levelP < params.alpha);   % OR into mask
                        else
                            warning('[%s] exp %d: catStats missing field "%s" — skipping level.', ...
                                stimType, ex, fName);
                        end
                    end

                    fprintf('  [%s] exp %d: %d / %d neurons pass per-category OR filter.\n', ...
                        stimType, ex, sum(orMask), nUnits);
                end % ~useGeneralFilter

            end % isCatMode detection block

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
                lockedPreBase = preBase;                                  % shared across all stimuli
            end

            % Lock bin edges per stim on first encounter
            if isempty(lockedEdges{s})
                lockedEdges{s} = 0 : params.binWidth : windowTotal;      % bin edges from 0 to windowTotal
                lockedNBins(s) = numel(lockedEdges{s}) - 1;              % number of bins
                tAxis{s}       = lockedEdges{s}(1:end-1);               % left edge per bin (ms)
                stimDurAll(s)  = rawStimDur_ms;                          % for xline in plot
                fprintf('  [%s] Locked window: preBase=%d ms, stimDur=%.0f ms, nBins=%d\n', ...
                    stimType, lockedPreBase, rawStimDur_ms, lockedNBins(s));
            end

            % --- Find responsive neurons ---
            % In category mode with per-level filtering (default), use
            % the OR mask from StatisticsPerNeuronPerCategory.
            % In all other cases, use general StatisticsPerNeuron p-values.
            if isCatMode && ~params.useGeneralFilter
                eNeurons = find(orMask);                                 % per-category OR filter
            else
                eNeurons = find(pvals < params.alpha);                   % general filter
            end

            if isempty(eNeurons)
                fprintf('  [%s] No responsive neurons in exp %d.\n', stimType, ex);
                continue
            end

            fprintf('  [%s] %d responsive neuron(s) in exp %d.\n', ...
                stimType, numel(eNeurons), ex);

            % ----------------------------------------------------------
            % Per-neuron PSTH
            % ----------------------------------------------------------
            nAppended = 0;                                               % count neurons actually added to raster
            for ni = 1:numel(eNeurons)

                u = eNeurons(ni);                                        % index into the good-unit list

                % ==============================================================
                % ALL-TRIALS BASELINE (category mode only)
                % ==============================================================
                % In category mode, compute baseline stats from ALL trials
                % before selecting the preferred level.  The pre-stimulus
                % period is stimulus-independent (the animal cannot predict
                % the upcoming category), so pooling all trials gives the
                % most stable baseline estimate.  In standard mode bMean/bStd
                % are left empty and computed later from the (single) PSTH.
                bMean = [];                                              % sentinel: compute later in standard mode
                bStd  = [];

                if isCatMode && params.zScore

                    % Build all-trials spike matrix (just for baseline)
                    MRall = BuildBurstMatrix( ...
                        goodU(:, u), ...                                 % spike identity for this unit
                        round(p_sort.t), ...                             % rounded sample timestamps
                        round(directimesSorted - lockedPreBase), ...     % ALL trial starts
                        round(windowTotal));                             % window length (ms)
                    MRall = squeeze(MRall);                              % remove singleton dims

                    % Build all-trials PSTH (only baseline portion matters)
                    nTrialsAll    = size(MRall, 1);                      % total number of trials
                    spikeTimesAll = repmat((1:size(MRall,2)), nTrialsAll, 1); % column indices
                    spikeTimesAll = spikeTimesAll(logical(MRall));       % keep spike positions only
                    countsAll     = histcounts(spikeTimesAll, lockedEdges{s}); % spike count per bin
                    allTrialPSTH  = (countsAll / (params.binWidth * nTrialsAll)) * 1000; % spk/s

                    % Extract baseline stats from the all-trials PSTH
                    baselineBins = tAxis{s} < lockedPreBase;             % pre-stimulus bin mask
                    bMean = mean(allTrialPSTH(baselineBins));            % all-trial baseline mean
                    bStd  = std(allTrialPSTH(baselineBins));             % all-trial baseline SD

                    % Early exit if baseline SD is zero — z-score is
                    % undefined regardless of which level we pick.
                    if bStd == 0
                        nDroppedZeroSD(s) = nDroppedZeroSD(s) + 1;
                        continue                                         % skip neuron
                    end
                end

                % ==========================================================
                % PREFERRED-CATEGORY MODE: determine preferred level and
                % build PSTH from only that level's trials
                % ==========================================================
                if isCatMode

                    % --- Compute mean post-onset rate per level ---
                    % Use ALL trials (before TakeTopPercentTrials) for a
                    % stable preference estimate.
                    nLevels      = numel(availableLevels);               % number of category levels
                    meanRatePerLevel = nan(1, nLevels);                  % pre-allocate rate per level

                    for li = 1:nLevels
                        levelMask = trialCatValues == availableLevels(li); % logical mask: trials of this level
                        levelOnsets = directimesSorted(levelMask);       % onset times for this level's trials

                        if isempty(levelOnsets)
                            continue                                     % no trials for this level
                        end

                        % Build binary spike matrix for this level's trials
                        MRtemp = BuildBurstMatrix( ...
                            goodU(:, u), ...                             % spike identity for this unit
                            round(p_sort.t), ...                         % rounded sample timestamps
                            round(levelOnsets - lockedPreBase), ...      % trial start = onset - baseline
                            round(windowTotal));                         % window length (ms)
                        MRtemp = squeeze(MRtemp);                        % remove singleton dims

                        % Handle single-trial case: ensure matrix is (1 x T)
                        if isvector(MRtemp) && numel(levelOnsets) == 1
                            MRtemp = MRtemp(:)';                         % force row vector
                        end

                        % Mean firing rate in the post-onset window (spk/ms -> spk/s)
                        postCols = (lockedPreBase + 1) : size(MRtemp, 2); % columns after stim onset
                        meanRatePerLevel(li) = mean(MRtemp(:, postCols), 'all') * 1000; % convert to spk/s
                    end

                    % --- Determine preferred level by argmax ---
                    [bestRate, bestIdx] = max(meanRatePerLevel);         % highest mean rate across levels

                    if isnan(bestRate)
                        continue                                         % all levels empty — skip neuron
                    end

                    prefLevel = availableLevels(bestIdx);                 % the preferred category level

                    % --- Build PSTH from preferred-level trials only ---
                    prefMask   = trialCatValues == prefLevel;            % logical mask for preferred level
                    prefOnsets = directimesSorted(prefMask);             % onset times of preferred trials

                    MRhist = BuildBurstMatrix( ...
                        goodU(:, u), ...                                 % spike identity for this unit
                        round(p_sort.t), ...                             % rounded sample timestamps
                        round(prefOnsets - lockedPreBase), ...           % trial start = onset - baseline
                        round(windowTotal));                             % window length (ms)
                    MRhist = squeeze(MRhist);                            % remove singleton dims

                    % Handle single-trial case
                    if isvector(MRhist) && numel(prefOnsets) == 1
                        MRhist = MRhist(:)';                             % force row vector
                    end

                else
                    % ==================================================
                    % STANDARD MODE: build PSTH from all trials
                    % ==================================================
                    prefLevel = NaN;                                      % not applicable in standard mode

                    % Binary spike matrix: trials x time at 1 ms resolution
                    MRhist = BuildBurstMatrix( ...
                        goodU(:, u), ...                                 % spike identity for this unit
                        round(p_sort.t), ...                             % rounded sample timestamps
                        round(directimesSorted - lockedPreBase), ...     % trial start = onset - baseline
                        round(windowTotal));                             % window length (ms, integer)
                    MRhist = squeeze(MRhist);                            % remove singleton dims
                end

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

                    if isempty(bMean)
                        % Standard mode: compute baseline from this PSTH
                        bMean = mean(neuronPSTH(baselineBins));          % baseline mean
                        bStd  = std(neuronPSTH(baselineBins));           % baseline SD
                    end
                    % Category mode: bMean/bStd already set from all-trials
                    % PSTH above (more stable, stimulus-independent estimate)

                    if bStd > 0
                        neuronPSTH = (neuronPSTH - bMean) / bStd;       % z-score normalisation
                    else
                        % FIX 5: count and log rather than silently skip
                        nDroppedZeroSD(s) = nDroppedZeroSD(s) + 1;
                        continue                                         % skip: undefined z-score
                    end
                end

                % --- Smooth PSTH if requested ---
                if params.smooth > 0
                    smoothBins = round(params.smooth / params.binWidth); % smoothing SD in bin units
                    neuronPSTH = smoothdata(neuronPSTH, 'gaussian', smoothBins); % smooth in-place
                end

                % --- Append neuron row to accumulators ---
                rasterAll{s}    = [rasterAll{s}; neuronPSTH];           % PSTH row
                phyAll{s}(end+1)    = goodPhyIDs(u);                    % Phy cluster ID for this unit
                expAll{s}(end+1)    = ex;                                % source experiment
                prefLevelAll{s}(end+1) = prefLevel;                     % preferred level (NaN if not catMode)
                nAppended = nAppended + 1;                               % count actually-appended neurons

                % Store recording depth
                if params.sortBy == "depth"
                    depthRow = depthTable.Experiment == ex & depthTable.Unit == u;
                    if any(depthRow)
                        depthAll{s}(end+1) = depthTable.Depth_um(depthRow); % matched depth
                    else
                        depthAll{s}(end+1) = NaN;                       % not found
                    end
                else
                    depthAll{s}(end+1) = NaN;                            % unused; keeps vector aligned
                end

            end % neuron loop

            % Report if any responsive neurons were dropped during PSTH building
            nDropped = numel(eNeurons) - nAppended;                      % neurons lost to zero-SD or empty levels
            if nDropped > 0
                fprintf('  [%s] exp %d: %d / %d responsive neurons dropped (zero-SD or empty preferred level).\n', ...
                    stimType, ex, nDropped, numel(eNeurons));
            end

        end % stimulus loop
    end % experiment loop

    % FIX 5: report zero-SD dropped neurons
    for s = 1:nStim
        if nDroppedZeroSD(s) > 0
            fprintf('  [%s] Dropped %d neuron(s) with zero baseline SD.\n', ...
                params.stimTypes(s), nDroppedZeroSD(s));
        end
    end

    % Report experiments skipped due to insufficient category levels
    if params.sortBy == "preferredCategory"
        for s = 1:nStim
            if nSkippedCat(s) > 0
                fprintf('  [%s] Skipped %d experiment(s) with <2 levels of "%s".\n', ...
                    params.stimTypes(s), nSkippedCat(s), params.splitCategory(s));
            end
        end
    end

    % FIX 9: verify accumulator alignment after all experiments
    for s = 1:nStim
        nRows = size(rasterAll{s}, 1);                                   % number of neuron rows
        assert(numel(expAll{s})       == nRows, 'expAll{%d} length mismatch.', s);
        assert(numel(phyAll{s})       == nRows, 'phyAll{%d} length mismatch.', s);
        assert(numel(depthAll{s})     == nRows, 'depthAll{%d} length mismatch.', s);
        assert(numel(prefLevelAll{s}) == nRows, 'prefLevelAll{%d} length mismatch.', s);
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
        S.(sprintf('%s_raster',    stimField)) = rasterAll{s};           % PSTH matrix
        S.(sprintf('%s_depth',     stimField)) = depthAll{s};            % depth vector
        S.(sprintf('%s_exp',       stimField)) = expAll{s};              % experiment ID vector
        S.(sprintf('%s_phy',       stimField)) = phyAll{s};              % Phy cluster ID vector
        S.(sprintf('%s_prefLevel', stimField)) = prefLevelAll{s};        % preferred level vector (NEW)
    end

    save([saveDir nameOfFile], '-struct', 'S');                          % write cache to disk
    fprintf('\nSaved raster data to:\n  %s\n', [saveDir nameOfFile]);

else
    % ------------------------------------------------------------------
    % Reload cached data
    % ------------------------------------------------------------------
    lockedEdges   = S.lockedEdges;                                       % restore bin edges
    lockedPreBase = S.lockedPreBase;                                     % restore baseline
    stimDurAll    = S.stimDurAll;                                        % restore stim durations

    rasterAll    = cell(1, numel(params.stimTypes));                     % pre-allocate
    depthAll     = cell(1, numel(params.stimTypes));
    expAll       = cell(1, numel(params.stimTypes));
    phyAll       = cell(1, numel(params.stimTypes));
    prefLevelAll = cell(1, numel(params.stimTypes));                     % NEW

    for s = 1:numel(params.stimTypes)
        stimField    = matlab.lang.makeValidName(params.stimTypes(s));   % valid struct field name
        rasterAll{s} = S.(sprintf('%s_raster', stimField));              % restore PSTH matrix
        depthAll{s}  = S.(sprintf('%s_depth',  stimField));              % restore depth vector
        expAll{s}    = S.(sprintf('%s_exp',    stimField));              % restore experiment IDs

        % phyAll may be absent in old caches — require recompute if needed
        phyField = sprintf('%s_phy', stimField);
        if isfield(S, phyField)
            phyAll{s} = S.(phyField);                                    % restore Phy IDs
        elseif params.sortBy == "spatialTuning" || params.sortBy == "preferredCategory"
            error(['Cache file lacks Phy IDs (old format). ' ...
                   'Re-run with params.overwrite = true.']);
        else
            phyAll{s} = nan(1, size(rasterAll{s}, 1));                   % fill NaN if unused
        end

        % prefLevelAll may be absent in old caches
        plField = sprintf('%s_prefLevel', stimField);
        if isfield(S, plField)
            prefLevelAll{s} = S.(plField);                               % restore preferred levels
        elseif params.sortBy == "preferredCategory"
            error(['Cache file lacks prefLevel data (old format). ' ...
                   'Re-run with params.overwrite = true.']);
        else
            prefLevelAll{s} = nan(1, size(rasterAll{s}, 1));             % fill NaN if unused
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

        % Pre-filter to this stimulus for speed.
        % Use abbrevToLegacyNames so matching works whether the tuning
        % table was built with abbreviations or legacy full names.
        legacyNames = abbrevToLegacyNames(params.stimTypes(s));          % e.g. ["MB","linearlyMovingBall"]
        stimMask = ismember(string(tuningTable.stimulus), legacyNames);  % match any known name variant
        subT     = tuningTable(stimMask, :);                             % subtable for this stim

        for k = 1:nNeu
            % Match by experiment AND Phy cluster ID
            row = subT.experimentNum == expAll{s}(k) & ...
                  subT.phyID         == phyAll{s}(k);
            if any(row)
                tuningAll{s}(k) = subT.(params.tuningIndexCol)(find(row, 1)); % tuning index value
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
% Also store level-group boundaries for plotting (preferredCategory mode)
levelGroupBounds = cell(1, numel(params.stimTypes));                     % {s}: struct with edges and labels

for s = 1:numel(params.stimTypes)

    data = rasterAll{s};                                                 % nNeurons x nBins
    if isempty(data); continue; end

    if params.sortBy == "peak"
        % --- Sort by peak-response latency ---
        postStimBins = tAxis{s} >= lockedPreBase;                        % post-onset mask

        % Local smoothed copy for peak detection only (avoids double-smoothing)
        if size(data, 2) > 100
            dataForSort = ConvBurstMatrix( ...
                data, fspecial('gaussian', [1 params.GaussianLength+20], 2), 'same'); % wider kernel for long windows
        else
            dataForSort = ConvBurstMatrix( ...
                data, fspecial('gaussian', [1 params.GaussianLength], 2), 'same'); % narrow kernel for short windows
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
            'MissingPlacement', 'last');                                  % most-tuned first (default)

    elseif params.sortBy == "preferredCategory"
        % --- Sort by preferred category level, then by response strength ---
        % If this stim slot has no category assigned (splitCategory(s) == ""),
        % fall through to unsorted (identity permutation).
        if strlength(params.splitCategory(s)) == 0
            sortIdx = 1:size(data, 1);                                   % no category for this stim: no reorder

        else
            % Primary: group neurons by preferred level (ascending level value)
            % Secondary: within each group, sort by mean post-onset response
            %            (descending, so strongest responder at top of group)

            levels   = prefLevelAll{s};                                  % preferred level per neuron
            nNeurons = size(data, 1);                                    % total neurons for this stim

            % Compute mean post-onset response for secondary sort
            postStimBins = tAxis{s} >= lockedPreBase;                    % post-onset mask (same as in peak sort)
            meanPostResp = mean(data(:, postStimBins), 2)';              % 1 x nNeurons: mean post-onset value

            % Get unique levels present (excluding NaN, which shouldn't occur
            % here but is handled defensively)
            uniqueLevels = unique(levels(~isnan(levels)));                % sorted ascending by default

            % Build sort index: iterate levels in order, within each group
            % sort by response strength
            sortIdx   = zeros(1, nNeurons);                              % pre-allocate
            cursor    = 0;                                               % running position counter
            groupInfo = struct('edges', [], 'labels', {{}});             % for plot annotations

            for li = 1:numel(uniqueLevels)
                groupMask    = levels == uniqueLevels(li);               % neurons preferring this level
                groupIndices = find(groupMask);                          % their row indices
                groupResp    = meanPostResp(groupMask);                  % their mean responses

                [~, withinOrder] = sort(groupResp, 'descend');           % strongest first within group
                nGroup = numel(groupIndices);                            % neurons in this group

                sortIdx(cursor+1 : cursor+nGroup) = groupIndices(withinOrder); % fill sorted indices

                % Record group boundary for plotting
                groupInfo.edges(end+1)  = cursor + nGroup;               % row index of last neuron in group
                groupInfo.labels{end+1} = sprintf('%s=%g', ...
                    params.splitCategory(s), uniqueLevels(li));          % label: e.g. "direction=0"

                cursor = cursor + nGroup;                                % advance cursor
            end

            % Handle any NaN-level neurons (shouldn't happen, but defensive)
            nanMask = isnan(levels);
            if any(nanMask)
                nanIndices = find(nanMask);                              % indices of NaN neurons
                nNan       = numel(nanIndices);
                sortIdx(cursor+1 : cursor+nNan) = nanIndices;            % append at end
                groupInfo.edges(end+1)  = cursor + nNan;
                groupInfo.labels{end+1} = 'unclassified';
                cursor = cursor + nNan;
            end

            levelGroupBounds{s} = groupInfo;                             % store for plotting
        end

    else
        % --- No reordering ---
        sortIdx = 1:size(data, 1);                                       % identity permutation
    end

    % Apply sort to all parallel vectors
    rasterAll{s}    = data(sortIdx, :);                                  % reorder PSTH rows
    depthAll{s}     = depthAll{s}(sortIdx);                              % reorder depths
    expAll{s}       = expAll{s}(sortIdx);                                % reorder experiment IDs
    phyAll{s}       = phyAll{s}(sortIdx);                                % reorder Phy IDs
    prefLevelAll{s} = prefLevelAll{s}(sortIdx);                          % reorder preferred levels

    if params.sortBy == "spatialTuning"
        tuningAll{s} = tuningAll{s}(sortIdx);                            % keep tuning vector aligned
    end

end

% =========================================================================
% PLOT
% =========================================================================

% Short display labels for each stimulus type (abbreviations are already
% concise, but the map allows custom labels if desired).
stimLegendMap = containers.Map( ...
    {'RG',  'MB',  'MBR', 'SDGs', 'SDGm', 'NI',  'NV',  'FFF'}, ...
    {'RG',  'MB',  'MBR', 'SDGs', 'SDGm', 'NI',  'NV',  'FFF'});

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

    % ------------------------------------------------------------------
    % Preferred-category group boundary lines and labels (NEW)
    % ------------------------------------------------------------------
    if params.sortBy == "preferredCategory" && ~isempty(levelGroupBounds{s})

        gInfo   = levelGroupBounds{s};
        nGroups = numel(gInfo.edges);

        ytPositions = zeros(1, nGroups);
        ytLabels    = cell(1, nGroups);

        prevEdge = 0;
        for gi = 1:nGroups
            edgeRow  = gInfo.edges(gi);
            nInGroup = edgeRow - prevEdge;

            ytPositions(gi) = (prevEdge + edgeRow) / 2;
            % Level value + count only, e.g. "129 (n=15)"
            rawLabel = gInfo.labels{gi};                                 % e.g. "size=129"
            levelVal = extractAfter(rawLabel, '=');                       % e.g. "129"
            ytLabels{gi} = sprintf('%s (n=%d)', levelVal, nInGroup);

            if gi < nGroups
                yline(ax, edgeRow + 0.5, 'k-', 'LineWidth', 1.5);
            end

            prevEdge = edgeRow;
        end

        set(ax, 'YTick', ytPositions, 'YTickLabel', ytLabels, ...
            'TickLength', [0 0]);
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
    xticks(ax, ticksSec);                                                % apply x tick positions
    xticklabels(ax, arrayfun(@(v) sprintf('%.2g', v), ticksSec, 'UniformOutput', false)); % formatted labels

    if s == 1
        if params.sortBy == "preferredCategory" && strlength(params.splitCategory(s)) > 0
            ylabel(ax, params.splitCategory(s), ...
                'FontName', 'helvetica', 'FontSize', 8);                 % category name as ylabel
        else
            ylabel(ax, 'Neuron #', 'FontName', 'helvetica', 'FontSize', 8);
        end
    end

    title(ax, sprintf('%s  (n=%d)', shortName, size(data,1)), ...
        'FontName', 'helvetica', 'FontSize', 8);                        % tile title with neuron count

    ax.FontName = 'helvetica';                                           % consistent font
    ax.FontSize = 8;                                                     % readable size
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
    cb.Label.String = 'Z-score';                                         % label for z-scored data
else
    cb.Label.String = 'Firing rate (spk/s)';                             % label for raw rates
end
cb.Label.FontName = 'helvetica';                                         % colorbar label font
cb.Label.FontSize = 8;
cb.FontName       = 'helvetica';                                         % tick font
cb.FontSize       = 8;
set(fig, 'Units', 'centimeters', 'Position', [20 20 9 12]);             % final figure position/size

sgtitle(sprintf('N = %d experiments', numel(exList)), ...
    'FontName', 'helvetica', 'FontSize', 10);                           % super-title

if params.PaperFig
    vs_first.printFig(fig, sprintf('Raster-%s-%s', stimLabel,params.splitCategory), PaperFig=params.PaperFig); % export for publication
end

end

% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function gt = detectGratingType(stimKey)
% detectGratingType  Return the GratingType parameter from a stimulus
%   abbreviation.  'moving' for SDGm, 'static' for SDGs, '' otherwise.
    switch stimKey
        case "SDGm"; gt = 'moving';                                      % moving (drifting) grating
        case "SDGs"; gt = 'static';                                      % static grating
        otherwise;   gt = '';                                            % not a grating stimulus
    end
end

function fn = resolveFieldName(stimKey, speedParam)
% resolveFieldName  Map a stimulus abbreviation + speed selector to the
%   sub-field name used in ResponseWindow / StatisticsPerNeuron structs.
%
%   stimKey    — stimulus abbreviation string (e.g. "MB", "SDGs")
%   speedParam — "max" or other speed selector
%   fn         — sub-field key (e.g. 'Speed1', 'Static', '')
    switch stimKey
        case {"MB","MBR"}
            if speedParam == "max"
                fn = 'Speed1';                                           % fastest speed condition
            else
                fn = 'Speed2';                                           % slower speed condition
            end
        case "SDGs"
            fn = 'Static';                                               % static grating sub-field
        case "SDGm"
            fn = 'Moving';                                               % moving grating sub-field
        otherwise
            fn = '';                                                      % flat struct (RG, NI, NV, FFF)
    end
end

function names = abbrevToLegacyNames(stimKey)
% abbrevToLegacyNames  Return the set of stimulus name strings that might
%   appear in external tables (e.g. tuning tables) for a given abbreviation.
%   Includes both the abbreviation itself and any legacy full names, so that
%   matching works regardless of which convention the table was built with.
    switch stimKey
        case "RG";    names = ["RG",   "rectGrid"];
        case "MB";    names = ["MB",   "linearlyMovingBall"];
        case "MBR";   names = ["MBR",  "linearlyMovingBar"];
        case "SDGs";  names = ["SDGs", "StaticGrating"];
        case "SDGm";  names = ["SDGm", "MovingGrating"];
        case "NI";    names = ["NI",   "naturalImage"];
        case "NV";    names = ["NV",   "naturalMovie"];
        case "FFF";   names = ["FFF",  "fullFieldFlash"];
        otherwise;    names = string(stimKey);                           % fallback: just the abbreviation
    end
end

function fName = levelToFieldName(catName, value)
% levelToFieldName  Build a valid MATLAB field name matching
%   StatisticsPerNeuronPerCategory's naming convention.
%   e.g. ('direction', 0)   -> 'direction_0'
%        ('size', 5)        -> 'size_5'
%        ('speed', 0.3)     -> 'speed_0p3'
%        ('offset', -1)     -> 'offset_neg1'
    fName = sprintf('%s_%g', lower(strtrim(catName)), value);            % base: 'cat_value'
    fName = strrep(fName, '.', 'p');                                     % decimal point -> 'p'
    fName = strrep(fName, '-', 'neg');                                   % minus sign    -> 'neg'
end