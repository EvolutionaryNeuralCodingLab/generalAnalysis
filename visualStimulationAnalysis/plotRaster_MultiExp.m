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
%   BUG 1  (fixed earlier) — Sort-by-peak double-smoothing.
%   BUG 2  (fixed earlier) — Wrong colour limits for zScore + gray.
%   BUG 3  (fixed earlier) — Stim-offset xline in wrong units.
%   FIX 4  — Trial selection now uses post-onset window only, avoiding
%            bias toward high-baseline trials.
%   FIX 5  — Zero-SD neurons are counted and logged instead of silently
%            dropped.
%   FIX 6  — TakeTopPercentTrials = 0 is handled explicitly.
%   FIX 7  — Baseline/binWidth alignment assertion added.
%   FIX 8  — Diverging colormap warns when climNeg = 0.
%   FIX 9  — Accumulator alignment assertion after experiment loop.
%   FIX A  — depthFile is now loaded outside the forloop guard, preventing
%            a crash when reloading from cache with sortBy="depth".
%   FIX B  — Figure size override removed; single nStim-scaled sizing.
%   FIX C  — printFig filename now joins splitCategory safely.
%   FIX D  — Cache key includes numel(exList) to reduce collision risk.
%   FIX E  — smoothdata window parameter documented correctly.
%   FIX F  — Dead gt assignment removed from main-loop switch preamble.
%   NEW    — sortBy = "spatialTuning": sort neurons by a column from an
%            external spatial-tuning table, matched via Phy cluster ID.
%   NEW    — phyAll accumulator tracks Phy cluster IDs for each neuron row.
%   NEW    — sortBy = "preferredCategory": group neurons by their preferred
%            category level (e.g. direction, size), show only
%            preferred-level trials in each PSTH row, secondary sort by
%            mean post-onset firing rate within each group.
%   NEW    — unionUnits mode: neurons are matched across all stimulus
%            panels so that row k is the same physical neuron everywhere.
%            The union of per-stim responsive sets is used; all stim types
%            must be present for each included experiment.
%   NEW    — "SDG" combined stimulus type: single-panel raster that shows
%            the static grating phase as the main stimulus window, with the
%            moving phase appearing as post-stimulus time.  A second xline
%            marks the static-to-moving transition.

arguments
    exList double                                                        % vector of experiment IDs to include
    params.stimTypes       (1,:) string  = ["RG","MB"]                   % stimulus abbreviations — one tile each (RG|MB|MBR|SDGs|SDGm|SDG|NI|NV|FFF)
    params.binWidth        double        = 10                            % PSTH bin width in ms
    params.smooth          double        = 0                             % Gaussian smoothing SD in ms (0 = off)
    params.statType        string        = "MaxPermuteTest"              % which statistics field to use
    params.speed           string        = "max"                         % ball-speed selector: "max" or other
    params.alpha           double        = 0.05                          % significance threshold
    params.postStim        double        = 0                             % post-stimulus window in ms
    params.preBase         double        = 200                           % pre-stimulus baseline in ms
    params.overwrite       logical       = false                         % if true, recompute even if cache exists
    params.TakeTopPercentTrials double   = 1                             % fraction (0,1] of trials to keep; [] or 0 = keep all
    params.zScore          logical       = true                          % z-score each neuron using its baseline
    params.sortBy          string        = "none"                        % "peak" | "depth" | "spatialTuning" | "preferredCategory" | "none"
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
    % --- Union responsive units mode ---
    params.unionUnits      logical       = false                         % if true, match neurons across stim panels (same row = same neuron)
end

% -------------------------------------------------------------------------
% Sanity check: baseline must be an integer multiple of bin width
% -------------------------------------------------------------------------
assert(mod(params.preBase, params.binWidth) == 0, ...                    % prevents misaligned baseline bin edges
    'preBase (%g ms) must be a multiple of binWidth (%g ms).', ...
    params.preBase, params.binWidth);

% -------------------------------------------------------------------------
% Validate unionUnits constraints
% -------------------------------------------------------------------------
if params.unionUnits
    % preferredCategory sorts differently per stimulus, which would break
    % the row-alignment guarantee.  Block this combination up front.
    assert(params.sortBy ~= "preferredCategory", ...
        'unionUnits = true is incompatible with sortBy = "preferredCategory" (per-stim sort would break row alignment).');

    % Union mode also requires that all stim types share the same neuron
    % pool, so per-stim category splitting is not meaningful.
    assert(all(strlength(params.splitCategory) == 0) || ...
           (isscalar(params.splitCategory) && params.splitCategory == ""), ...
        'unionUnits = true is incompatible with splitCategory (different filtering per stim breaks row alignment).');

    % Union across stim types needs ≥2 types.  Exception: a single "SDG"
    % tile unions across the static and moving *phases* internally.
    isSingleSDG = numel(params.stimTypes) == 1 && params.stimTypes == "SDG";
    assert(numel(params.stimTypes) >= 2 || isSingleSDG, ...
        'unionUnits = true requires at least 2 stimTypes, or a single "SDG" (phase-union).');
end

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
% (FIX A: loaded unconditionally outside forloop to prevent scope crash
%  when reloading from cache)
% -------------------------------------------------------------------------
depthFile  = '';                                                         % initialise empty; only populated if needed
depthTable = [];                                                         % initialise empty
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
p = [p 'lizards'];                                                      % include 'lizards' folder

if ~exist([p '\Combined_lizard_analysis'], 'dir')                        % create output dir if absent
    cd(p)                                                                % change to parent dir
    mkdir Combined_lizard_analysis                                       % create sub-directory
end
saveDir = [p '\Combined_lizard_analysis'];                               % output directory

stimLabel  = strjoin(params.stimTypes, '-');                             % e.g. "RG-MB"

% --- Include splitCategory in cache filename to avoid collisions ---
% Join all per-stim categories into a compact suffix, e.g. "_cat-size-direction"
if params.sortBy == "preferredCategory"
    nonEmpty = params.splitCategory(strlength(params.splitCategory) > 0); % only non-empty categories
    catSuffix = sprintf('_cat-%s', strjoin(nonEmpty, '-'));               % e.g. "_cat-size-direction"
else
    catSuffix = '';                                                       % no suffix for other sort modes
end

% --- Include union flag in cache filename ---
if params.unionUnits
    unionSuffix = '_union';                                              % distinguish union cache from standard
else
    unionSuffix = '';                                                     % no suffix for standard mode
end

% FIX D: include numel(exList) to reduce cache key collisions when
% experiment lists share the same first/last ID but differ in between.
nameOfFile = sprintf('\\Ex_%d-%d_n%d_Raster_%s%s%s.mat', ...
    exList(1), exList(end), numel(exList), stimLabel, catSuffix, unionSuffix); % cache filename

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
    prefLevelAll = cell(1, nStim);                                       % preferred category level per neuron
    respGroupAll = cell(1, nStim);                                       % union group label per neuron ("both", "X only", etc.)

    for s = 1:nStim
        rasterAll{s}    = [];                                            % initialise empty PSTH matrix
        depthAll{s}     = [];                                            % initialise empty depth vector
        expAll{s}       = [];                                            % initialise empty exp-ID vector
        phyAll{s}       = [];                                            % initialise empty Phy-ID vector
        prefLevelAll{s} = [];                                            % initialise empty preferred-level vector
        respGroupAll{s} = string.empty;                                  % initialise empty string vector
    end

    % Counter for neurons dropped due to zero-SD baseline
    nDroppedZeroSD = zeros(1, nStim);                                    % per-stim counter

    % Counter for experiments skipped due to insufficient category levels
    nSkippedCat = zeros(1, nStim);                                       % per-stim counter

    % Counter for experiments skipped in union mode (missing stim type)
    nSkippedUnion = 0;                                                   % scalar counter

    % --- Shared time-axis variables ---
    lockedPreBase = [];                                                  % baseline duration (ms) — locked on first exp
    lockedEdges   = cell(1, nStim);                                      % bin edges (ms) per stimulus
    lockedNBins   = zeros(1, nStim);                                     % bin count per stimulus
    tAxis         = cell(1, nStim);                                      % left bin edge (ms) per stimulus
    stimDurAll    = zeros(1, nStim);                                     % stim duration (ms) per stim for xline
    movingOnsetAll = nan(1, nStim);                                      % ms from stim onset to moving phase (SDG only; NaN elsewhere)

    % -----------------------------------------------------------------
    % Pre-scan: find shortest stimulus duration per type across all exps
    % -----------------------------------------------------------------
    minStimDur     = inf(1, nStim);                                      % initialise with inf
    minMovingOnset = inf(1, nStim);                                      % SDG only: shortest static phase across exps (= transition point)

    for ei = 1:nExp
        for s = 1:nStim
            try
                NPtmp = loadNPclassFromTable(exList(ei));                % load NP object
                stimKey = params.stimTypes(s);                           % current stim abbreviation

                % --- Build analysis object from abbreviation ---
                switch stimKey
                    case "RG";   objTmp = rectGridAnalysis(NPtmp);                 % receptive-field grid
                    case "MB";   objTmp = linearlyMovingBallAnalysis(NPtmp);       % moving ball
                    case "MBR";  objTmp = linearlyMovingBarAnalysis(NPtmp);        % moving bar
                    case {"SDGs","SDGm","SDG"}
                        objTmp = StaticDriftingGratingAnalysis(NPtmp);             % grating (all variants)
                    case "NI";   objTmp = imageAnalysis(NPtmp);                    % natural image
                    case "NV";   objTmp = movieAnalysis(NPtmp);                    % natural movie
                    case "FFF";  objTmp = fullFieldFlashAnalysis(NPtmp);           % full-field flash
                    otherwise;   error('Unknown stimType abbreviation: %s', stimKey);
                end
                NRtmp = objTmp.ResponseWindow;                           % response-window struct

                % Resolve fieldName (same logic as main loop)
                fn = resolveFieldName(stimKey, params.speed);            % sub-field key for this stim type

                try
                    dur = NRtmp.(fn).stimDur;                            % sub-field duration
                catch
                    dur = NRtmp.stimDur;                                 % flat struct fallback
                end

                % For combined SDG, the stimulus spans both phases:
                % total duration = static phase + moving phase.
                % Also track the static duration separately as the
                % transition point for the xline.
                if stimKey == "SDG"
                    staticDur_ms  = objTmp.VST.static_time * 1000;       % static phase duration in ms
                    movingDur_ms  = NRtmp.Moving.stimDur;                % moving phase duration in ms
                    dur           = staticDur_ms + movingDur_ms;         % total stimulus = both phases
                    minMovingOnset(s) = min(minMovingOnset(s), staticDur_ms); % track shortest static phase (= transition point)
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

        % =============================================================
        % UNION MODE: pre-scan all stim types for this experiment
        % =============================================================
        % When unionUnits is true we need to:
        %   (a) verify every stim type can be loaded for this experiment,
        %   (b) find responsive units per stim,
        %   (c) compute the union.
        % If any stim type fails to load, skip the entire experiment.

        % Pre-allocate per-stim containers for this experiment.
        % These are populated in the pre-scan (union mode) or in the
        % main per-stim loop (standard mode).
        expObjs        = cell(1, nStim);                                 % analysis objects
        expGoodU       = cell(1, nStim);                                 % good-unit identity matrices
        expPsort       = cell(1, nStim);                                 % p_sort structs
        expGoodPhyIDs  = cell(1, nStim);                                 % Phy cluster IDs for good units
        expFieldNames  = strings(1, nStim);                              % resolved sub-field names
        expStartStim   = zeros(1, nStim);                                % ms offset for stim onset
        expStats       = cell(1, nStim);                                 % statistics structs
        expPvals       = cell(1, nStim);                                 % p-value vectors
        expC           = cell(1, nStim);                                 % condition matrices
        expDirectimes  = cell(1, nStim);                                 % onset times (ms)
        expResponsive  = cell(1, nStim);                                 % responsive-unit index vectors

        skipExperiment = false;                                          % flag: skip if any stim fails (union mode)

        for s = 1:nStim

            stimType = params.stimTypes(s);                              % current stim abbreviation

            % --- Build stimulus-specific analysis object ---
            try
                switch stimType
                    case "RG"
                        obj = rectGridAnalysis(NP);                      % receptive-field grid stimulus
                    case "MB"
                        obj = linearlyMovingBallAnalysis(NP);            % moving ball stimulus
                    case "MBR"
                        obj = linearlyMovingBarAnalysis(NP);             % moving bar stimulus
                    case {"SDGs","SDGm","SDG"}
                        obj = StaticDriftingGratingAnalysis(NP);         % grating stimulus (all variants)
                    case "NI"
                        obj = imageAnalysis(NP);                         % natural image stimulus
                    case "NV"
                        obj = movieAnalysis(NP);                         % natural movie stimulus
                    case "FFF"
                        obj = fullFieldFlashAnalysis(NP);                % full-field flash stimulus
                    otherwise
                        error('Unknown stimType abbreviation: %s', stimType);
                end
                NeuronResp = obj.ResponseWindow;                         % response-window struct
            catch ME
                warning('Could not build %s for exp %d: %s', stimType, ex, ME.message);
                if params.unionUnits
                    skipExperiment = true;                                % union mode: skip entire experiment
                    break                                                % exit inner loop immediately
                end
                continue                                                 % standard mode: skip this stim/exp
            end

            expObjs{s} = obj;                                            % store analysis object

            % --- Select statistics struct ---
            if params.statType == "BootstrapPerNeuron"
                Stats = obj.BootstrapPerNeuron;                          % bootstrap p-values
            else
                Stats = obj.StatisticsPerNeuron;                         % default: permutation test
            end
            expStats{s} = Stats;                                         % store for later use

            % --- Resolve sub-field name and stim-onset offset ---
            fieldName = resolveFieldName(stimType, params.speed);        % sub-field key for this stim type
            expFieldNames(s) = fieldName;                                % store resolved name

            startStim = 0;                                               % ms offset for stim onset
            if stimType == "SDGm"
                startStim = obj.VST.static_time * 1000;                  % moving phase onset (s -> ms)
            end
            % SDG combined: onset is start of static phase, so startStim = 0.
            expStartStim(s) = startStim;                                 % store onset offset

            % --- Convert Phy sorting to tIc format ---
            p_sort = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
            label  = string(p_sort.label');                              % quality label per unit
            goodU  = p_sort.ic(:, label == 'good');                      % keep only curated 'good' units
            expGoodU{s}  = goodU;                                        % store good-unit matrix
            expPsort{s}  = p_sort;                                       % store full p_sort struct

            % --- Extract Phy cluster IDs for good units ---
            goodPhyIDs = p_sort.phy_ID(label == 'good');                 % Phy cluster IDs matching goodU columns
            expGoodPhyIDs{s} = goodPhyIDs;                               % store Phy IDs

            % --- General response p-values (StatisticsPerNeuron) ---
            try
                pvals = Stats.(fieldName).pvalsResponse;                 % stim-specific p-values
            catch
                pvals = Stats.pvalsResponse;                             % flat struct fallback
            end
            expPvals{s} = pvals;                                         % store p-value vector

            % --- Stimulus onset times and condition matrix ---
            try
                C = NeuronResp.(fieldName).C;                            % condition matrix (sub-field)
            catch
                C = NeuronResp.C;                                        % condition matrix (flat)
            end
            expC{s} = C;                                                 % store condition matrix
            expDirectimes{s} = C(:, 1)' + startStim;                    % onset times in ms

            % --- Find responsive neurons (general filter) ---
            expResponsive{s} = find(pvals < params.alpha);               % indices of responsive units

        end % stim pre-scan loop

        % --- Union mode: skip experiment if any stim failed ---
        if params.unionUnits && skipExperiment
            nSkippedUnion = nSkippedUnion + 1;                           % count skipped experiment
            fprintf('  Skipping experiment %d (missing stim type in union mode).\n', ex);
            continue                                                     % skip to next experiment
        end

        % --- Union mode: verify shared neuron pool and compute union ---
        % Also classify each neuron into responsiveness groups:
        %   - For single SDG:  "both", "static", "moving"
        %   - For multi-stim:  "both", "<stimA>", "<stimB>"
        expNeuronGroup = strings(1, 0);                                  % per-unit group label (populated below if union mode)

        if params.unionUnits

            isSingleSDG = nStim == 1 && params.stimTypes == "SDG";       % special case: phase-union

            if ~isSingleSDG
                % --- Multi-stim: verify shared neuron pool ---
                refPhyIDs = sort(expGoodPhyIDs{1});                      % reference: good Phy IDs from first stim
                for s = 2:nStim
                    if ~isequal(sort(expGoodPhyIDs{s}), refPhyIDs)
                        warning(['Experiment %d: good-unit Phy IDs differ between %s and %s. ' ...
                                 'Spike sorting may differ. Skipping experiment.'], ...
                            ex, params.stimTypes(1), params.stimTypes(s));
                        skipExperiment = true;                           % flag to skip this experiment
                        break
                    end
                end
                if skipExperiment
                    nSkippedUnion = nSkippedUnion + 1;                   % count skipped
                    continue                                             % skip to next experiment
                end
            end

            % =============================================================
            % Compute union and classify neurons into groups
            % =============================================================
            nGoodUnits = numel(expPvals{1});                             % number of good units (shared across stim types)

            if isSingleSDG
                % --- SDG phase-union: union across static and moving phases ---
                % expPvals{1} holds Static p-values (from resolveFieldName).
                % Also extract Moving p-values from the same analysis object.
                staticPvals = expPvals{1};                               % already stored by pre-scan
                try
                    movingPvals = expStats{1}.Moving.pvalsResponse;      % moving-phase p-values
                catch
                    movingPvals = expStats{1}.pvalsResponse;             % flat struct fallback
                end

                staticResp = staticPvals(:) < params.alpha;              % logical: responsive to static phase
                movingResp = movingPvals(:) < params.alpha;              % logical: responsive to moving phase
                unionMask  = staticResp | movingResp;                    % union: responsive to either phase
                unionIdx   = find(unionMask);                            % indices of union neurons

                % Classify each good unit (indexed by position in goodU)
                expNeuronGroup = strings(1, nGoodUnits);                 % pre-allocate
                expNeuronGroup(staticResp & movingResp)  = "both";       % responsive to both phases
                expNeuronGroup(staticResp & ~movingResp) = "static";    % static phase only
                expNeuronGroup(~staticResp & movingResp) = "moving";    % moving phase only
                % Units not in the union keep "" (will never be queried)

                fprintf('  SDG phase-union exp %d: %d static, %d moving, %d both, %d total union.\n', ...
                    ex, sum(staticResp & ~movingResp), sum(~staticResp & movingResp), ...
                    sum(staticResp & movingResp), numel(unionIdx));

            else
                % --- Multi-stim union: union across stimulus types ---
                % Build a logical matrix: nGoodUnits x nStim
                respMatrix = false(nGoodUnits, nStim);                   % pre-allocate
                for s = 1:nStim
                    respMatrix(expResponsive{s}, s) = true;              % mark responsive units per stim
                end

                unionMask = any(respMatrix, 2);                          % responsive to at least one stim
                unionIdx  = find(unionMask);                             % indices of union neurons

                % Classify each unit
                nRespStim = sum(respMatrix, 2);                          % how many stim types each unit responds to
                expNeuronGroup = strings(1, nGoodUnits);                 % pre-allocate

                % Neurons responsive to ALL stim types → "both" (or "all" for 3+)
                expNeuronGroup(nRespStim >= 2) = "both";

                % Neurons responsive to exactly one stim type → "<stimName> only"
                for s = 1:nStim
                    onlyThisStim = respMatrix(:, s) & nRespStim == 1;    % responsive to this stim only
                    expNeuronGroup(onlyThisStim) = params.stimTypes(s); % e.g. "MB"
                end

                fprintf('  Union exp %d: %d both, ', ex, sum(nRespStim >= 2));
                for s = 1:nStim
                    fprintf('%d %s only, ', sum(respMatrix(:,s) & nRespStim==1), params.stimTypes(s));
                end
                fprintf('%d total union.\n', numel(unionIdx));
            end

            if isempty(unionIdx)
                fprintf('  No responsive neurons in any stim type/phase for exp %d.\n', ex);
                continue                                                 % nothing to add
            end

            % Override each stim's responsive set with the shared union
            for s = 1:nStim
                expResponsive{s} = unionIdx;                             % same neuron list for every panel
            end
        end

        % =============================================================
        % Per-stim PSTH building (shared by both modes)
        % =============================================================
        for s = 1:nStim

            stimType  = params.stimTypes(s);                             % current stim abbreviation
            fieldName = expFieldNames(s);                                % resolved sub-field name
            obj       = expObjs{s};                                      % analysis object for this stim

            if isempty(obj)
                continue                                                 % object failed to load (standard mode: already warned)
            end

            % Retrieve pre-computed variables from the pre-scan
            goodU          = expGoodU{s};                                % good-unit identity matrix
            p_sort         = expPsort{s};                                % full p_sort struct
            goodPhyIDs     = expGoodPhyIDs{s};                           % Phy cluster IDs
            pvals          = expPvals{s};                                % p-value vector
            directimesSorted = expDirectimes{s};                         % onset times (ms)
            C              = expC{s};                                    % condition matrix

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
                NeuronResp = obj.ResponseWindow;                         % re-fetch for column name access
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
                % parameters matching catColNames.  catIdx is 1-based into
                % catColNames.  Offset by 1 for the time column.
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
                windowTotal   = preBase + rawStimDur_ms + params.postStim; % full window in ms
            else
                rawStimDur_ms = params.postStim;                         % fixed window
                windowTotal   = preBase + params.postStim;               % full window in ms
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

                % Store moving-phase onset for SDG combined type.
                % Use the pre-scanned minimum static duration across
                % experiments so the transition xline is consistent
                % with the truncated window.
                if stimType == "SDG"
                    movingOnsetAll(s) = minMovingOnset(s);               % static-to-moving transition (ms)
                end

                fprintf('  [%s] Locked window: preBase=%d ms, stimDur=%.0f ms, nBins=%d\n', ...
                    stimType, lockedPreBase, rawStimDur_ms, lockedNBins(s));
            end

            % --- Find responsive neurons ---
            % In union mode, expResponsive was already overwritten with
            % the shared union set.  In standard mode it holds the
            % per-stim responsive indices.
            if isCatMode && ~params.useGeneralFilter && ~params.unionUnits
                eNeurons = find(orMask);                                 % per-category OR filter
            else
                eNeurons = expResponsive{s};                             % general filter or union set
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
                % most stable baseline estimate.
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
                        % In union mode, we MUST append a NaN row to
                        % preserve row alignment.
                        if params.unionUnits
                            rasterAll{s}    = [rasterAll{s}; nan(1, lockedNBins(s))]; % NaN placeholder row
                            phyAll{s}(end+1)    = goodPhyIDs(u);         % Phy ID (for alignment tracking)
                            expAll{s}(end+1)    = ex;                    % experiment ID
                            prefLevelAll{s}(end+1) = NaN;                % no preferred level
                            depthAll{s}(end+1)  = NaN;                   % no depth
                            respGroupAll{s}(end+1) = expNeuronGroup(u);  % union group label
                            nAppended = nAppended + 1;                   % count as appended (placeholder)
                        end
                        continue                                         % skip neuron
                    end
                end

                % ==========================================================
                % PREFERRED-CATEGORY MODE: determine preferred level and
                % build PSTH from only that level's trials
                % ==========================================================
                if isCatMode

                    % --- Compute mean post-onset rate per level ---
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
                        % In union mode, append NaN row to preserve alignment
                        if params.unionUnits
                            rasterAll{s}    = [rasterAll{s}; nan(1, lockedNBins(s))]; % NaN placeholder row
                            phyAll{s}(end+1)    = goodPhyIDs(u);
                            expAll{s}(end+1)    = ex;
                            prefLevelAll{s}(end+1) = NaN;
                            depthAll{s}(end+1)  = NaN;
                            respGroupAll{s}(end+1) = expNeuronGroup(u);  % union group label
                            nAppended = nAppended + 1;
                        end
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
                        % In union mode, append NaN row to preserve alignment
                        if params.unionUnits
                            rasterAll{s}    = [rasterAll{s}; nan(1, lockedNBins(s))]; % NaN placeholder row
                            phyAll{s}(end+1)    = goodPhyIDs(u);
                            expAll{s}(end+1)    = ex;
                            prefLevelAll{s}(end+1) = prefLevel;
                            depthAll{s}(end+1)  = NaN;
                            respGroupAll{s}(end+1) = expNeuronGroup(u);  % union group label
                            nAppended = nAppended + 1;
                        end
                        continue                                         % skip: undefined z-score
                    end
                end

                % --- Smooth PSTH if requested ---
                if params.smooth > 0
                    % FIX E: smoothdata('gaussian', N) treats N as window
                    % WIDTH, not SD.  Convert: width ≈ 6*SD ensures ~99.7%
                    % of the kernel mass is captured.
                    smoothSD    = params.smooth / params.binWidth;        % smoothing SD in bin units
                    smoothWidth = max(3, round(6 * smoothSD));           % window width in bins (at least 3)
                    neuronPSTH  = smoothdata(neuronPSTH, 'gaussian', smoothWidth); % smooth in-place
                end

                % --- Append neuron row to accumulators ---
                rasterAll{s}    = [rasterAll{s}; neuronPSTH];           % PSTH row
                phyAll{s}(end+1)    = goodPhyIDs(u);                    % Phy cluster ID for this unit
                expAll{s}(end+1)    = ex;                                % source experiment
                prefLevelAll{s}(end+1) = prefLevel;                     % preferred level (NaN if not catMode)
                % Store union responsiveness group label
                if params.unionUnits && numel(expNeuronGroup) >= u
                    respGroupAll{s}(end+1) = expNeuronGroup(u);          % e.g. "both", "static", "MB"
                else
                    respGroupAll{s}(end+1) = "";                         % not in union mode or not classified
                end
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

    % Report experiments skipped in union mode
    if params.unionUnits && nSkippedUnion > 0
        fprintf('  Union mode: skipped %d experiment(s) with missing stim types.\n', nSkippedUnion);
    end

    % FIX 9: verify accumulator alignment after all experiments
    for s = 1:nStim
        nRows = size(rasterAll{s}, 1);                                   % number of neuron rows
        assert(numel(expAll{s})       == nRows, 'expAll{%d} length mismatch.', s);
        assert(numel(phyAll{s})       == nRows, 'phyAll{%d} length mismatch.', s);
        assert(numel(depthAll{s})     == nRows, 'depthAll{%d} length mismatch.', s);
        assert(numel(prefLevelAll{s}) == nRows, 'prefLevelAll{%d} length mismatch.', s);
        if params.unionUnits
            assert(numel(respGroupAll{s}) == nRows, 'respGroupAll{%d} length mismatch.', s);
        end
    end

    % UNION MODE: verify row counts match across all stim types.
    % This is the fundamental invariant: same row = same neuron.
    if params.unionUnits
        refRows = size(rasterAll{1}, 1);                                 % row count for first stim
        for s = 2:nStim
            assert(size(rasterAll{s}, 1) == refRows, ...
                'Union mode row-count mismatch: stim %d has %d rows, stim 1 has %d.', ...
                s, size(rasterAll{s}, 1), refRows);
        end

        % Also verify Phy IDs match across panels (same neuron in each row)
        for s = 2:nStim
            if ~isequal(phyAll{s}, phyAll{1}) || ~isequal(expAll{s}, expAll{1})
                error(['Union mode: Phy IDs or experiment IDs differ between stim 1 and stim %d. ' ...
                       'Row alignment is broken.'], s);
            end
        end
        fprintf('  Union mode verified: %d neurons, row alignment intact.\n', refRows);
    end

    % ------------------------------------------------------------------
    % Save processed data to disk
    % ------------------------------------------------------------------
    S.expList         = exList;                                          % for validation on reload
    S.lockedEdges     = lockedEdges;                                     % per-stim bin edges
    S.lockedPreBase   = lockedPreBase;                                   % shared baseline
    S.stimDurAll      = stimDurAll;                                      % per-stim duration for xline
    S.movingOnsetAll  = movingOnsetAll;                                  % static-to-moving transition (SDG only)
    S.params          = params;                                          % full parameter set

    for s = 1:numel(params.stimTypes)
        stimField = matlab.lang.makeValidName(params.stimTypes(s));      % valid struct field name
        S.(sprintf('%s_raster',    stimField)) = rasterAll{s};           % PSTH matrix
        S.(sprintf('%s_depth',     stimField)) = depthAll{s};            % depth vector
        S.(sprintf('%s_exp',       stimField)) = expAll{s};              % experiment ID vector
        S.(sprintf('%s_phy',       stimField)) = phyAll{s};              % Phy cluster ID vector
        S.(sprintf('%s_prefLevel', stimField)) = prefLevelAll{s};        % preferred level vector
        S.(sprintf('%s_respGroup', stimField)) = respGroupAll{s};        % union responsiveness group label
    end

    save([saveDir nameOfFile], '-struct', 'S');                          % write cache to disk
    fprintf('\nSaved raster data to:\n  %s\n', [saveDir nameOfFile]);

else
    % ------------------------------------------------------------------
    % Reload cached data
    % ------------------------------------------------------------------
    lockedEdges    = S.lockedEdges;                                      % restore bin edges
    lockedPreBase  = S.lockedPreBase;                                    % restore baseline
    stimDurAll     = S.stimDurAll;                                       % restore stim durations

    % Restore movingOnsetAll if present, otherwise NaN (old cache)
    if isfield(S, 'movingOnsetAll')
        movingOnsetAll = S.movingOnsetAll;                               % restore static-to-moving transition
    else
        movingOnsetAll = nan(1, numel(params.stimTypes));                % old cache: no SDG combined data
    end

    rasterAll    = cell(1, numel(params.stimTypes));                     % pre-allocate
    depthAll     = cell(1, numel(params.stimTypes));
    expAll       = cell(1, numel(params.stimTypes));
    phyAll       = cell(1, numel(params.stimTypes));
    prefLevelAll = cell(1, numel(params.stimTypes));
    respGroupAll = cell(1, numel(params.stimTypes));

    for s = 1:numel(params.stimTypes)
        stimField    = matlab.lang.makeValidName(params.stimTypes(s));   % valid struct field name
        rasterAll{s} = S.(sprintf('%s_raster', stimField));              % restore PSTH matrix
        depthAll{s}  = S.(sprintf('%s_depth',  stimField));              % restore depth vector
        expAll{s}    = S.(sprintf('%s_exp',    stimField));              % restore experiment IDs

        % phyAll may be absent in old caches — require recompute if needed
        phyField = sprintf('%s_phy', stimField);
        if isfield(S, phyField)
            phyAll{s} = S.(phyField);                                    % restore Phy IDs
        elseif params.sortBy == "spatialTuning" || params.sortBy == "preferredCategory" || params.unionUnits
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

        % respGroupAll may be absent in old caches
        rgField = sprintf('%s_respGroup', stimField);
        if isfield(S, rgField)
            respGroupAll{s} = S.(rgField);                               % restore group labels
        elseif params.unionUnits
            error(['Cache file lacks respGroup data (old format). ' ...
                   'Re-run with params.overwrite = true.']);
        else
            respGroupAll{s} = repmat("", 1, size(rasterAll{s}, 1));      % fill empty if unused
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

        legacyNames = abbrevToLegacyNames(params.stimTypes(s));          % e.g. ["MB","linearlyMovingBall"]
        stimMask = ismember(string(tuningTable.stimulus), legacyNames);  % match any known name variant
        subT     = tuningTable(stimMask, :);                             % subtable for this stim

        for k = 1:nNeu
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
% In union mode, sorting must be computed ONCE and applied identically
% to all stimulus panels to preserve row alignment.
% Strategy: compute sortIdx from the FIRST non-empty stim (or a combined
% metric), then apply it to all panels.

levelGroupBounds = cell(1, numel(params.stimTypes));                     % {s}: struct with edges and labels

if params.unionUnits
    % ------------------------------------------------------------------
    % UNION SORT: primary sort by responsiveness group, secondary sort
    % by sortBy within each group.  One sort order applied to all panels.
    % ------------------------------------------------------------------
    % Group order: "both" first, then individual-stim/phase groups
    % (alphabetical).  Within each group, neurons are sorted by the
    % selected sortBy criterion (peak, depth, spatialTuning, or identity).

    refStim  = find(~cellfun(@isempty, rasterAll), 1);                   % first non-empty stim index
    nNeurons = size(rasterAll{refStim}, 1);                              % shared neuron count
    groups   = respGroupAll{refStim};                                    % group label per neuron (shared across panels)

    % Determine unique groups, with "both" forced to position 1
    uniqueGroups = unique(groups);                                       % alphabetical order
    uniqueGroups = uniqueGroups(uniqueGroups ~= "");                     % remove empty strings (shouldn't occur)
    hasBoth = ismember("both", uniqueGroups);                            % check if "both" is present
    if hasBoth
        uniqueGroups(uniqueGroups == "both") = [];                       % remove "both" temporarily
        uniqueGroups = ["both", sort(uniqueGroups)];                     % prepend "both", sort the rest
    else
        uniqueGroups = sort(uniqueGroups);                               % just sort alphabetically
    end

    % --- Compute secondary sort metric (within-group ordering) ---
    % Pre-compute whatever the sortBy criterion needs, so we can apply
    % it independently within each group.
    if params.sortBy == "peak"
        % Average PSTH across panels for peak detection
        avgData = zeros(nNeurons, size(rasterAll{refStim}, 2));          % accumulator
        nContrib = 0;                                                    % contributing panel count
        for s = 1:numel(params.stimTypes)
            if ~isempty(rasterAll{s})
                nBins = min(size(avgData, 2), size(rasterAll{s}, 2));    % common bin count
                avgData(:, 1:nBins) = avgData(:, 1:nBins) + rasterAll{s}(:, 1:nBins);
                nContrib = nContrib + 1;
            end
        end
        avgData = avgData / max(nContrib, 1);                            % mean across panels

        postStimBins = tAxis{refStim} >= lockedPreBase;                  % post-onset bin mask

        if size(avgData, 2) > 100
            dataForSort = ConvBurstMatrix( ...
                avgData, fspecial('gaussian', [1 params.GaussianLength+20], 2), 'same');
        else
            dataForSort = ConvBurstMatrix( ...
                avgData, fspecial('gaussian', [1 params.GaussianLength], 2), 'same');
        end

        [~, secondaryMetric] = max(dataForSort(:, postStimBins), [], 2); % peak bin per neuron
        secondaryDirection   = 'ascend';                                 % early-peaking first

    elseif params.sortBy == "depth"
        secondaryMetric    = depthAll{refStim}(:);                       % depth per neuron
        secondaryDirection = 'ascend';                                   % shallowest first

    elseif params.sortBy == "spatialTuning"
        secondaryMetric    = tuningAll{refStim}(:);                      % tuning index per neuron
        secondaryDirection = char(params.tuningSortOrder);                % user-specified order

    else
        secondaryMetric    = (1:nNeurons)';                              % identity (preserve original order)
        secondaryDirection = 'ascend';
    end

    % --- Build sortIdx: iterate groups, sort within each ---
    sortIdx   = zeros(1, nNeurons);                                      % pre-allocate
    cursor    = 0;                                                       % running position
    groupInfo = struct('edges', [], 'labels', {{}});                     % for plot annotations

    for gi = 1:numel(uniqueGroups)
        groupMask    = groups == uniqueGroups(gi);                       % neurons in this group
        groupIndices = find(groupMask);                                  % their row indices
        nGroup       = numel(groupIndices);                              % count

        if nGroup == 0; continue; end                                    % skip empty groups

        % Secondary sort within group
        groupMetric = secondaryMetric(groupIndices);                     % metric for this group's neurons
        [~, withinOrder] = sort(groupMetric, secondaryDirection, 'MissingPlacement', 'last');

        sortIdx(cursor+1 : cursor+nGroup) = groupIndices(withinOrder);   % fill sorted indices

        % Record group boundary for plotting
        groupInfo.edges(end+1)  = cursor + nGroup;                       % last row of this group
        groupInfo.labels{end+1} = char(uniqueGroups(gi));                % group label (e.g. "both", "static")

        cursor = cursor + nGroup;                                        % advance cursor
    end

    % Store group boundaries for all panels (identical since union mode)
    for s = 1:numel(params.stimTypes)
        levelGroupBounds{s} = groupInfo;                                 % reuse same struct for plot annotation
    end

    % Apply the SAME sort to ALL panels
    for s = 1:numel(params.stimTypes)
        if isempty(rasterAll{s}); continue; end
        rasterAll{s}    = rasterAll{s}(sortIdx, :);                      % reorder PSTH rows
        depthAll{s}     = depthAll{s}(sortIdx);                          % reorder depths
        expAll{s}       = expAll{s}(sortIdx);                            % reorder experiment IDs
        phyAll{s}       = phyAll{s}(sortIdx);                            % reorder Phy IDs
        prefLevelAll{s} = prefLevelAll{s}(sortIdx);                      % reorder preferred levels
        respGroupAll{s} = respGroupAll{s}(sortIdx);                      % reorder group labels
        if params.sortBy == "spatialTuning"
            tuningAll{s} = tuningAll{s}(sortIdx);                        % keep tuning vector aligned
        end
    end

else
    % ------------------------------------------------------------------
    % STANDARD SORT: independent per-stim sort order
    % ------------------------------------------------------------------
    for s = 1:numel(params.stimTypes)

        data = rasterAll{s};                                             % nNeurons x nBins
        if isempty(data); continue; end

        if params.sortBy == "peak"
            postStimBins = tAxis{s} >= lockedPreBase;                    % post-onset mask

            % Local smoothed copy for peak detection only
            if size(data, 2) > 100
                dataForSort = ConvBurstMatrix( ...
                    data, fspecial('gaussian', [1 params.GaussianLength+20], 2), 'same');
            else
                dataForSort = ConvBurstMatrix( ...
                    data, fspecial('gaussian', [1 params.GaussianLength], 2), 'same');
            end

            [~, peakBin] = max(dataForSort(:, postStimBins), [], 2);     % peak column per neuron
            [~, sortIdx] = sort(peakBin);                                % early-peaking first

        elseif params.sortBy == "depth"
            [~, sortIdx] = sort(depthAll{s}, 'ascend');                  % shallowest first

        elseif params.sortBy == "spatialTuning"
            [~, sortIdx] = sort(tuningAll{s}, params.tuningSortOrder, ...
                'MissingPlacement', 'last');                              % most-tuned first

        elseif params.sortBy == "preferredCategory"
            if strlength(params.splitCategory(s)) == 0
                sortIdx = 1:size(data, 1);                               % no category: no reorder
            else
                levels   = prefLevelAll{s};                              % preferred level per neuron
                nNeurons = size(data, 1);                                % total neurons for this stim

                postStimBins = tAxis{s} >= lockedPreBase;                % post-onset mask
                meanPostResp = mean(data(:, postStimBins), 2)';          % mean post-onset value per neuron

                uniqueLevels = unique(levels(~isnan(levels)));            % unique levels (ascending)

                sortIdx   = zeros(1, nNeurons);                          % pre-allocate
                cursor    = 0;                                           % running position
                groupInfo = struct('edges', [], 'labels', {{}});         % for plot annotations

                for li = 1:numel(uniqueLevels)
                    groupMask    = levels == uniqueLevels(li);           % neurons preferring this level
                    groupIndices = find(groupMask);                      % their row indices
                    groupResp    = meanPostResp(groupMask);              % their mean responses

                    [~, withinOrder] = sort(groupResp, 'descend');       % strongest first
                    nGroup = numel(groupIndices);                        % neurons in this group

                    sortIdx(cursor+1 : cursor+nGroup) = groupIndices(withinOrder);

                    groupInfo.edges(end+1)  = cursor + nGroup;           % group boundary
                    groupInfo.labels{end+1} = sprintf('%s=%g', ...
                        params.splitCategory(s), uniqueLevels(li));

                    cursor = cursor + nGroup;                            % advance cursor
                end

                % Handle any NaN-level neurons defensively
                nanMask = isnan(levels);
                if any(nanMask)
                    nanIndices = find(nanMask);
                    nNan       = numel(nanIndices);
                    sortIdx(cursor+1 : cursor+nNan) = nanIndices;
                    groupInfo.edges(end+1)  = cursor + nNan;
                    groupInfo.labels{end+1} = 'unclassified';
                    cursor = cursor + nNan;
                end

                levelGroupBounds{s} = groupInfo;                         % store for plotting
            end

        else
            sortIdx = 1:size(data, 1);                                   % identity permutation
        end

        % Apply sort to all parallel vectors
        rasterAll{s}    = data(sortIdx, :);
        depthAll{s}     = depthAll{s}(sortIdx);
        expAll{s}       = expAll{s}(sortIdx);
        phyAll{s}       = phyAll{s}(sortIdx);
        prefLevelAll{s} = prefLevelAll{s}(sortIdx);
        respGroupAll{s} = respGroupAll{s}(sortIdx);                      % reorder group labels (empty in standard mode)

        if params.sortBy == "spatialTuning"
            tuningAll{s} = tuningAll{s}(sortIdx);
        end

    end
end

% =========================================================================
% PLOT
% =========================================================================

% Short display labels for each stimulus type
stimLegendMap = containers.Map( ...
    {'RG',  'MB',  'MBR', 'SDGs', 'SDGm', 'SDG', 'NI',  'NV',  'FFF'}, ...
    {'RG',  'MB',  'MBR', 'SDGs', 'SDGm', 'SDG', 'NI',  'NV',  'FFF'});

nStim = numel(params.stimTypes);

% ------------------------------------------------------------------
% Global colour limits — computed once, shared across all tiles
% ------------------------------------------------------------------
allValues = [];
for s = 1:nStim
    if ~isempty(rasterAll{s})
        % In union mode, NaN rows are placeholders — exclude from colour scaling
        vals = rasterAll{s}(:)';                                         % flatten PSTH matrix
        allValues = [allValues, vals(~isnan(vals))]; %#ok<AGROW>         % pool non-NaN values
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
% FIX B: single figure sizing, width scales with panel count
set(fig, 'Units', 'centimeters', 'Position', [5 5 5*nStim + 2, 10]);

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

    % In union mode, NaN rows render as the lowest colour by default.
    % Set NaN to transparent (white) by using 'AlphaData'.
    if params.unionUnits
        nanMask = all(isnan(data), 2);                                   % rows that are entirely NaN
        if any(nanMask)
            alphaMap = ones(size(data));                                  % fully opaque by default
            alphaMap(nanMask, :) = 0;                                    % make NaN rows transparent
            set(findobj(ax, 'Type', 'Image'), 'AlphaData', alphaMap);    % apply alpha mask
        end
    end

    % ------------------------------------------------------------------
    % Depth-bin boundary lines (depth sort only)
    % ------------------------------------------------------------------
    if params.sortBy == "depth" && ~isempty(depthAll{s}) && ~isempty(depthFile)

        D2 = load(depthFile);                                            % reload bin edges (FIX A: depthFile now always defined)
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
    % Group boundary lines and labels (preferredCategory or union mode)
    % ------------------------------------------------------------------
    showGroupBounds = ~isempty(levelGroupBounds{s}) && ...
        (params.sortBy == "preferredCategory" || params.unionUnits);

    if showGroupBounds

        gInfo   = levelGroupBounds{s};
        nGroups = numel(gInfo.edges);

        ytPositions = zeros(1, nGroups);
        ytLabels    = cell(1, nGroups);

        prevEdge = 0;
        for gi = 1:nGroups
            edgeRow  = gInfo.edges(gi);
            nInGroup = edgeRow - prevEdge;

            ytPositions(gi) = (prevEdge + edgeRow) / 2;                  % centre of group for tick label
            rawLabel = gInfo.labels{gi};                                 % e.g. "size=129" or "both"

            % Format label: for preferredCategory labels contain "=",
            % for union-mode labels are plain strings (e.g. "both").
            if contains(rawLabel, '=')
                levelVal = extractAfter(rawLabel, '=');                   % e.g. "129"
            else
                levelVal = rawLabel;                                     % plain label (e.g. "both", "static")
            end
            ytLabels{gi} = sprintf('%s (n=%d)', levelVal, nInGroup);     % label with count

            if gi < nGroups
                yline(ax, edgeRow + 0.5, 'k-', 'LineWidth', 1.5);       % group boundary line
            end

            prevEdge = edgeRow;
        end

        set(ax, 'YTick', ytPositions, 'YTickLabel', ytLabels, ...
            'TickLength', [0 0], 'YTickLabelRotation', 90);              % vertical tick labels
    end

    % --- Stim onset / offset lines (seconds) ---
    xline(ax, 0,                  'k', 'LineWidth', 1.5);             % onset at t = 0 s
    xline(ax, stimDurAll(s)/1000, 'k', 'LineWidth', 1.5);             % offset in seconds

    % --- SDG combined: additional line at static-to-moving transition ---
    % For the combined SDG type, stimDurAll(s) = static + moving (full
    % stimulus).  The onset xline (t=0) and offset xline (stimDurAll/1000)
    % already bracket the full stimulus.  This third xline marks the
    % boundary between the two phases within the stimulus window.
    if params.stimTypes(s) == "SDG" && ~isnan(movingOnsetAll(s))
        xline(ax, movingOnsetAll(s)/1000, 'r--', 'LineWidth', 1.0 ...
           , 'LabelHorizontalAlignment', 'left', ...
            'FontSize', 6);                                              % red dashed line at static→moving transition
    end

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

sgtitle(sprintf('N = %d experiments', numel(exList)), ...
    'FontName', 'helvetica', 'FontSize', 10);                           % super-title

% --- Export figure for publication if requested ---
if params.PaperFig
    % FIX C: safely join splitCategory into a single string for filename
    if any(strlength(params.splitCategory) > 0)
        catStr = strjoin(params.splitCategory(strlength(params.splitCategory) > 0), '-');
    else
        catStr = 'all';                                                  % fallback when no category is set
    end
    if params.unionUnits
        catStr = [catStr '_union'];                                      % mark union mode in filename
    end
    vs_first.printFig(fig, sprintf('Raster-%s-%s', stimLabel, catStr), PaperFig=params.PaperFig);
end

end

% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function gt = detectGratingType(stimKey)
% detectGratingType  Return the GratingType parameter from a stimulus
%   abbreviation.  'moving' for SDGm, 'static' for SDGs/SDG, '' otherwise.
    switch stimKey
        case "SDGm";           gt = 'moving';                            % moving (drifting) grating
        case {"SDGs", "SDG"};  gt = 'static';                            % static grating (or combined viewed from static onset)
        otherwise;             gt = '';                                   % not a grating stimulus
    end
end

function fn = resolveFieldName(stimKey, speedParam)
% resolveFieldName  Map a stimulus abbreviation + speed selector to the
%   sub-field name used in ResponseWindow / StatisticsPerNeuron structs.
%
%   stimKey    — stimulus abbreviation string (e.g. "MB", "SDGs", "SDG")
%   speedParam — "max" or other speed selector
%   fn         — sub-field key (e.g. 'Speed1', 'Static', '')
    switch stimKey
        case {"MB","MBR"}
            if speedParam == "max"
                fn = 'Speed1';                                           % fastest speed condition
            else
                fn = 'Speed2';                                           % slower speed condition
            end
        case {"SDGs", "SDG"}
            fn = 'Static';                                               % static grating sub-field (SDG combined uses static onset)
        case "SDGm"
            fn = 'Moving';                                               % moving grating sub-field
        otherwise
            fn = '';                                                      % flat struct (RG, NI, NV, FFF)
    end
end

function names = abbrevToLegacyNames(stimKey)
% abbrevToLegacyNames  Return the set of stimulus name strings that might
%   appear in external tables (e.g. tuning tables) for a given abbreviation.
    switch stimKey
        case "RG";    names = ["RG",   "rectGrid"];
        case "MB";    names = ["MB",   "linearlyMovingBall"];
        case "MBR";   names = ["MBR",  "linearlyMovingBar"];
        case "SDGs";  names = ["SDGs", "StaticGrating"];
        case "SDGm";  names = ["SDGm", "MovingGrating"];
        case "SDG";   names = ["SDG",  "SDGs", "StaticGrating"];         % combined maps to static for tuning lookup
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