function plotPSTH_MultiExp(exList, params)
% plotPSTH_MultiExp  Compute and plot population PSTHs across experiments.
%
%   plotPSTH_MultiExp(exList)              — default parameters
%   plotPSTH_MultiExp(exList, Name=Value)  — override any parameter
%
%   Computes a peri-stimulus time histogram for each experiment, then plots
%   the grand-average PSTH ± SEM across experiments.  Supports multiple
%   stimulus types, optional depth-bin stratification, and optional
%   within-stimulus category splits (e.g. one PSTH per ball size).
%
%   STIMULUS TYPE ABBREVIATIONS
%       MB   — linearlyMovingBall   (linearlyMovingBallAnalysis)
%       MBR  — linearlyMovingBar    (linearlyMovingBarAnalysis)
%       RG   — rectGrid             (rectGridAnalysis)
%       SDGm — StaticDriftingGrating, moving phase
%       SDGs — StaticDriftingGrating, static phase
%       NV   — natural video        (movieAnalysis)
%       NI   — natural images       (imageAnalysis)
%       FFF  — fullFieldFlash       (fullFieldFlashAnalysis)
%
%   KEY PARAMETERS
%       stimTypes        — which stimulus analyses to include (abbreviations)
%       requireAllStims  — if true, skip experiments that lack ANY of the
%                          requested stimTypes; ensures a fully matched
%                          population across stim types (recommended for
%                          cross-stim comparisons in publication figures).
%                          NOTE: stimulus PRESENCE is now verified even when
%                          splitBy == "" (the no-split case), and the session
%                          containing each stimulus is discovered across
%                          sessions [0, 1, 2] rather than assuming session 0.
%       splitBy          — category variable to split within each stim type
%                          (e.g. "size", "direction").  "" = no split.
%                          Experiments with <2 levels are automatically skipped.
%       splitLevels      — numeric vector of specific levels to use (e.g. [5 10 20]).
%                          Empty = use all available levels.  Experiments
%                          missing any of the requested levels are skipped.
%       binWidth         — PSTH bin width in ms
%       smooth           — Gaussian smoothing window in ms (0 = none)
%       TakeTopPercentTrials — fraction of trials to keep (1 = all trials)
%       byDepth          — stratify neurons by cortical depth
%
%   See the 'arguments' block below for the full parameter list and defaults.

% -------------------------------------------------------------------------
% Input validation via MATLAB arguments block
% -------------------------------------------------------------------------
arguments
    exList  double                                                         % vector of experiment IDs to include
    params.stimTypes       (1,:) string  = ["RG", "MB"]                   % stimulus types to include (abbreviations)
    params.requireAllStims logical       = false                           % if true, skip experiments that lack ANY of the requested stimTypes (recommended for cross-stim comparisons)
    params.splitBy         string        = ""                              % category variable for within-stim split; "" = no split
    params.splitLevels     double        = []                              % specific category levels to use (e.g. [5 10 20]); empty = all available
    params.binWidth        double        = 50                              % PSTH bin width in ms
    params.smooth          double        = 0                              % Gaussian smoothing window in ms (0 = no smoothing)
    params.statType        string        = "maxPermuteTest"               % statistical test used for per-neuron p-values
    params.speed           string        = "max"                          % speed condition selector for MB/MBR stimuli
    params.alpha           double        = 0.05                           % significance threshold for neuron responsiveness
    params.shadeSTD        logical       = true                           % shade ±SEM around the mean PSTH line
    params.postStim        double        = 500                            % duration after stimulus onset to include (ms)
    params.preBase         double        = 200                            % pre-stimulus baseline duration (ms)
    params.overwrite       logical       = false                          % if true, recompute PSTHs even when a saved file exists
    params.overwriteResponse logical     = false                          % if true, force recompute of ResponseWindow
    params.overwriteStats  logical       = false                          % if true, force recompute of per-neuron statistics
    params.useCategoryPvals logical      = false                          % if true and splitBy is active, use per-category p-values (OR across levels) instead of general per-neuron p-values
    params.nBootCategory   double        = 10000                          % bootstrap iterations for StatisticsPerNeuronPerCategory
    params.TakeTopPercentTrials double   = 1                              % fraction of trials to keep (1 = all; values <1 bias PSTH amplitudes — see note below)
    params.zScore          logical       = false                          % z-score each neuron's PSTH to its own pre-stimulus baseline
    params.PaperFig        logical       = false                          % export publication-quality figure via printFig
    params.byDepth         logical       = false                          % stratify neurons into 3 cortical depth bins
    params.unionResponsive logical       = true                           % if true, include neurons responsive to ANY stimType (OR-union across stim types)
end

% -------------------------------------------------------------------------
% NOTE ON TakeTopPercentTrials (default = 1 = all trials)
% -------------------------------------------------------------------------
% Selecting the top N% of trials by mean spike count inflates PSTH
% amplitudes and biases the response profile.  For a publication figure
% this parameter should remain at 1 (all trials) unless there is a
% specific, pre-registered reason (e.g. attention gating in a behaving
% animal).
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Guard: splitBy and byDepth together create too many lines to read
% -------------------------------------------------------------------------
if params.splitBy ~= "" && params.byDepth
    error(['splitBy and byDepth cannot both be active — the resulting ', ...
           'combinatorial line count is unreadable.  Use one at a time.']);
end

% -------------------------------------------------------------------------
% Guard: unionResponsive requires ≥2 stim types to be meaningful
% -------------------------------------------------------------------------
if params.unionResponsive && numel(params.stimTypes) < 2
    warning('unionResponsive has no effect with a single stimType — ignoring.');
    params.unionResponsive = false;                                        % disable to avoid misleading output
end

% -------------------------------------------------------------------------
% Guard: unionResponsive + useCategoryPvals is logically ambiguous
% -------------------------------------------------------------------------
if params.unionResponsive && params.useCategoryPvals
    error(['unionResponsive and useCategoryPvals cannot both be true. ', ...
           'The union pre-pass uses general per-neuron p-values (statType). ', ...
           'Use one mode at a time.']);
end

% -------------------------------------------------------------------------
% Guard: requireAllStims only makes sense with ≥2 stim types
% -------------------------------------------------------------------------
if params.requireAllStims && numel(params.stimTypes) < 2
    warning('requireAllStims has no effect with a single stimType — ignoring.');
    params.requireAllStims = false;                                        % disable to avoid unnecessary experiment exclusions
end

% -------------------------------------------------------------------------
% Load depth-bin info if byDepth is requested
% -------------------------------------------------------------------------
if params.byDepth
    % Path to the precomputed depth table produced by getNeuronDepths()
    depthFile = 'W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\NeuronDepths.mat';
    if ~exist(depthFile, 'file')
        error('NeuronDepths.mat not found. Run getNeuronDepths() first.');
    end
    D             = load(depthFile);                                       % load the depth struct from disk
    depthTable    = D.depthTable;                                          % table with columns Experiment, Unit, Depth_um
    depthBinEdges = D.depthBinEdges;                                       % 4-element vector of bin boundaries in µm
    nDepthBins    = 3;                                                     % three cortical depth bins: shallow / middle / deep
    fprintf('Depth bins loaded:\n');
    fprintf('  Bin 1 (shallow): %.0f – %.0f µm\n', depthBinEdges(1), depthBinEdges(2));
    fprintf('  Bin 2 (middle) : %.0f – %.0f µm\n', depthBinEdges(2), depthBinEdges(3));
    fprintf('  Bin 3 (deep)   : %.0f – %.0f µm\n', depthBinEdges(3), depthBinEdges(4));
else
    nDepthBins = 1;                                                        % no depth stratification: treat all neurons as one bin
end

% -------------------------------------------------------------------------
% Build save directory path using the first experiment as a reference
% -------------------------------------------------------------------------
NP_first  = loadNPclassFromTable(exList(end));                               % load NP class for the first experiment to extract paths
vs_first  = linearlyMovingBallAnalysis(NP_first);                         % construct a ball-analysis object just to access getAnalysisFileName

basePath  = extractBefore(vs_first.getAnalysisFileName, 'lizards');        % trim the path at 'lizards' to get the shared root
basePath  = [basePath 'lizards'];                                           % re-append 'lizards' to form the correct base
saveDir   = fullfile(basePath, 'Combined_lizard_analysis');                % combined output directory
if ~exist(saveDir, 'dir')
    mkdir(saveDir);                                                        % create directory if it doesn't already exist
end

 extractBefore(vs_first.getAnalysisFileName, 'lizards');

% ---- Construct the cache filename (encodes key parameter choices) ----
stimLabel   = strjoin(params.stimTypes, '-');                              % e.g. 'RG-MB'
depthSuffix = '';
if params.byDepth;  depthSuffix = '_byDepth'; end                         % append depth tag if stratifying by depth

splitSuffix = '';
if params.splitBy ~= ""; splitSuffix = ['_by' char(params.splitBy)]; end  % append split variable name to filename

if ~isempty(params.splitLevels)
    % Build a string of the requested levels, e.g. '_lvl5_10_20'
    lvlStr = strjoin(arrayfun(@(v) sprintf('%g',v), params.splitLevels, ...
        'UniformOutput', false), '_');
    splitSuffix = [splitSuffix '_lvl' lvlStr];
end

allStimSuffix = '';
if params.requireAllStims; allStimSuffix = '_allStims'; end               % tag the cache file so it is distinct from the unfiltered version

nameOfFile = sprintf('Ex_%d-%d_Combined_PSTHs_%s%s%s%s.mat', ...
    exList(1), exList(end), stimLabel, splitSuffix, depthSuffix, allStimSuffix);
fullSavePath = fullfile(saveDir, nameOfFile);

% -------------------------------------------------------------------------
% Decide: recompute or load from disk
% -------------------------------------------------------------------------
% The cache filename only encodes a few parameters (stim / split / depth /
% allStims).  Many other parameters (binWidth, alpha, statType, preBase,
% postStim, zScore, TakeTopPercentTrials, speed, unionResponsive,
% useCategoryPvals) ALSO change the stored psthAll but are NOT in the
% filename.  We therefore verify them explicitly via computationParamsMatch
% on load, so changing one of them and rerunning forces a recompute instead
% of silently loading a stale result.
% -------------------------------------------------------------------------
forloop = true;                                                            % assume we need to recompute
if exist(fullSavePath, 'file') == 2 && ~params.overwrite
    S = load(fullSavePath);                                                % load the cached struct
    if isequal(S.expList, exList) && isfield(S, 'params') && ...
            computationParamsMatch(S.params, params)
        fprintf('Loading saved PSTHs from:\n  %s\n', fullSavePath);
        forloop = false;                                                   % cached data matches — skip the computation loop
    else
        fprintf('Experiment list or computation params changed — recomputing.\n'); % stale cache; recompute
    end
end

% =========================================================================
%  MAIN COMPUTATION LOOP (skip if loaded from disk)
% =========================================================================
if forloop

    nStim = numel(params.stimTypes);                                       % number of stimulus types requested
    nExp  = numel(exList);                                                 % number of experiments in the list

    % =====================================================================
    %  DISCOVERY PASS — find valid sessions and category levels
    % =====================================================================
    % For each experiment × stim type:
    %   - If splitBy == "" : probe sessions [0, 1, 2] for stimulus PRESENCE
    %     (non-empty VST) and record the session it lives in (or -1).
    %     This is what lets requireAllStims work in the no-split case.
    %   - If splitBy ~= "" : additionally require the category column to
    %     have ≥2 usable levels (or all of splitLevels).
    %
    % Results are stored in sessionMap so the main loop can skip invalid
    % experiments without re-searching.
    %
    %   sessionMap(ei, s)  = session to use (-1 = skip this exp×stim)
    %   catLabelsAll{s}    = string array of category labels for stim s
    % =====================================================================

    sessionMap   = zeros(nExp, nStim);                                     % pre-allocate session map; 0 = default session, -1 = skip
    catLabelsAll = cell(nStim, 1);                                         % will store global category label strings per stim type

    if params.splitBy ~= ""
        fprintf('\nDiscovering categories (splitBy = "%s") ...\n', params.splitBy);
    else
        fprintf('\nDiscovering stimulus presence per experiment ...\n');   % no split: presence/session discovery only
    end

    for s = 1:nStim                                                        % iterate over each stimulus type
        stimKey = params.stimTypes(s);                                     % current stimulus abbreviation (e.g. "MB")

        if params.splitBy == ""
            % ---- No split: still probe each experiment for stimulus
            %      PRESENCE and record the session it lives in (or -1). ----
            % FIX: previously this branch unconditionally set sessionMap to 0
            % for every experiment, so missing stimuli were never detected and
            % requireAllStims had no effect.  It also forced session 0 only.
            catLabelsAll{s} = "all";                                       % single pseudo-category — no splitting

            for ei = 1:nExp
                try
                    NP_tmp = loadNPclassFromTable(exList(ei));             % load NP class for this experiment
                catch
                    sessionMap(ei, s) = -1;                               % experiment could not be loaded — mark skip
                    continue
                end

                % Probe sessions [0, 1, 2]; return first session that has this stim
                sess = findStimSession(NP_tmp, stimKey, params.overwriteResponse);
                sessionMap(ei, s) = sess;                                 % -1 if stim absent in all sessions

                if sess < 0
                    fprintf('  Exp %d [%s]: stimulus not present in any session — will skip.\n', ...
                        exList(ei), stimKey);
                else
                    fprintf('  Exp %d [%s]: stimulus present in session %d.\n', ...
                        exList(ei), stimKey, sess);
                end
            end
        else
            % ---- Split requested: scan each experiment for a valid session ----
            allLevelsFound = [];                                           % accumulate all category levels seen across experiments

            for ei = 1:nExp
                try
                    NP_tmp = loadNPclassFromTable(exList(ei));             % load NP class for this experiment
                catch
                    sessionMap(ei, s) = -1;                               % experiment could not be loaded — mark as skip
                    continue
                end

                % Try sessions [0, 1, 2]; return first with ≥2 valid levels
                [~, sess, levels] = findValidSession( ...
                    NP_tmp, stimKey, params.speed, params.splitBy, ...
                    params.splitLevels, params.overwriteResponse);

                if sess < 0
                    sessionMap(ei, s) = -1;                               % no usable session found for this stim — mark skip
                    fprintf('  Exp %d [%s]: no session with ≥2 levels of "%s" — will skip.\n', ...
                        exList(ei), stimKey, params.splitBy);
                else
                    sessionMap(ei, s) = sess;                             % store the valid session number
                    allLevelsFound = [allLevelsFound; levels(:)];         %#ok<AGROW> accumulate levels seen in this experiment
                    fprintf('  Exp %d [%s]: session %d has levels [%s]\n', ...
                        exList(ei), stimKey, sess, num2str(levels(:)', '%g '));
                end
            end

            % Determine the global set of category levels across all experiments
            uniqueVals = unique(allLevelsFound);                           % sorted unique levels seen across all experiments

            % If user requested specific levels, restrict to the intersection
            if ~isempty(params.splitLevels)
                uniqueVals = intersect(uniqueVals, params.splitLevels(:));
            end

            if numel(uniqueVals) < 2
                % Not enough global levels to split — fall back to unsplit mode.
                % NOTE: in this fall-back we set sessions to 0 (presence was
                % NOT probed for the unsplit case here).  If you intend to
                % combine this fall-back with requireAllStims, prefer calling
                % the function with splitBy == "" directly so presence is probed.
                fprintf('  [%s] splitBy="%s" has only %d global level — falling back to unsplit.\n', ...
                    stimKey, params.splitBy, numel(uniqueVals));
                catLabelsAll{s}  = "all";                                  % treat as unsplit
                sessionMap(:, s) = 0;                                      % reset all to default session
            else
                catLabelsAll{s} = string(uniqueVals(:)');                  % store the final category labels as a string row vector
                fprintf('  [%s] final category levels: %s\n', ...
                    stimKey, strjoin(catLabelsAll{s}, ', '));
            end
        end
    end

    % =====================================================================
    %  EXPERIMENT-LEVEL FILTER: skip experiments missing any stim type
    %  (only applies when requireAllStims = true)
    % =====================================================================
    % By default, an experiment missing stim A but having stim B will still
    % contribute to stim B's PSTH.  When requireAllStims = true, we exclude
    % the entire experiment so that every line in the final plot reflects
    % exactly the same population of experiments.  This is the principled
    % choice for cross-stim comparisons in publication figures.
    %
    % In the split case, an experiment is also excluded if it has the stim
    % but fewer than 2 usable category levels (sessionMap was set to -1
    % above) — i.e. "missing stim" and "missing levels" both trigger
    % exclusion.  State this explicitly in Methods when reporting N.
    % =====================================================================

    if params.requireAllStims
        skipExp = false(nExp, 1);                                          % logical flag: should this experiment be skipped entirely?

        for ei = 1:nExp
            if any(sessionMap(ei, :) < 0)                                 % if any stim type has no valid session for this experiment...
                skipExp(ei) = true;                                        % ...mark the experiment for exclusion
                fprintf('  requireAllStims: Exp %d excluded (missing ≥1 stim type).\n', exList(ei));
            end
        end

        % Force session map to -1 for all stim types in excluded experiments
        for ei = 1:nExp
            if skipExp(ei)
                sessionMap(ei, :) = -1;                                    % mark all stim types as invalid for this experiment
            end
        end

        nKept = sum(~skipExp);                                             % number of experiments that pass the filter
        fprintf('  requireAllStims: %d / %d experiments kept (all stim types present).\n', ...
            nKept, nExp);

        if nKept == 0
            error('requireAllStims: no experiments have all requested stim types [%s].', ...
                strjoin(params.stimTypes, ', '));                           % abort early — nothing to plot
        end
    end

    % ----- Find the maximum number of categories across stim types -------
    maxCats = max(cellfun(@numel, catLabelsAll));                           % largest nCats across all stim types

    % ----- Pre-allocate psthAll: cell array of per-experiment mean traces ----
    % Dimensions: nStim × nDepthBins × maxCats
    % Each cell accumulates one row per CONTRIBUTING experiment (mean firing
    % rate vector).  Skipped experiments simply do not append a row; because
    % the plotting stage reduces each condition independently (it strips
    % all-NaN rows and never relies on cross-stim row alignment), we no longer
    % insert NaN placeholder rows for skipped experiments.
    psthAll = cell(nStim, nDepthBins, maxCats);

    % ----- Provenance: why each experiment did or did not contribute -------
    % An experiment can silently fail to append a row at five different points
    % below, so the stored PSTH alone cannot tell you which sessions the
    % reported N refers to. Record a reason per experiment ("" = contributed)
    % and save it alongside the traces so the N in the figure title is
    % auditable after the fact.
    dropReason = strings(nExp, 1);

    % ----- Time-axis parameters: locked on the first successful experiment ----
    lockedPreBase = [];                                                    % pre-stimulus baseline duration, locked once
    lockedNBins   = [];                                                    % number of time bins, locked once
    lockedEdges   = [];                                                    % bin edge vector in ms, locked once

    % =====================================================================
    %  MAIN EXPERIMENT LOOP
    % =====================================================================
    for ei = 1:nExp

        ex = exList(ei);                                                   % current experiment ID
        fprintf('\n=== Experiment %d (%d/%d) ===\n', ex, ei, nExp);

        % ---- Check if this experiment was excluded by requireAllStims ---
        if params.requireAllStims && all(sessionMap(ei, :) < 0)
            fprintf('  Skipping exp %d (excluded by requireAllStims).\n', ex);
            dropReason(ei) = "requireAllStims: stim type absent";           % record why
            continue                                                       % advance to the next experiment (no row appended)
        end

        % ---- Load NP class for this experiment --------------------------
        try
            NP = loadNPclassFromTable(ex);                                 % load the Neuropixels data class
        catch ME
            warning('Could not load experiment %d: %s', ex, ME.message);
            dropReason(ei) = "experiment failed to load";                   % record why
            continue
        end

        % =================================================================
        %  Union-responsive pre-pass
        %
        %  For each stim type, find responsive neurons (p < alpha) using
        %  per-neuron p-values, then take the OR-union across stim types.
        %  The resulting index vector (unionENeurons) is applied to ALL
        %  stim types in the main loop, so every PSTH line reflects
        %  the same set of neurons within each experiment.
        %
        %  Rationale: spike sorting is stimulus-agnostic — within ONE
        %  recording all stim objects share the same spikeSortingFolder,
        %  so goodU is identical across stim types and neuron indices are
        %  directly comparable.  (If stim objects in a recording could ever
        %  point at different sorting folders, this index reuse would break;
        %  that assumption should be verified for the dataset.)
        %
        %  Note: when requireAllStims = true, this pre-pass only runs on
        %  experiments that passed the completeness filter above, so the
        %  union is always built from ALL requested stim types (no stim
        %  will contribute an empty set due to a missing session).
        % =================================================================
        if params.unionResponsive
            eNeuronsPerStim = cell(1, nStim);                              % one responsive-neuron index vector per stim type

            for su = 1:nStim
                stimTypeU = params.stimTypes(su);                          % stimulus abbreviation for this union pass
                sessU     = sessionMap(ei, su);                            % pre-validated session for this exp × stim

                if sessU < 0                                               % stim type missing or excluded — contributes no neurons
                    eNeuronsPerStim{su} = [];
                    continue
                end

                % ---- Build analysis object for this stim type -----------
                try
                    objU = buildStimObject(NP, stimTypeU, sessU);         % construct stim-specific analysis object
                catch
                    eNeuronsPerStim{su} = [];
                    continue
                end
                if isempty(objU)
                    eNeuronsPerStim{su} = [];
                    continue
                end

                % ---- Check stimulus was actually presented in this session ----
                try
                    stimMissingU = isempty(objU.VST);                     % VST = visual stimulus table; empty means stim not run
                catch
                    stimMissingU = true;
                end
                if stimMissingU
                    eNeuronsPerStim{su} = [];
                    continue
                end

                % ---- Ensure ResponseWindow is computed ------------------
                try
                    objU.ResponseWindow('overwrite', params.overwriteResponse);
                catch
                    eNeuronsPerStim{su} = [];
                    continue
                end

                % ---- Load statistics struct -----------------------------
                try
                    if params.statType == "BootstrapPerNeuron"
                        StatsU = objU.BootstrapPerNeuron;
                    elseif params.statType == "maxPermuteTest"
                        StatsU = objU.StatisticsPerNeuron;                % permutation-test p-values per neuron
                    else
                        StatsU = objU.ShufflingAnalysis;
                    end
                catch
                    eNeuronsPerStim{su} = [];
                    continue
                end

                % ---- Extract p-values using the correct field name ------
                [fieldNameU, ~] = getFieldAndOffset(objU, stimTypeU, params.speed);
                try
                    pvalsU = StatsU.(fieldNameU).pvalsResponse;            % e.g. Stats.Speed1.pvalsResponse
                catch
                    try
                        pvalsU = StatsU.pvalsResponse;                     % flat struct fallback (some older stim types)
                    catch
                        eNeuronsPerStim{su} = [];
                        continue
                    end
                end

                eNeuronsPerStim{su} = find(pvalsU < params.alpha);         % indices of neurons passing the significance threshold
            end

            % ---- Take the OR-union across all stim types ----------------
            unionENeurons = [];                                            % will hold sorted unique neuron indices
            for su = 1:nStim
                unionENeurons = union(unionENeurons, eNeuronsPerStim{su}); % union() returns sorted, unique indices
            end

            fprintf('  Union-responsive: %d neuron(s) responsive to ≥1 of [%s] in exp %d.\n', ...
                numel(unionENeurons), strjoin(params.stimTypes, ', '), ex);

            % ---- If the union is empty, skip this experiment entirely ---
            if isempty(unionENeurons)
                fprintf('  No neurons responsive to any stim type — skipping exp %d.\n', ex);
                dropReason(ei) = "no union-responsive neurons";             % record why
                continue                                                   % advance to the next experiment (no row appended)
            end
        end

        % =================================================================
        %  Per-experiment trace buffer
        % =================================================================
        % Instead of appending each stim's mean trace to psthAll inside the
        % stim loop, we first collect every stim's experiment-mean trace into
        % expBuf and record which (stim, bin, cat) slots produced valid data
        % in validBuf.  The append to psthAll happens ONCE after the stim
        % loop, governed by an inclusion rule (see below).  This is what lets
        % unionResponsive enforce matched experiments across stim lines:
        % a kept experiment contributes to every stim line for a given
        % (bin, cat), or to none — so the per-line N is necessarily equal.
        expBuf   = cell(nStim, nDepthBins, maxCats);                       % per-stim experiment-mean traces (or NaN row)
        validBuf = false(nStim, nDepthBins, maxCats);                      % true where a real (non-all-NaN) trace was produced

        % =================================================================
        %  Loop over stimulus types
        % =================================================================
        for s = 1:nStim

            stimType = params.stimTypes(s);                                % current stimulus abbreviation

            % ---- Check session map: skip if no valid session was found --
            sess = sessionMap(ei, s);
            if sess < 0
                fprintf('  [%s] Skipping exp %d (no valid session).\n', stimType, ex);
                continue
            end

            % ---- Construct the analysis object using the selected session ----
            try
                obj = buildStimObject(NP, stimType, sess);                % build stim-specific analysis object
            catch ME
                warning('Could not build %s (session %d) for exp %d: %s', ...
                    stimType, sess, ex, ME.message);
                continue
            end

            if isempty(obj)
                continue
            end

            % ---- Check that the stimulus was actually presented ---------
            % The constructor may succeed even if this protocol was not run;
            % VST being empty is the reliable indicator.
            try
                stimMissing = isempty(obj.VST);
            catch
                stimMissing = true;
            end
            if stimMissing
                fprintf('  [%s] Stimulus not present in exp %d — skipping.\n', stimType, ex);
                continue
            end

            % ---- Ensure ResponseWindow is computed ----------------------
            try
                obj.ResponseWindow('overwrite', params.overwriteResponse);
            catch ME
                warning('  [%s] ResponseWindow failed for exp %d: %s — skipping.', ...
                    stimType, ex, ME.message);
                continue
            end

            % ---- Load statistics (p-values per neuron) ------------------
            % IMPORTANT: only load per-stim Stats when we will actually use
            % it for neuron selection (i.e. NOT unionResponsive).  When
            % unionResponsive is true the per-stim selection is discarded and
            % replaced by the precomputed union, so loading Stats here is
            % pointless — and worse, a Stats failure for ONE stim would
            % `continue` and drop that stim's row for an experiment the union
            % already kept, making the per-line N unequal (e.g. MB n=16 but
            % SDGm n=14).  Skipping the load keeps every kept experiment
            % contributing to every stim line.
            Stats = [];
            if ~params.unionResponsive
                try
                    if params.statType == "BootstrapPerNeuron"
                        Stats = obj.BootstrapPerNeuron;
                    elseif params.statType == "maxPermuteTest"
                        Stats = obj.StatisticsPerNeuron;                   % preferred: permutation-based p-value per neuron
                    else
                        Stats = obj.ShufflingAnalysis;
                    end
                catch ME
                    warning('  [%s] Statistics failed for exp %d: %s — skipping.', ...
                        stimType, ex, ME.message);
                    continue
                end
            end

            % ---- Determine field name and stim-onset offset -------------
            [fieldName, startStim] = getFieldAndOffset(obj, stimType, params.speed);

            % ---- Get sorted good-unit spike data ------------------------
            p_sort = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder); % load Phy-curated sorting results
            label  = string(p_sort.label');                                % string array of unit labels ('good', 'mua', etc.)
            goodU  = p_sort.ic(:, label == "good");                        % restrict to manually curated 'good' units only

            % ---- Extract condition matrix C and stimulus onset times ----
            C = getConditionMatrix(obj, stimType, params.speed);           % each row = one trial; C(:,1) = onset time in ms
            directimesSorted = C(:, 1)' + startStim;                      % shift onset times by the intra-stimulus offset (e.g. SDGm starts after static phase)

            % ---- Lock the time-axis on the first valid experiment -------
            % All experiments must share the same bin edges for the grand average to be valid.
            preBase     = params.preBase;                                  % pre-stimulus baseline window in ms
            windowTotal = preBase + params.postStim;                       % total window duration per trial in ms

            if isempty(lockedPreBase)
                lockedPreBase = preBase;                                   % store the baseline duration for this analysis run
                lockedEdges   = 0 : params.binWidth : windowTotal;         % bin edges from 0 to windowTotal in steps of binWidth
                lockedNBins   = numel(lockedEdges) - 1;                   % number of time bins
                fprintf('  Locked window: preBase=%d ms, postStim=%d ms, nBins=%d\n', ...
                    lockedPreBase, params.postStim, lockedNBins);
            end

            % ---- Determine whether a category split is active -----------
            nCats         = numel(catLabelsAll{s});                        % number of category levels for this stim type
            isSplitActive = params.splitBy ~= "" && ~isequal(catLabelsAll{s}, "all"); % true = split is in use

            % ---- Select responsive neurons -----------------------------
            % Priority order:
            %   (0) unionResponsive : use the precomputed OR-union across
            %       stim types.  The experiment was already skipped above if
            %       this union was empty, so it is guaranteed non-empty here
            %       and IDENTICAL across all stim types → equal per-line N.
            %       This branch deliberately does not touch per-stim Stats.
            %   (a) useCategoryPvals + active split : OR across per-category
            %       p-values.
            %   (b) Default : overall per-neuron p-values (statType).
            if params.unionResponsive
                % ---- Mode (0): shared union set across stim types -------
                eNeurons = unionENeurons;                                  % same set for every stim type in this experiment
                fprintf('  [%s] Using union-responsive set: %d neuron(s) in exp %d.\n', ...
                    stimType, numel(eNeurons), ex);
            elseif isSplitActive && params.useCategoryPvals
                % ---- Mode (a): per-category p-values, OR across levels --
                try
                    catStats = obj.StatisticsPerNeuronPerCategory( ...
                        'compareCategory', char(params.splitBy), ...
                        'nBoot',           params.nBootCategory, ...
                        'overwrite',       params.overwriteStats);
                catch ME
                    warning('  [%s] StatisticsPerNeuronPerCategory failed for exp %d: %s — skipping.', ...
                        stimType, ex, ME.message);
                    continue
                end

                % Start with a scalar false; first OR will broadcast to a logical vector
                orMask = false;                                            % will expand to [nNeurons × 1] on the first iteration

                for ci = 1:nCats
                    levelVal = str2double(catLabelsAll{s}(ci));            % numeric value of this category level
                    fName    = levelToFieldName(params.splitBy, levelVal); % build field name, e.g. 'size_5'
                    if isfield(catStats, fName)
                        orMask = orMask | (catStats.(fName).pvalsResponse(:) < params.alpha); % OR in the mask for this level
                    end
                end
                eNeurons = find(orMask);                                   % indices of neurons significant for at least one level

                fprintf('  [%s] Using per-category p-values: %d responsive neuron(s) in exp %d.\n', ...
                    stimType, numel(eNeurons), ex);
            else
                % ---- Mode (b): general per-neuron p-values --------------
                try
                    pvals = Stats.(fieldName).pvalsResponse;               % p-values from the main statistics struct
                catch
                    pvals = Stats.pvalsResponse;                           % flat struct fallback
                end
                eNeurons = find(pvals < params.alpha);                     % indices of neurons below the significance threshold
            end

            if isempty(eNeurons)
                fprintf('  [%s] No responsive neurons in exp %d.\n', stimType, ex);
                continue
            end

            % Per-stim count logging (the union branch already printed above)
            if ~params.unionResponsive && ~(isSplitActive && params.useCategoryPvals)
                fprintf('  [%s] %d responsive neuron(s) in exp %d.\n', ...
                    stimType, numel(eNeurons), ex);
            end

            % ---- Extract per-trial category column (only if splitting) --
            catValues = [];
            if isSplitActive
                catCol    = getCategoryColumn(obj, stimType, params.speed, params.splitBy); % extract trial-by-trial category values
                catValues = catCol(:)';                                    % row vector of category value per trial
            end

            % ==============================================================
            %  Per-neuron PSTH computation
            % ==============================================================
            psthRateNeurons = NaN(numel(eNeurons), lockedNBins, nCats);   % [nNeurons × nBins × nCats] — NaN-initialised
            neuronBinIdx    = zeros(numel(eNeurons), 1);                   % depth-bin assignment for each neuron (0 = no depth data)
            baseMeanNeuron  = nan(numel(eNeurons), 1);                     % per-neuron baseline mean (spk/s) for z-scoring
            baseStdNeuron   = nan(numel(eNeurons), 1);                     % per-neuron baseline SD   (spk/s) for z-scoring

            for ni = 1:numel(eNeurons)
                u = eNeurons(ni);                                          % unit index into goodU

                % ---- Per-neuron baseline for z-scoring ------------------
                % Across-trial baseline statistics, pooled over ALL trials
                % (the pre-stimulus period is stimulus-independent, so this is
                % more stable than a per-category estimate and identical across
                % the split levels). Matches sdBase in StatisticsPerNeuron and
                % does NOT shrink with trial count, unlike the SD across PSTH bins.
                if params.zScore
                    MRbase = BuildBurstMatrix( ...
                        goodU(:, u), ...                                   % spike times for this unit
                        round(p_sort.t), ...                               % all sample times (rounded to ms)
                        round(directimesSorted - lockedPreBase), ...       % ALL trial starts (onset - baseline)
                        round(windowTotal));                               % window length in ms
                    MRbase = reshape(MRbase, numel(directimesSorted), []); % [nTrials × time]; reshape avoids singleton collapse
                    perTrialBase = mean(MRbase(:, 1:lockedPreBase), 2) * 1000; % [nTrials × 1] baseline rate (spk/s)
                    baseMeanNeuron(ni) = mean(perTrialBase);               % across-trial baseline mean
                    baseStdNeuron(ni)  = std(perTrialBase);                % across-trial baseline SD
                end

                % ---- Assign neuron to depth bin -------------------------
                if params.byDepth
                    depthRow  = depthTable.Experiment == ex & depthTable.Unit == u; % find this unit's row in the depth table
                    if ~any(depthRow)
                        neuronBinIdx(ni) = 0;                              % unit not found in depth table — will be excluded from averaging
                        continue
                    end
                    unitDepth = depthTable.Depth_um(depthRow);             % cortical depth of this unit in µm
                    unitDepth = unitDepth(1);                              % guard: take the first match if the table has duplicate (Experiment,Unit) rows
                    if unitDepth <= depthBinEdges(2)
                        neuronBinIdx(ni) = 1;                              % shallow bin
                    elseif unitDepth <= depthBinEdges(3)
                        neuronBinIdx(ni) = 2;                              % middle bin
                    else
                        neuronBinIdx(ni) = 3;                              % deep bin
                    end
                else
                    neuronBinIdx(ni) = 1;                                  % all neurons go to the single bin when byDepth is false
                end

                % ---- Build PSTH for each category level -----------------
                for ci = 1:nCats

                    % Build the trial selection mask for this category level
                    if ~isSplitActive
                        trialMask = true(size(directimesSorted));          % all trials belong to the single category
                    else
                        targetVal = str2double(catLabelsAll{s}(ci));       % numeric value of the current category level
                        trialMask = ismembertol(catValues, targetVal, 1e-3); % logical mask: trials matching this level
                    end
                    catOnsets = directimesSorted(trialMask);               % onset times (ms) for the selected trials

                    if isempty(catOnsets)
                        psthRateNeurons(ni, :, ci) = NaN(1, lockedNBins); % no trials for this level — leave as NaN
                        continue
                    end

                    % Build binary spike matrix: rows = trials, columns = ms within window
                    MRhist = BuildBurstMatrix( ...
                        goodU(:, u), ...                                   % spike times for this unit
                        round(p_sort.t), ...                               % all sample times (rounded to ms)
                        round(catOnsets - lockedPreBase), ...              % trial start = onset minus baseline
                        round(windowTotal));                               % window length in ms
                    MRhist = squeeze(MRhist);                              % remove the singleton neuron dimension

                    % FIX: guard the single-trial case.  squeeze() of a
                    % [1 × 1 × W] array yields a [W × 1] column vector, which
                    % would make the per-bin column indexing below address the
                    % wrong dimension.  Force a [nTrials × W] orientation.
                    if numel(catOnsets) == 1
                        MRhist = reshape(MRhist, 1, []);                   % single trial → guarantee a 1 × W row vector
                    end

                    % ---- Optional: keep only top-N% of trials by mean spike count ----
                    % WARNING: this inflates PSTH amplitudes and biases the response
                    % profile — see note above.  Only activate with a principled reason.
                    if ~isempty(params.TakeTopPercentTrials) && params.TakeTopPercentTrials < 1
                        MeanTrial  = mean(MRhist, 2);                      % mean spike count per trial
                        [~, ind]   = sort(MeanTrial, 'descend');           % sort trials from highest to lowest spike count
                        nKeep      = max(1, round(numel(MeanTrial) * params.TakeTopPercentTrials)); % number of trials to retain
                        MRhist     = MRhist(ind(1:nKeep), :);             % keep only the top-N% trials
                    end

                    nTrials = size(MRhist, 1);                             % number of trials (after optional trimming)

                    % ---- Compute PSTH by direct bin summation -----------
                    % Sum spikes within each bin, divide by nTrials and bin
                    % width (converted to seconds) to get spikes/second.
                    counts = zeros(1, lockedNBins);                        % per-bin spike counts, summed across trials
                    for bi = 1:lockedNBins
                        msStart = lockedEdges(bi) + 1;                    % first ms index of this bin (1-based)
                        msEnd   = lockedEdges(bi + 1);                    % last ms index of this bin
                        counts(bi) = sum(MRhist(:, msStart:msEnd), 'all'); % total spikes in this bin across all trials
                    end
                    psthRateNeurons(ni, :, ci) = (counts / nTrials) / (params.binWidth / 1000); % convert to firing rate in spk/s

                end % category loop
            end % neuron loop

            % ---- Report neurons dropped for lacking depth data ----------
            if params.byDepth
                nNoDepth = sum(neuronBinIdx == 0);
                if nNoDepth > 0
                    fprintf('  [%s] %d neuron(s) excluded in exp %d (not in depth table).\n', ...
                        stimType, nNoDepth, ex);
                end
            end

          % ==============================================================
            %  Optional: z-score each neuron's PSTH to its own baseline
            %  Uses the across-trial baseline statistics computed per neuron
            %  in the loop above (one baseline per neuron, shared across all
            %  category levels) — comparable across stimuli with different
            %  trial counts, unlike the SD across PSTH bins.
            % ==============================================================
            if params.zScore
                for ni = 1:size(psthRateNeurons, 1)
                    bMean = baseMeanNeuron(ni);                            % per-neuron baseline mean (spk/s)
                    bStd  = baseStdNeuron(ni);                             % per-neuron baseline SD   (spk/s)
                    for ci = 1:nCats
                        trace = psthRateNeurons(ni, :, ci);               % firing-rate vector for this neuron × category
                        if all(isnan(trace)); continue; end               % skip neuron×cat with no data
                        if bStd > 0
                            psthRateNeurons(ni, :, ci) = (trace - bMean) / bStd; % z-score to across-trial baseline
                        else
                            psthRateNeurons(ni, :, ci) = NaN;             % zero baseline variance — z-score undefined; exclude
                        end
                    end
                end
            end

            % ==============================================================
            %  Average across neurons per depth bin × category, then store
            %  the experiment-level mean trace in the per-experiment buffer.
            %  (The actual append to psthAll happens after the stim loop, so
            %  that inclusion can be made all-or-nothing across stim lines.)
            % ==============================================================
            for b = 1:nDepthBins
                binNeurons = neuronBinIdx == b;                            % logical mask: neurons assigned to depth bin b

                if ~any(binNeurons)
                    % No neurons in this depth bin for this experiment —
                    % leave the buffer slot as a NaN row (validBuf stays false).
                    for ci = 1:nCats
                        expBuf{s, b, ci} = NaN(1, lockedNBins);
                    end
                    continue
                end

                for ci = 1:nCats
                    catData = psthRateNeurons(binNeurons, :, ci);          % [nBinNeurons × nBins] firing rates for this category

                    if all(isnan(catData), 'all')
                        % All neurons NaN for this category — record a NaN row;
                        % validBuf stays false so this slot counts as "no data".
                        expBuf{s, b, ci} = NaN(1, lockedNBins);
                    else
                        psthExp          = mean(catData, 1, 'omitnan');    % mean across neurons (grand-mean for this experiment)
                        expBuf{s, b, ci} = psthExp(:)';                    % store as a row vector
                        validBuf(s, b, ci) = true;                         % mark this slot as having real data
                    end
                end

                fprintf('    [%s] Depth bin %d: %d neuron(s) in exp %d.\n', ...
                    stimType, b, sum(binNeurons), ex);
            end

        end % stim-type loop

        % =================================================================
        %  Commit this experiment's traces to psthAll
        % =================================================================
        % If the time axis was never locked, no stim produced any data for
        % this experiment — nothing to commit.
        if isempty(lockedNBins)
            fprintf('  Exp %d produced no data for any stim — not committed.\n', ex);
            dropReason(ei) = "no usable trace for any stim";                % record why
            continue
        end

        if params.unionResponsive
            % ---- Matched inclusion: per (bin, cat), commit a row to EVERY
            %      stim line only if ALL stims produced a valid trace for
            %      that slot.  This guarantees equal per-line N: every line
            %      reflects exactly the same set of experiments.  If any stim
            %      is missing data for a slot, the experiment is dropped from
            %      ALL stim lines for that slot (and a warning names it). ----
            for b = 1:nDepthBins
                for ci = 1:maxCats
                    % Consider only stim types that actually define this cat index
                    stimsThisSlot = find(arrayfun(@(s) ci <= numel(catLabelsAll{s}), 1:nStim));
                    if isempty(stimsThisSlot); continue; end

                    allValid = all(validBuf(stimsThisSlot, b, ci));        % every relevant stim produced data?

                    if allValid
                        for s = stimsThisSlot
                            psthAll{s, b, ci} = appendOrInit(psthAll{s, b, ci}, expBuf{s, b, ci});
                        end
                    else
                        missing = stimsThisSlot(~validBuf(stimsThisSlot, b, ci));
                        warning(['  Exp %d (bin %d, cat %d): stim(s) [%s] produced no ', ...
                                 'usable trace; dropping this experiment from ALL stim ', ...
                                 'lines for this condition to keep N matched.'], ...
                                 ex, b, ci, strjoin(cellstr(params.stimTypes(missing)), ', '));
                        % Record why, without overwriting an earlier reason
                        if dropReason(ei) == ""
                            dropReason(ei) = "matched-inclusion: no trace for " + ...
                                strjoin(params.stimTypes(missing), '/');
                        end
                    end
                end
            end
        else
            % ---- Independent inclusion: each stim line includes this
            %      experiment iff that stim produced a valid trace.  Lines
            %      may therefore have different N (expected without union). ----
            for s = 1:nStim
                for b = 1:nDepthBins
                    for ci = 1:numel(catLabelsAll{s})
                        if validBuf(s, b, ci)
                            psthAll{s, b, ci} = appendOrInit(psthAll{s, b, ci}, expBuf{s, b, ci});
                        end
                    end
                end
            end
        end

    end % experiment loop

    % =====================================================================
    %  Save results to disk
    % =====================================================================
    % =====================================================================
    %  Experiment-inclusion summary
    % =====================================================================
    % The figure title reports N contributing experiments, but five separate
    % code paths above can drop one. Print the accounting and store it so the
    % reported N can be justified later without re-running.
    contributed = (dropReason == "");                                      % experiments that appended a row
    S.expUsed    = exList(contributed);                                    % experiment IDs behind the reported N
    S.expDropped = exList(~contributed);                                   % experiment IDs excluded
    S.dropReason = dropReason;                                             % per-experiment reason ("" = contributed)

    fprintf('\n===== EXPERIMENT INCLUSION SUMMARY =====\n');
    fprintf('  %d of %d experiments contributed to the PSTH.\n', ...
        sum(contributed), nExp);
    if any(~contributed)
        fprintf('  Excluded:\n');
        reasons = unique(dropReason(~contributed));                        % group the excluded by reason
        for r = 1:numel(reasons)
            hit = dropReason == reasons(r);                                % experiments sharing this reason
            fprintf('    [%s] n=%d: exp %s\n', reasons(r), sum(hit), ...
                strjoin(string(exList(hit)), ', '));
        end
    end
    fprintf('========================================\n');

    S.expList        = exList;                                             % store experiment list for cache-validity check on reload
    S.lockedEdges    = lockedEdges;                                        % bin edge vector, needed to reconstruct tAxis on reload
    S.lockedPreBase  = lockedPreBase;                                      % baseline duration, needed for tAxisPlot alignment
    S.params         = params;                                             % full parameter struct for provenance / reproducibility and cache validation

    % Store catLabelsAll — needed to reconstruct category loop on reload
    S.catLabelsAll   = catLabelsAll;

    % Flatten psthAll cell array into named fields for safe .mat storage
    for s = 1:numel(params.stimTypes)
        stimField = matlab.lang.makeValidName(params.stimTypes(s));        % make the stim abbreviation a valid struct field name
        for b = 1:nDepthBins
            for ci = 1:numel(catLabelsAll{s})
                fieldKey = sprintf('%s_bin%d_cat%d', stimField, b, ci);   % unique field per stim × bin × category
                S.(fieldKey) = psthAll{s, b, ci};                         % [nExp × nBins] matrix for this condition
            end
        end
    end

    save(fullSavePath, '-struct', 'S');                                    % save the struct fields as top-level variables in the .mat file
    fprintf('\nSaved PSTHs to:\n  %s\n', fullSavePath);

else
    % =================================================================
    %  Reload psthAll from the saved struct S (cache hit)
    % =================================================================
    lockedEdges   = S.lockedEdges;                                         % recover bin edges
    lockedPreBase = S.lockedPreBase;                                       % recover baseline duration
    catLabelsAll  = S.catLabelsAll;                                        % recover category labels

    maxCats = max(cellfun(@numel, catLabelsAll));                           % maximum number of categories across stim types
    psthAll = cell(numel(params.stimTypes), nDepthBins, maxCats);          % re-allocate cell array

    for s = 1:numel(params.stimTypes)
        stimField = matlab.lang.makeValidName(params.stimTypes(s));        % re-derive the valid field name
        for b = 1:nDepthBins
            for ci = 1:numel(catLabelsAll{s})
                fieldKey = sprintf('%s_bin%d_cat%d', stimField, b, ci);   % re-construct the field key
                if isfield(S, fieldKey)
                    psthAll{s, b, ci} = S.(fieldKey);                     % load the [nExp × nBins] matrix
                else
                    warning('Field "%s" not found in saved file.', fieldKey);
                    psthAll{s, b, ci} = [];                                % leave empty if field is missing
                end
            end
        end
    end
end

% =========================================================================
%  PLOTTING
% =========================================================================

tAxis     = lockedEdges(1:end-1);                                          % left bin edges within the analysis window (ms)
tAxisPlot = tAxis - lockedPreBase;                                         % shift so time 0 = stimulus onset

% ---- Colour palette and legend label maps --------------------------------
nStim      = numel(params.stimTypes);                                      % number of stimulus types (re-read from params in case forloop was skipped)
baseColors = lines(nStim);                                                 % default MATLAB colour cycle, one colour per stim type

% Map stimulus abbreviations to readable legend labels
stimLegendMap = containers.Map( ...
    {'MB', 'MBR', 'RG', 'SDGm', 'SDGs', 'NV', 'NI', 'FFF'}, ...
    {'MB', 'MBR', 'RG', 'SDGm', 'SDGs', 'NV', 'NI', 'FFF'});

depthShades = [0.05, 0.45, 0.78];                                         % brightness multipliers for depth bins (shallow = brightest)
binLabels   = {'shallow', 'middle', 'deep'};                              % human-readable depth-bin labels

% ---- First pass: smooth traces, compute mean & SEM, find global ylim ---
yMax = -Inf;                                                               % running maximum across all plotted conditions
yMin =  Inf;                                                               % running minimum across all plotted conditions

maxCatsPlot = max(cellfun(@numel, catLabelsAll));                          % re-derive for plotting (in case forloop was skipped)
meanStore   = cell(nStim, nDepthBins, maxCatsPlot);                        % mean PSTH traces keyed by [stim × bin × cat]
semStore    = cell(nStim, nDepthBins, maxCatsPlot);                        % SEM traces keyed by [stim × bin × cat]
nExpStore   = zeros(nStim, nDepthBins, maxCatsPlot);                       % number of valid experiments per condition

for s = 1:nStim
    for b = 1:nDepthBins
        for ci = 1:numel(catLabelsAll{s})
            data = psthAll{s, b, ci};                                      % [nExp × nBins] matrix for this condition
            if isempty(data); continue; end

            validRows = ~all(isnan(data), 2);                              % exclude experiments that contributed all-NaN rows
            data      = data(validRows, :);
            if isempty(data); continue; end

            nValid = size(data, 1);                                        % number of experiments with real data for this condition
            nExpStore(s, b, ci) = nValid;

            % Smooth each experiment's trace individually BEFORE computing
            % the mean and SEM.  This preserves trial-to-trial variability
            % in the SEM estimate while applying the same smoothing uniformly.
            if params.smooth > 0
                smoothBins = max(1, round(params.smooth / params.binWidth)); % convert ms→bins; guard against rounding to 0
                for ri = 1:nValid
                    data(ri, :) = smoothdata(data(ri, :), 'gaussian', smoothBins); % Gaussian-weighted kernel smoothing
                end
            end

            meanTrace = mean(data, 1, 'omitnan');                          % grand mean across experiments
            semTrace  = std(data, 0, 1, 'omitnan') / sqrt(nValid);        % SEM across experiments

            meanStore{s, b, ci} = meanTrace;                              % store for plotting
            semStore{s, b, ci}  = semTrace;

            yMax = max(yMax, max(meanTrace + semTrace));                   % update global y-axis maximum
            yMin = min(yMin, min(meanTrace - semTrace));                   % update global y-axis minimum
        end
    end
end

% Guard: if no condition produced finite data, fall back to a default range
% so that ylim() does not error on [Inf, -Inf].
if ~isfinite(yMax) || ~isfinite(yMin)
    warning('No finite PSTH data to plot — using default y-limits.');
    yMin = 0; yMax = 1;
end

% Add 10% padding to the y-axis range
yPad = (yMax - yMin) * 0.1;
if yPad == 0; yPad = max(abs(yMax), 1) * 0.1; end                         % guard against a zero-width range (single flat trace)
if params.zScore
    yLims = [yMin - yPad, yMax + yPad];                                   % z-scored data can be negative — allow full range
else
    yLims = [max(0, yMin - yPad), yMax + yPad];                           % firing rates cannot be negative — clamp lower bound at 0
end

% ---- Create figure -------------------------------------------------------
fig = figure;
set(fig, 'Units', 'centimeters', 'Position', [5 5 10 7]);                 % set figure size in cm for reproducible layout
ax  = axes(fig);                                                           % create a single axes in the figure
hold(ax, 'on');                                                            % allow multiple plot calls without clearing

legendHandles = [];                                                        % accumulate plot handles for the legend
legendLabels  = {};                                                        % accumulate matching label strings

% ---- Build per-stim colour maps for category splits --------------------
catColorMaps = cell(nStim, 1);                                            % one [nCats × 3] colour matrix per stim type

% Perceptually distinct base colours (blue, orange, purple)
basePalette = [
    0.0000    0.4470    0.7410   % blue
    0.8500    0.3250    0.0980   % orange
    0.4940    0.1840    0.5560   % purple
    ];

for s = 1:nStim

    nc   = numel(catLabelsAll{s});                                         % number of category levels for this stim type
    cmap = zeros(nc, 3);                                                   % colour matrix for this stim type

    for ci = 1:nc
        baseIdx  = mod(ci-1, size(basePalette,1)) + 1;                    % cycle through base colours
        shadeIdx = floor((ci-1) / size(basePalette,1));                   % increment shade tier after one full cycle
        baseC    = basePalette(baseIdx,:);                                 % selected base colour

        if shadeIdx == 0
            newC = baseC;                                                  % first cycle: full-saturation colour
        else
            newC = baseC + (1-baseC)*0.45;                                % subsequent cycles: lighter tint
        end

        cmap(ci,:) = min(max(newC,0),1);                                  % clamp to [0,1] to avoid out-of-range RGB values
    end

    catColorMaps{s} = cmap;                                               % store colour map for this stim type
end

% ---- Plot each stim × depth bin × category condition --------------------
for s = 1:nStim

    stimKey = char(params.stimTypes(s));                                   % current stimulus abbreviation as char for map lookup
    if isKey(stimLegendMap, stimKey)
        shortName = stimLegendMap(stimKey);                                % readable label from the map
    else
        shortName = stimKey;                                               % fallback: use the abbreviation directly
    end

    for b = 1:nDepthBins
        for ci = 1:numel(catLabelsAll{s})

            meanPSTH = meanStore{s, b, ci};                               % grand-mean trace for this condition
            semPSTH  = semStore{s, b, ci};
            if isempty(meanPSTH); continue; end                           % skip conditions with no valid data

            nValid = nExpStore(s, b, ci);                                 % number of experiments contributing to this condition

            % ---- Choose line colour and legend text --------------------
            isSplitHere = params.splitBy ~= "" && ~isequal(catLabelsAll{s}, "all"); % true = this stim type is being split

            if params.byDepth
                % Depth-stratified: darken the base colour by depth bin
                lineColor   = baseColors(s,:) * (1 - depthShades(b));
                legendLabel = sprintf('%s %s (%.0f–%.0f µm, n=%d)', ...
                    shortName, binLabels{b}, ...
                    depthBinEdges(b), depthBinEdges(b+1), nValid);
            elseif isSplitHere
                % Category split: each level gets a distinct colour from catColorMaps
                lineColor   = catColorMaps{s}(ci, :);
                legendLabel = sprintf('%s (n=%d)', catLabelsAll{s}(ci), nValid);
            else
                % Default: one colour per stim type
                lineColor   = baseColors(s,:);
                legendLabel = sprintf('%s (n=%d)', shortName, nValid);
            end

            % ---- SEM shading -------------------------------------------
            if params.shadeSTD && nValid > 1
                upper = meanPSTH + semPSTH;                               % upper bound of the shaded region
                lower = meanPSTH - semPSTH;                               % lower bound
                xFill = [tAxisPlot(:)', fliplr(tAxisPlot(:)')];           % x coordinates: forward then backward for closed polygon
                yFill = [upper(:)',      fliplr(lower(:)')];              % y coordinates: upper then lower (reversed)
                fill(ax, xFill, yFill, lineColor, ...
                    'FaceAlpha', 0.08, 'EdgeColor', 'none');              % semi-transparent filled polygon, no border
            end

            % ---- Mean PSTH line ----------------------------------------
            h = plot(ax, tAxisPlot(:)', meanPSTH(:)', ...
                'Color', lineColor, 'LineWidth', 1.5);                    % plot mean as a solid line

            legendHandles(end+1) = h;                                     %#ok<AGROW> collect handle for legend
            legendLabels{end+1}  = legendLabel;                           %#ok<AGROW> collect label for legend

        end % category loop
    end % depth-bin loop
end % stim-type loop

% ---- Reference lines at stimulus onset and offset -----------------------
xline(ax, 0,               'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off'); % vertical dashed line at t = 0 (stim onset)
xline(ax, params.postStim, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off'); % vertical dashed line at stim offset

% ---- Axis labels and formatting -----------------------------------------
if params.zScore
    yLabel = 'Z-score';                                                    % z-scored mode: dimensionless
else
    yLabel = 'Firing rate [spk/s]';                                       % default: spikes per second
end

xlabel(ax, 'Time re. stim onset [ms]', 'FontName', 'Helvetica', 'FontSize', 8);
ylabel(ax, yLabel,                      'FontName', 'Helvetica', 'FontSize', 8);
xlim(ax, [tAxisPlot(1), tAxisPlot(end)]);                                 % set x range to match the locked window
ylim(ax, yLims);                                                          % apply the globally computed y limits

legend(legendHandles, legendLabels, ...
    'Location', 'best', ...
    'FontName', 'Helvetica', ...
    'FontSize', 7);

ax.FontName       = 'Helvetica';                                           % set all axis text to Helvetica
ax.FontSize       = 8;
ax.YAxis.FontSize = 8;
ax.XAxis.FontSize = 8;
hold(ax, 'off');

% ---- Title: report the ACTUAL number of contributing experiments --------
% nExpStore holds the true per-condition experiment counts; the maximum over
% conditions is the number of experiments that contributed at least one line.
% Reporting numel(exList) here would overstate N after skips / filtering.
nContrib = max(nExpStore(:));
if params.requireAllStims
    titleStr = sprintf('N = %d experiments (all stim types present)', nContrib);
else
    titleStr = sprintf('N = %d experiments', nContrib);
end
title(ax, titleStr, 'FontName', 'Helvetica', 'FontSize', 10);

% ---- Export publication figure if requested -----------------------------
if params.PaperFig
    stimStr = strjoin(params.stimTypes, '-');                              % join stim abbreviations for the filename
    vs_first.printFig(fig, sprintf('PSTH-%s%s%s%s', stimStr, splitSuffix, depthSuffix, allStimSuffix), ...
        PaperFig = params.PaperFig);
end

end  % end of main function


% #########################################################################
%                       LOCAL HELPER FUNCTIONS
% #########################################################################


function validSession = findStimSession(NP, stimKey, overwriteRW)
% findStimSession  Probe sessions [0, 1, 2] and return the first session in
%   which the stimulus was actually presented (non-empty VST) and whose
%   ResponseWindow is computable.  Returns -1 if the stimulus is absent from
%   all sessions.
%
%   Used for the no-split case (splitBy == ""), where category levels are
%   irrelevant but stimulus PRESENCE must still be verified — this is what
%   allows requireAllStims to exclude experiments missing a stim type, and
%   ensures the correct session is selected rather than always assuming 0.

validSession = -1;                                                         % default: stimulus not found in any session

for sess = [0, 1, 2]
    candidate = buildStimObject(NP, stimKey, sess);                       % attempt to build analysis object for this session
    if isempty(candidate)
        continue                                                           % construction failed — try next session
    end

    % Was the stimulus actually presented in this session?
    try
        if isempty(candidate.VST)
            continue                                                       % VST empty → stim not run in this session
        end
    catch
        continue
    end

    % Must be computable downstream (the main loop relies on this)
    try
        candidate.ResponseWindow('overwrite', overwriteRW);
    catch
        continue                                                           % ResponseWindow failed — try next session
    end

    validSession = sess;                                                   % first valid session found
    return
end
end


function [obj, validSession, levels] = findValidSession(NP, stimKey, speedParam, splitBy, splitLevels, overwriteRW)
% findValidSession  Try sessions [0, 1, 2] for a stimulus and return the
%   first session whose category column (splitBy) has ≥2 usable levels.
%   If splitLevels is non-empty, ALL requested levels must be present.
%
%   INPUTS
%       NP           — loaded NP class for this experiment
%       stimKey      — stimulus abbreviation (e.g. "MB", "RG")
%       speedParam   — "max" or other speed selector
%       splitBy      — category column name (e.g. "size")
%       splitLevels  — specific levels required (numeric vector, or [])
%       overwriteRW  — logical: force recompute of ResponseWindow
%
%   OUTPUTS
%       obj          — analysis object for the valid session (or [])
%       validSession — session number (0, 1, or 2), or -1 if none found
%       levels       — unique category levels found in the chosen session

obj          = [];                                                         % default: no valid object found
validSession = -1;                                                         % default: no valid session found
levels       = [];                                                         % default: no levels found

for sess = [0, 1, 2]
    candidate = buildStimObject(NP, stimKey, sess);                       % attempt to build analysis object for this session
    if isempty(candidate)
        continue                                                           % object construction failed — try next session
    end

    % Check that the stimulus was actually presented in this session
    try
        if isempty(candidate.VST)
            continue                                                       % VST empty → stim not run
        end
    catch
        continue
    end

    % Ensure ResponseWindow is computed (required to access C and colNames)
    try
        candidate.ResponseWindow('overwrite', overwriteRW);
    catch
        continue                                                           % ResponseWindow computation failed — try next session
    end

    % Extract condition matrix and column names from ResponseWindow
    rw = candidate.ResponseWindow;
    [C, colNames] = getCmatrixLocal(rw, stimKey, speedParam);             % C: trial × parameter matrix; colNames: parameter names only
    if isempty(C) || isempty(colNames)
        continue                                                           % condition matrix unavailable — try next session
    end

    % Locate the requested category column by name
    catIdx = find(strcmpi(colNames, splitBy), 1);                         % case-insensitive match against parameter names
    if isempty(catIdx)
        continue                                                           % splitBy column not present in this stim type
    end

    % Column layout: colNames(k) → C(:, k+1)  because C(:,1) = onset times
    catColIdx   = catIdx + 1;                                             % offset by 1 to account for the onset column at C(:,1)
    rawCol      = C(:, catColIdx);                                        % raw category values for all trials
    rawCol      = rawCol(~isnan(rawCol));                                  % remove NaN rows (incomplete trials)
    availLevels = uniquetol(rawCol, 1e-6);                                % unique levels with floating-point tolerance

    % If specific levels were requested, verify that ALL are present
    if ~isempty(splitLevels)
        allPresent = true;
        for lv = splitLevels(:)'
            if ~any(abs(availLevels - lv) < 1e-6)
                allPresent = false;                                        % at least one requested level is absent
                break
            end
        end
        if ~allPresent
            continue                                                       % required levels missing — try next session
        end
        useLevels = splitLevels(:);                                       % restrict to the requested subset
    else
        useLevels = availLevels;                                          % use all available levels
    end

    % Need at least 2 levels to make a split meaningful
    if numel(useLevels) < 2
        continue
    end

    % Valid session found — store results and return
    obj          = candidate;
    validSession = sess;
    levels       = useLevels;
    return
end
end


function [C, colNames] = getCmatrixLocal(rw, stimKey, speedParam)
% getCmatrixLocal  Extract condition matrix C and parameter column names
%   from a ResponseWindow struct.
%
%   colNames = stimulus-parameter names only (rw.colNames{1}(5:end)).
%   C(:,1)   = onset times (ms).
%   C(:, k+1) = parameter column for colNames(k).

C        = [];
colNames = {};

% colNames{1}(1:4) are internal bookkeeping columns; parameters start at index 5
try
    allColNames = rw.colNames{1};
    colNames    = allColNames(5:end);                                      % strip the 4 fixed bookkeeping columns
catch
    return                                                                 % colNames unavailable — cannot proceed
end

% Select C from the correct sub-field based on stimulus type
switch stimKey
    case {"MB", "MBR"}
        % Moving-ball/bar stimuli: choose speed field
        if speedParam == "max"
            fld = 'Speed1';                                                % 'Speed1' = fastest speed condition
        else
            fld = 'Speed2';
        end
        if isfield(rw, fld)
            C = rw.(fld).C;
        else
            % Fall back to the last speed field found
            speedFields = fieldnames(rw);
            speedFields = speedFields(startsWith(speedFields, 'Speed'));
            if ~isempty(speedFields)
                C = rw.(speedFields{end}).C;
            end
        end

    case "SDGm"
        % Drifting-grating moving phase
        if isfield(rw, 'C')
            C = rw.C;
        elseif isfield(rw, 'Moving') && isfield(rw.Moving, 'C')
            C = rw.Moving.C;
        end

    case "SDGs"
        % Static grating phase
        if isfield(rw, 'C')
            C = rw.C;
        elseif isfield(rw, 'Static') && isfield(rw.Static, 'C')
            C = rw.Static.C;
        end

    otherwise
        % Generic fallback: use top-level C field
        if isfield(rw, 'C')
            C = rw.C;
        end
end
end


function obj = buildStimObject(NP, stimKey, session)
% buildStimObject  Construct the analysis object for a stimulus key,
%   optionally with a specific session number.
%
%   obj = buildStimObject(NP, "MB", 0)   — default (no Session arg)
%   obj = buildStimObject(NP, "MB", 1)   — Session=1
%   obj = buildStimObject(NP, "MB", 2)   — Session=2
%
%   Returns [] if construction fails or stimKey is unrecognised.

if nargin < 3; session = 0; end                                           % default to session 0 if not specified

obj = [];
try
    % SDGm and SDGs are both served by StaticDriftingGratingAnalysis
    switch stimKey
        case {"SDGm", "SDGs"},  ctorKey = "SDG";                          % map both moving and static phases to the same constructor key
        otherwise,              ctorKey = stimKey;
    end

    if session == 0
        % Default session: call constructor without a Session argument
        switch ctorKey
            case "MB",  obj = linearlyMovingBallAnalysis(NP);
            case "MBR", obj = linearlyMovingBarAnalysis(NP);
            case "RG",  obj = rectGridAnalysis(NP);
            case "SDG", obj = StaticDriftingGratingAnalysis(NP);
            case "NV",  obj = movieAnalysis(NP);
            case "NI",  obj = imageAnalysis(NP);
            case "FFF", obj = fullFieldFlashAnalysis(NP);
            otherwise,  error('Unknown stimulus key: "%s".', stimKey);
        end
    else
        % Non-default session: pass Session as a name-value pair
        switch ctorKey
            case "MB",  obj = linearlyMovingBallAnalysis(NP, 'Session', session);
            case "MBR", obj = linearlyMovingBarAnalysis(NP, 'Session', session);
            case "RG",  obj = rectGridAnalysis(NP, 'Session', session);
            case "SDG", obj = StaticDriftingGratingAnalysis(NP, 'Session', session);
            case "NV",  obj = movieAnalysis(NP, 'Session', session);
            case "NI",  obj = imageAnalysis(NP, 'Session', session);
            case "FFF", obj = fullFieldFlashAnalysis(NP, 'Session', session);
            otherwise,  error('Unknown stimulus key: "%s".', stimKey);
        end
    end
catch ME
    fprintf('  Could not create %s Session=%d: %s\n', stimKey, session, ME.message);
    obj = [];                                                              % return empty on any constructor failure
end
end


function [fieldName, startStim] = getFieldAndOffset(obj, stimKey, speedParam)
% getFieldAndOffset  Return the response-table field name and the
%   intra-stimulus onset offset (ms) for the given stimulus abbreviation.
%
%   For SDGm, startStim accounts for the static phase that precedes the
%   moving phase; the offset shifts trial onsets to the moving-phase start.

startStim = 0;                                                            % default: stimulus onset coincides with the trial onset

switch stimKey
    case "MB"
        if speedParam == "max"; fieldName = 'Speed1'; else; fieldName = 'Speed2'; end
    case "MBR"
        if speedParam == "max"; fieldName = 'Speed1'; else; fieldName = 'Speed2'; end
    case "SDGs"
        fieldName = 'Static';                                             % static-grating phase sub-field
    case "SDGm"
        fieldName = 'Moving';                                             % moving-grating phase sub-field
        startStim = obj.VST.static_time * 1000;                          % static phase duration (s→ms) as onset offset
    otherwise
        fieldName = 'Speed1';                                             % generic fallback — valid for RG, NV, NI, FFF
end
end


function C = getConditionMatrix(obj, stimType, speedParam)
% getConditionMatrix  Extract the condition matrix C from ResponseWindow.
%
%   C(:,1) = trial onset times (ms).
%   C(:,2:end) = stimulus parameters.

[fieldName, ~] = getFieldAndOffset(obj, stimType, speedParam);            % get the correct sub-field name
NeuronResp     = obj.ResponseWindow;                                      % full ResponseWindow struct

try
    C = NeuronResp.(fieldName).C;                                         % retrieve C from the speed/phase sub-field
catch
    C = NeuronResp.C;                                                     % fallback: top-level C (some stim types)
end
end


function catCol = getCategoryColumn(obj, stimType, speedParam, splitBy)
% getCategoryColumn  Extract the per-trial category column from C by
%   matching splitBy against the ResponseWindow parameter column names.
%
%   Column layout (ResponseWindow):
%     colNames{1}(1:4)  = internal bookkeeping (discarded)
%     colNames{1}(5:end) = stimulus-parameter names
%     C(:,1)   = onset times
%     C(:,2:)  = parameter columns
%     → paramNames(k) maps to C(:, k+1)

responseParams = obj.ResponseWindow;                                       % ResponseWindow struct

allColNames = responseParams.colNames{1};                                  % full column name list including bookkeeping headers
paramNames  = allColNames(5:end);                                          % parameter names only (strip the first 4 bookkeeping columns)

matchIdx = find(strcmpi(paramNames, splitBy), 1);                         % case-insensitive search for splitBy column
if isempty(matchIdx)
    error(['splitBy = "%s" does not match any column in colNames.\n' ...
           '  Available: %s'], splitBy, strjoin(string(paramNames), ', '));
end

colIdxInC = matchIdx + 1;                                                  % add 1 because C(:,1) = onset times (paramNames(1) → C(:,2))

C = getConditionMatrix(obj, stimType, speedParam);                        % retrieve the full condition matrix
catCol = C(:, colIdxInC);                                                  % extract the matching column
end


function arr = appendOrInit(arr, newRow)
% appendOrInit  Append newRow to arr, or initialise arr as newRow if empty.
%
%   Used to accumulate per-experiment PSTH rows without pre-allocating
%   the exact number of experiments.
if isempty(arr)
    arr = newRow;                                                          % first experiment: initialise the matrix
else
    arr = [arr; newRow];                                                   % subsequent experiments: vertical concatenation
end
end


function tf = computationParamsMatch(a, b)
% computationParamsMatch  True only if all parameters that affect the STORED
%   psthAll are identical between two parameter structs.  This is the cache-
%   validity check used on reload: changing any of these and rerunning the
%   same exList must force a recompute rather than load a stale result.
%
%   Excluded fields are plotting-only or IO-only and do NOT change psthAll:
%     smooth (applied at plot time), shadeSTD, PaperFig,
%     overwrite, overwriteResponse, overwriteStats.

flds = {'stimTypes', 'requireAllStims', 'splitBy', 'splitLevels', ...
        'binWidth', 'statType', 'speed', 'alpha', 'postStim', 'preBase', ...
        'useCategoryPvals', 'nBootCategory', 'TakeTopPercentTrials', ...
        'zScore', 'byDepth', 'unionResponsive'};

tf = true;
for k = 1:numel(flds)
    f = flds{k};
    if ~isfield(a, f) || ~isfield(b, f) || ~isequal(a.(f), b.(f))
        tf = false;                                                        % any mismatch (or missing field) → cache is stale
        return
    end
end
end


function fName = levelToFieldName(catName, value)
% levelToFieldName  Build a valid struct field name matching the convention
%   used by StatisticsPerNeuronPerCategory.
%
%   Examples:
%     levelToFieldName("size", 5)     →  'size_5'
%     levelToFieldName("speed", 0.3)  →  'speed_0p3'
%     levelToFieldName("dir", -90)    →  'dir_neg90'

fName = sprintf('%s_%g', lower(strtrim(char(catName))), value);           % base name, e.g. 'size_5' or 'speed_0.3'
fName = strrep(fName, '.', 'p');                                          % replace decimal point with 'p' (invalid in field names)
fName = strrep(fName, '-', 'neg');                                        % replace minus sign with 'neg' (invalid in field names)
end