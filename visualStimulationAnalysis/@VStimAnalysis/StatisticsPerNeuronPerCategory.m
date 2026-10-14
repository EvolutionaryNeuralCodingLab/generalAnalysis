function results = StatisticsPerNeuronPerCategory(obj, params)
% StatisticsPerNeuronPerCategory - Per-category statistical analysis of
% neuronal responses.
%
% For a specified stimulus category (e.g. 'size', 'direction', 'speed',
% 'luminosity'), this function:
%
%   1. Tests responsiveness separately for each category level using a
%      sign-flip permutation test (H0: mean response = baseline).
%
%   2. Tests whether responses differ ACROSS levels using a permutation-based
%      one-way ANOVA F-test (omnibus test).
%
%   3. Performs pairwise comparisons between all level pairs using a
%      two-sample permutation test, with FDR correction across all pairs
%      and neurons.
%
% When `CategoryMaximized` is supplied, for each neuron and each level of
% the compared category the trials are restricted to the level of the
% maximizing category that produces the largest mean response. The
% resulting `Diff_k` matrices then have a different number of valid trials
% per neuron and are NaN-padded; ALL three statistical tests below have
% been rewritten to be NaN-aware on a per-neuron basis.
%
% Note on the maximizing option: this is intended for between-level
% effect-size comparisons within a single neuron (e.g. "is size 3 more
% effective than size 1 at the best luminosity?"), NOT for unbiased
% claims about overall responsiveness, since the same trials are used
% to select bestLevel and to test against zero.
%
% Output is saved to a separate file: analysisFileName_categoryname.mat
%
% Category names are matched case-insensitively to responseParams.colNames.
% For linearlyMovingBall, comparing 'speed' receives special handling since
% Speed1 and Speed2 are stored in separate struct fields.
%
% Usage:
%   results = obj.StatisticsPerNeuronPerCategory('compareCategory', 'size')
%   results = obj.StatisticsPerNeuronPerCategory('compareCategory', 'size', ...
%               'CategoryMaximized', 'luminosity', 'overwrite', true)
%
% Reference for permutation tests:
%   Nichols & Holmes (2002) Human Brain Mapping 15:1-25
%
% Changelog (current revision):
%   * Fixed NaN propagation in the sign-flip null (matrix multiplication
%     was making the null distribution all-NaN whenever a column had any
%     NaN padding, which silently broke the maximized branch).
%   * Fixed normalization mismatch (observed used `omitmissing`, null
%     divided by full padded n).
%   * Rewrote computeFstat to compute df_between / df_within and group
%     means PER NEURON (previously `nk` was a scalar that collapsed
%     across neurons, driving F to ~0 in the maximized branch).
%   * Fixed baseline window so it always ends `BaselineBuffer` ms before
%     stimulus onset (previously it could overlap the stimulus when the
%     inter-trial interval was short).
%   * Removed dead variable `idxCompare`.
%   * Added a warning when `CategoryMaximized` is set in a speed-compare
%     run (where it is currently not honoured).

% -------------------------------------------------------------------------
% Argument block. arguments (Input) is the modern MATLAB way to declare a
% Name=Value interface with defaults and (optionally) per-argument
% validators. All fields end up under the struct `params`.
% -------------------------------------------------------------------------
arguments (Input)
    obj                                                                        % parent stimulus-analysis object (provides data accessors)
    params.compareCategory       = ''        % category name to compare (case-insensitive)
    params.nBoot                 = 10000     % permutation iterations
    params.BaseRespWindow        = 1000      % ms response window from stimulus onset (also used as max baseline window)
    params.BaselineBuffer        = 200       % ms buffer before stimulus onset for baseline — avoids contamination from off-responses of the preceding stimulus or anticipatory activity
    params.overwrite             = false     % recompute even if a saved file exists
    params.randomSeed            = 42        % fixed seed for reproducibility
    params.ApplyFDR              = false     % Benjamini-Hochberg FDR correction across pairs × neurons
    params.MovingWindowDuration  = 300       % ms sliding window for moving-ball per-trial peak response (response only; never applied to baseline). Only used when stimulus is linearlyMovingBall.
    params.GratingType           = "moving"  % grating phase to analyse when stimulus is StaticDriftingGrating
    params.CategoryMaximized     = ''        % for each neuron, pick the level of this category that maximizes mean response before running per-level statistics
    params.Lock                  = {}        % {name,value,...} pre-filter: lock other factors to fixed values before splitting by compareCategory (mirrors applyCategoryFocus/StatisticsPerNeuron); {} = no lock
end

% -------------------------------------------------------------------------
% Validate input. compareCategory is mandatory.
% -------------------------------------------------------------------------
if isempty(strtrim(params.compareCategory))                                    % treat empty / whitespace as "not specified"
    error('params.compareCategory must be specified (e.g. ''size'', ''direction'').');
end

% -------------------------------------------------------------------------
% Output file: append the (lowercased, trimmed) category name to the base
% analysis filename so each category gets its own .mat.
% -------------------------------------------------------------------------
outputFile = strrep(obj.getAnalysisFileName, '.mat', ...
    ['_' lower(strtrim(params.compareCategory)) '.mat']);                      % e.g. ..._size.mat

% Short-circuit if a saved result exists, the caller did not request
% overwriting, AND the saved computation params match the requested ones.
% BUG FIX: previously this loaded/returned any existing file unconditionally,
% so a stale cache computed under a different Lock (or nBoot, BaseRespWindow,
% etc.) would be silently returned. Mirrors StatisticsPerNeuron's guard.
if isfile(outputFile) && ~params.overwrite
    cached = load(outputFile);                                                 % peek at the saved result
    if isfield(cached, 'params') && computationParamsMatch(cached.params, params)
        if nargout == 1
            fprintf('Loading saved results from file.\n');
            results = cached;                                                 % params match -> return cached
        else
            fprintf('Analysis already exists (use overwrite option to recalculate).\n');
        end
        return
    else
        fprintf('Computation params changed — recomputing.\n');               % stale cache -> fall through
    end
end

% -------------------------------------------------------------------------
% Optional CategoryMaximized validation.
% -------------------------------------------------------------------------
useMaxCategory = ~isempty(strtrim(params.CategoryMaximized));                  % flag: true if maximizing is requested

if useMaxCategory && strcmpi(strtrim(params.CategoryMaximized), ...
                             strtrim(params.compareCategory))
    % Maximizing OVER the category being compared would collapse it away.
    error('CategoryMaximized must be different from compareCategory.');
end

% -------------------------------------------------------------------------
% Reproducibility.
% -------------------------------------------------------------------------
rng(params.randomSeed);                                                        % seed MATLAB's RNG before any randi/randperm calls

% -------------------------------------------------------------------------
% Load spike-sorted somatic units. `convertPhySorting2tIc` returns
% Phy/Kilosort output reshaped into the `tIc` ("times-by-cluster") layout
% used by BuildBurstMatrix.
% -------------------------------------------------------------------------
p        = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);          % p.t = spike times, p.ic = [4 × nUnits] index array, p.label = Phy quality label
label    = string(p.label');                                                   % convert label cell -> string column for easy comparison
goodU    = p.ic(:, label == 'good');                                           % keep only manually curated "good" (somatic) units
nNeurons = size(goodU, 2);                                                     % number of neurons that pass the quality filter

% Defensive: if no good units survived curation, bail out gracefully.
if isempty(goodU)
    warning('%s has no somatic neurons.', obj.dataObj.recordingName);
    results = [];
    return
end

responseParams = obj.ResponseWindow;                                           % pre-computed per-stimulus metadata (trial timing matrix C, stim duration, column names)

% -------------------------------------------------------------------------
% Make sure diode triggers (the photodiode-based stimulus onset times) are
% synced to the ephys clock. If the sync object is missing we (re)build the
% intermediate session-time and diode-trigger files and try again.
% -------------------------------------------------------------------------
try
    obj.getSyncedDiodeTriggers;
catch
    obj.getSessionTime("overwrite", true);
    obj.getDiodeTriggers("extractionMethod", 'digitalTriggerDiode', 'overwrite', true);
    obj.getSyncedDiodeTriggers;
end

% -------------------------------------------------------------------------
% Identify stimulus type and set per-stimulus flags used below.
% -------------------------------------------------------------------------
isMovingBall = isequal(obj.stimName, 'linearlyMovingBall') || ...               % moving-ball and moving-bar share the same processing path
               isequal(obj.stimName, 'linearlyMovingBar');
isGratingMov = isequal(obj.stimName, 'StaticDriftingGrating') && params.GratingType == "moving"; % drifting (moving) grating phase
isGratingStat= isequal(obj.stimName, 'StaticDriftingGrating') && params.GratingType == "static"; % static grating phase (declared for symmetry; not used downstream)
isSpeedComp  = isMovingBall && strcmpi(strtrim(params.compareCategory), 'speed');                 % special case: comparing speeds in moving-ball

% -------------------------------------------------------------------------
% Warn (but do not error) if the caller asked for maximization on a speed
% comparison — the speed-compare branch ignores CategoryMaximized.
% -------------------------------------------------------------------------
if useMaxCategory && isSpeedComp
    warning(['CategoryMaximized="%s" is currently ignored when ' ...
             'compareCategory="speed". Proceeding without maximization.'], ...
             params.CategoryMaximized);
    useMaxCategory = false;
end

% -------------------------------------------------------------------------
% Guard: Lock must not name the category being compared (or maximized) —
% locking it to a single value would collapse the very comparison requested.
% -------------------------------------------------------------------------
if ~isempty(params.Lock)
    lockNames = string(params.Lock(1:2:end));                                  % every other cell is a factor name
    assert(~any(strcmpi(lockNames, strtrim(params.compareCategory))), ...
        'Lock cannot include the category being compared ("%s").', params.compareCategory);
    if useMaxCategory
        assert(~any(strcmpi(lockNames, strtrim(params.CategoryMaximized))), ...
            'Lock cannot include the CategoryMaximized category ("%s").', params.CategoryMaximized);
    end
end

% -------------------------------------------------------------------------
% Get C matrix, trial times, stimulus duration and category column names.
% C is [nTrials × nCols]; C(:,1) is stimulus onset time, C(:,2:end) are the
% category labels, in the order given by colNames{1}(5:end).
% -------------------------------------------------------------------------
if isMovingBall
    nSpeeds  = numel(unique(obj.VST.speed));                                   % how many speed conditions were presented
    % colNames are identical across speed sub-structs; pull from Speed1.
    % `(5:end)` drops the first 4 columns of the colNames vector (the
    % bookkeeping columns: trial idx, onset, offset, etc.) leaving only
    % the actual stimulus parameter names. Needed in BOTH branches below —
    % the speed-compare branch uses it later to Lock-filter each speed's trials.
    colNames = responseParams.colNames{1}(5:end);

    if isSpeedComp
        % Speed comparison: each speed is itself a level, processed below
        % in its own dedicated block.
        if nSpeeds < 2
            fprintf(['Only one speed condition found in %s. ' ...
                'Cannot compare speeds.\n'], obj.stimName);
            results = [];
            return
        end

        nLevels = nSpeeds;
        levels  = (1:nSpeeds)';                                                % integer levels 1..nSpeeds (Speed1, Speed2, ...)
        fprintf('Comparing %d speed conditions for %s.\n', nSpeeds, obj.stimName);
    else
        % C is overwritten with a speed-pooled version inside the response
        % matrix block. Use Speed1 here only so catIdx / cCol can be found.
        C = responseParams.Speed1.C;
    end

elseif isGratingMov
    % Moving (drifting) phase of grating: shift trial onsets by the static
    % pre-period so timing aligns with the drift onset.
    C        = responseParams.C;
    C(:,1)   = C(:,1) + obj.VST.static_time*1000;                              % static_time is in seconds; convert to ms
    colNames = responseParams.colNames{1}(5:end);

else
    % All other stimuli (rectGrid, static grating phase, etc.)
    C        = responseParams.C;
    colNames = responseParams.colNames{1}(5:end);
end

% -------------------------------------------------------------------------
% Lock filter (non-speed-compare branches only): pre-restrict the reference
% C to locked-in trials BEFORE catIdx/levels are derived from it, so the
% existing "nLevels < 2" guard below also protects against a Lock that
% happens to eliminate a whole level. The speed-compare branch applies the
% identical filter per-speed later, right before its response matrices are built.
% -------------------------------------------------------------------------
if ~isempty(params.Lock) && ~isSpeedComp
    lockKeep = applyCategoryFocus(C, colNames, "all", params.Lock);            % keepMask only; grouping unused here
    C        = C(lockKeep, :);                                                 % restrict reference C to locked-in trials
end

% -------------------------------------------------------------------------
% Find the column of `compareCategory` inside C (and of CategoryMaximized
% if requested). Note: colNames{k} corresponds to C(:, k+1), because C(:,1)
% is stimulus onset.
% -------------------------------------------------------------------------
if ~isSpeedComp
    catIdx = find(strcmpi(colNames, strtrim(params.compareCategory)));         % case-insensitive lookup

    if isempty(catIdx)
        fprintf(['Category "%s" not found in this stimulus.\n' ...
                 'Available categories: %s\n'], ...
                 params.compareCategory, strjoin(colNames, ', '));
        results = [];
        return
    end

    % Optional maximizing category.
    if useMaxCategory
        maxIdx = find(strcmpi(colNames, strtrim(params.CategoryMaximized)));

        if isempty(maxIdx)
            fprintf(['CategoryMaximized "%s" not found in this stimulus.\n' ...
                     'Proceeding without maximizing.\n'], ...
                     params.CategoryMaximized);
            useMaxCategory = false;                                            % silently disable
        else
            maxCol    = maxIdx + 1;                                            % +1 because C(:,1) is onset time
            maxLevels = unique(C(:, maxCol));                                  % candidate levels of the maximizing category

            fprintf('Maximizing across "%s" with %d levels.\n', ...
                params.CategoryMaximized, numel(maxLevels));
        end
    end

    cCol    = catIdx + 1;                                                      % column index of compareCategory in C
    levels  = unique(C(:, cCol));                                              % unique level values [nLevels × 1]
    nLevels = numel(levels);

    if nLevels < 2
        fprintf(['Only one level found for category "%s" in %s. ' ...
                 'Nothing to compare.\n'], params.compareCategory, obj.stimName);
        results = [];
        return
    end

    fprintf('Comparing %d levels of "%s": [%s]\n', nLevels, params.compareCategory, ...
        num2str(levels', '%.4g '));
end

% =========================================================================
% Build response and baseline matrices, compute per-level Diff
% =========================================================================
%
% Sub-function: compute baseline window start offset and duration once so
% the same logic is reused everywhere. The baseline window always ENDS
% `BaselineBuffer` ms before stimulus onset and is at most `BaseRespWindow`
% ms long; if the inter-trial gap is too short, it shrinks accordingly.
%
%   onset - offsetMs ---------- onset - BaselineBuffer ---------- onset
%   |<--------- duration --------->|<------ buffer ------>|
%
preStim   = obj.VST.interTrialDelay*1000 - params.BaselineBuffer;              % ms of pre-stim time available for the baseline window itself
baseDur   = min(preStim, params.BaseRespWindow);                               % actual baseline duration in ms
baseStart = baseDur + params.BaselineBuffer;                                   % offset before onset where the baseline window STARTS
if baseDur <= 0
    error(['Negative or zero baseline duration: interTrialDelay=%.3f s, ' ...
           'BaselineBuffer=%d ms. Reduce BaselineBuffer or check timing.'], ...
           obj.VST.interTrialDelay, params.BaselineBuffer);
end

if isSpeedComp
    % ---------------------------------------------------------------------
    % Speed comparison branch: each speed is a level. Build per-level Diff
    % for each speed in turn.
    % ---------------------------------------------------------------------
    allDiff      = cell(nSpeeds, 1);                                           % cell{k} -> [nTrials_k × nNeurons] Diff for level k
    allBaselines = cell(nSpeeds, 1);                                           % cell{k} -> [nTrials_k × nNeurons] baseline for level k

    for s = 1:nSpeeds
        fName_s      = sprintf('Speed%d', s);                                  % e.g. 'Speed1'
        C_s          = responseParams.(fName_s).C;                             % this speed's condition matrix
        trialTimes_s = C_s(:,1)';                                              % row vector of trial onset times for this speed
        stimDur_s    = responseParams.(fName_s).stimDur;                       % stimulus duration for this speed (ms)

        % Lock filter: restrict this speed's trials to the locked-in subset
        % before building its response/baseline matrices. Applied per-speed
        % (rather than once) because each speed has its own C/trial ordering.
        if ~isempty(params.Lock)
            lockKeep_s   = applyCategoryFocus(C_s, colNames, "all", params.Lock);  % keepMask for this speed
            trialTimes_s = trialTimes_s(lockKeep_s);                           % keep only locked-in trials
        end

        % Response: sliding window over the full stimulus duration (peak in
        % time per trial), because the ball can cross the RF at different
        % times for different speeds and positions.
        Mr_s = BuildBurstMatrix(goodU, round(p.t), ...
            round(trialTimes_s), round(stimDur_s));                            % [nTrials × nNeurons × stimDur_s] binned spike matrix

        if ~isempty(params.MovingWindowDuration)
            assert(size(Mr_s,3) >= params.MovingWindowDuration, ...
                'Speed%d stimulus duration (%d ms) shorter than MovingWindowDuration (%d ms).', ...
                s, size(Mr_s,3), params.MovingWindowDuration);

            mrMov_s     = movmean(Mr_s, params.MovingWindowDuration, 3, ...
                                  'Endpoints', 'discard');                     % moving mean along time, discard edge-truncated windows
            responses_s = max(mrMov_s, [], 3);                                 % [nTrials_s × nNeurons] peak window rate (spikes/ms)
        else
            responses_s = mean(Mr_s, 3);                                       % fallback: simple mean across full stimulus duration
        end

        % Baseline: fixed window ending `BaselineBuffer` ms before onset.
        Mb_s = BuildBurstMatrix(goodU, round(p.t), ...
            round(trialTimes_s - baseStart), baseDur);                         % start `baseStart` ms before onset, span `baseDur` ms

        baselines_s     = mean(Mb_s, 3);                                       % [nTrials_s × nNeurons]
        allDiff{s}      = responses_s - baselines_s;                           % response-minus-baseline per trial
        allBaselines{s} = baselines_s;

        fprintf('Speed%d: using %d ms sliding window over %d ms stimulus.\n', ...
            s, params.MovingWindowDuration, round(stimDur_s));
    end

    % Pooled baseline SD across ALL trials and speeds (one number per neuron).
    sdBase = std(vertcat(allBaselines{:}), 0, 1);                              % [1 × nNeurons]

else
    % ---------------------------------------------------------------------
    % Standard comparison: build full Diff once, split by category level.
    % ---------------------------------------------------------------------

    if isMovingBall
        % Moving ball at a non-speed category: pool trials across ALL
        % speeds (each speed has its own stimDur, so we still process per
        % speed and then concatenate the per-speed response/baseline blocks).
        nSpeeds = numel(unique(obj.VST.speed));

        % Lock filter (per speed): compute once here and reuse below, so the
        % pooled C_all and the per-speed response/baseline blocks stay
        % row-aligned to the identical locked-in trial subset.
        lockKeepBySpeed = cell(nSpeeds, 1);

        % Concatenate C matrices from every speed so the pooled trial
        % ordering matches the pooled response matrix.
        C_all = [];
        for s = 1:nSpeeds
            fName_s = sprintf('Speed%d', s);
            C_s     = responseParams.(fName_s).C;                              % this speed's condition matrix
            if ~isempty(params.Lock)
                lockKeepBySpeed{s} = applyCategoryFocus(C_s, colNames, "all", params.Lock);  % keepMask for this speed
                C_s                = C_s(lockKeepBySpeed{s}, :);                % restrict to locked-in trials
            else
                lockKeepBySpeed{s} = true(size(C_s, 1), 1);                    % no lock -> keep everything
            end
            C_all   = [C_all; C_s];                                            %#ok<AGROW> ok-for-grow: nSpeeds is small
        end
        C    = C_all;                                                          % overwrite C with pooled (and locked) version
        cCol = catIdx + 1;                                                     % column index recomputed defensively (same value)

        % Rebuild unique levels from the pooled C.
        levels  = unique(C(:, cCol));
        nLevels = numel(levels);

        % BUG FIX guard: a Lock can eliminate an entire level's trials (the
        % earlier nLevels<2 check ran on unpooled, unlocked Speed1 data only,
        % before Lock even existed) — re-check here against the true pooled,
        % locked level set actually used for the comparison below.
        if nLevels < 2
            fprintf(['Only %d level(s) of "%s" remain in %s after applying Lock. ' ...
                     'Nothing to compare.\n'], nLevels, params.compareCategory, obj.stimName);
            results = [];
            return
        end

        responsesList = cell(nSpeeds, 1);
        baselinesList = cell(nSpeeds, 1);

        for s = 1:nSpeeds
            fName_s      = sprintf('Speed%d', s);
            trialTimes_s = responseParams.(fName_s).C(:,1)';
            trialTimes_s = trialTimes_s(lockKeepBySpeed{s});                   % same locked-in subset as C_all above
            stimDur_s    = responseParams.(fName_s).stimDur;

            % Response: full stimulus duration with sliding window.
            Mr_s = BuildBurstMatrix(goodU, round(p.t), ...
                round(trialTimes_s), round(stimDur_s));

            if ~isempty(params.MovingWindowDuration)
                assert(size(Mr_s,3) >= params.MovingWindowDuration, ...
                    'Speed%d: stimulus (%d ms) shorter than MovingWindowDuration (%d ms).', ...
                    s, size(Mr_s,3), params.MovingWindowDuration);

                mrMov_s          = movmean(Mr_s, params.MovingWindowDuration, 3, ...
                                           'Endpoints', 'discard');
                responsesList{s} = max(mrMov_s, [], 3);                        % [nTrials_s × nNeurons]
            else
                responsesList{s} = mean(Mr_s, 3);
            end

            % Baseline: fixed window — same logic for all speeds.
            Mb_s = BuildBurstMatrix(goodU, round(p.t), ...
                round(trialTimes_s - baseStart), baseDur);
            baselinesList{s} = mean(Mb_s, 3);

            fprintf('Speed%d: %d ms sliding window over %d ms, %d trials pooled.\n', ...
                s, params.MovingWindowDuration, round(stimDur_s), size(Mr_s,1));
        end

        % Concatenate across speeds: [nTotalTrials × nNeurons].
        responsesFull = vertcat(responsesList{:});
        baselinesFull = vertcat(baselinesList{:});
        DiffFull      = responsesFull - baselinesFull;

        % Pooled baseline SD across all trials and speeds.
        sdBase = std(baselinesFull, 0, 1);                                     % [1 × nNeurons]

    else
        % Fixed-window response for everything that isn't a moving ball.
        trialTimes = C(:,1)';

        MrFull        = BuildBurstMatrix(goodU, round(p.t), ...
            round(trialTimes), params.BaseRespWindow);
        responsesFull = mean(MrFull, 3);                                       % [nTrials × nNeurons]

        % Baseline: window ending `BaselineBuffer` ms before onset.
        MbFull        = BuildBurstMatrix(goodU, round(p.t), ...
            round(trialTimes - baseStart), baseDur);
        baselinesFull = mean(MbFull, 3);                                       % [nTrials × nNeurons]

        DiffFull      = responsesFull - baselinesFull;                         % [nTrials × nNeurons]

        % Pooled baseline SD across all trials.
        sdBase = std(baselinesFull, 0, 1);                                     % [1 × nNeurons]
    end

    % ---------------------------------------------------------------------
    % Split DiffFull by category level (shared code path for moving-ball
    % non-speed and all other stimuli). When CategoryMaximized is set, the
    % per-neuron best level of the maximizing category is chosen and the
    % resulting Diff_k is NaN-padded so downstream code sees a regular
    % rectangular array; ALL stats below are NaN-aware.
    % ---------------------------------------------------------------------
    allDiff = cell(nLevels, 1);

    for k = 1:nLevels
        maskCompare = C(:, cCol) == levels(k);                                 % logical mask: trials at this level of compareCategory

        if ~useMaxCategory
            % Plain version: keep every trial at this level.
            allDiff{k} = DiffFull(maskCompare, :);
            fprintf('Level %g: %d trials\n', levels(k), sum(maskCompare));

        else
            % Maximizing version: per neuron, pick the maxCategory level
            % with the largest mean response, then keep only those trials.
            Diff_k = nan(sum(maskCompare), nNeurons);                          % allocate NaN-padded matrix

            for n = 1:nNeurons
                bestResp  = -inf;                                              % best mean Diff seen so far for this neuron
                bestLevel = nan;                                               % which level of maxCategory produced it

                for m = 1:numel(maxLevels)
                    maskTmp = maskCompare & (C(:, maxCol) == maxLevels(m));    % trials at this (compare, max) combination
                    if ~any(maskTmp)
                        continue                                               % no trials in this cell
                    end
                    muTmp = mean(DiffFull(maskTmp, n), 1);                     % mean Diff for neuron n in this cell
                    if muTmp > bestResp
                        bestResp  = muTmp;
                        bestLevel = maxLevels(m);
                    end
                end

                finalMask = maskCompare & (C(:, maxCol) == bestLevel);         % select trials at the best maxCategory level
                tmp = DiffFull(finalMask, n);                                  % Diff values for this neuron at its best cell
                Diff_k(1:numel(tmp), n) = tmp;                                 % top-fill column n; remaining rows stay NaN
            end

            allDiff{k} = Diff_k;
            fprintf('Level %g: maximizing "%s"\n', levels(k), params.CategoryMaximized);
        end
    end
end

% =========================================================================
% Per-level responsiveness test
% Sign-flip permutation: H0: mean(Diff) = 0 for this level.
% NaN-aware so the maximized branch works correctly.
% One-tailed (excitatory): p = proportion of null >= observed.
% Z-score: bias-corrected by subtracting null mean, normalised by sdBase.
% =========================================================================
pValPerLevel  = nan(nLevels, nNeurons);                                        % [nLevels × nNeurons] p-values per level vs baseline
zPerLevel     = nan(nLevels, nNeurons);                                        % [nLevels × nNeurons] z-scores
obsStatLevels = nan(nLevels, nNeurons);                                        % [nLevels × nNeurons] observed mean Diff

for k = 1:nLevels
    Diff_k    = allDiff{k};                                                    % [nTrials_k × nNeurons], possibly NaN-padded
    nTrials_k = size(Diff_k, 1);                                               % padded row count (same across neurons)

    % Per-neuron valid trial counts; zero-fill NaNs so they contribute 0
    % to both the observed sum and the sign-flipped null sum. Sign × 0 = 0,
    % so randomly flipping the sign of a padded position doesn't affect the
    % null statistic.
    isValid    = ~isnan(Diff_k);                                               % [nTrials_k × nNeurons]
    nValid_k   = sum(isValid, 1);                                              % [1 × nNeurons]
    Diff_k0    = Diff_k;
    Diff_k0(~isValid) = 0;                                                     % zero-filled copy used for matrix arithmetic

    % Observed mean per neuron (each neuron divided by its own valid n).
    ObsStat_k          = sum(Diff_k0, 1) ./ max(nValid_k, 1);                  % [1 × nNeurons]
    ObsStat_k(nValid_k == 0) = NaN;                                            % neurons with no valid trials -> NaN (not 0)
    obsStatLevels(k,:) = ObsStat_k;

    % Sign-flip null. signs_k is ±1 with each column an independent draw.
    signs_k    = 2*randi(2, nTrials_k, params.nBoot) - 3;                      % [nTrials_k × nBoot], values in {-1, +1}
    nullSum_k  = signs_k' * Diff_k0;                                           % [nBoot × nNeurons], NaN-free because of zero-fill
    nullDist_k = nullSum_k ./ max(nValid_k, 1);                                % broadcast: divide each column by its own valid n
    nullDist_k(:, nValid_k == 0) = NaN;                                        % keep empty-neuron columns as NaN

    % One-tailed p-value: proportion of null >= observed.
    pValPerLevel(k,:) = mean(nullDist_k >= ObsStat_k, 1, 'omitnan');

    % Bias-corrected z-score.
    nullMean_k                = mean(nullDist_k, 1, 'omitnan');                % [1 × nNeurons]
    z_k                       = (ObsStat_k - nullMean_k) ./ sdBase;            % [1 × nNeurons]
    z_k(sdBase == 0)          = 0;                                             % avoid 0/0 for silent neurons
    z_k(nValid_k == 0)        = NaN;                                           % preserve missingness
    zPerLevel(k,:)            = z_k;
end

% =========================================================================
% Omnibus test: permutation one-way ANOVA F-test
% H0: mean response is equal across all levels.
% Permutes level labels nBoot times. computeFstat (helper at bottom)
% computes the F-statistic per neuron with full NaN-awareness.
% =========================================================================
DiffPooled   = vertcat(allDiff{:});                                            % [nTotalTrials × nNeurons]
nTotalTrials = size(DiffPooled, 1);

% Level label vector: integer k repeated nTrials_k times for each level.
levelLabelsVec = cell2mat(arrayfun(@(k) ...
    k * ones(size(allDiff{k},1), 1), ...
    (1:nLevels)', 'UniformOutput', false));                                    % [nTotalTrials × 1]

% Observed F per neuron.
F_obs = computeFstat(DiffPooled, levelLabelsVec, levels);                      % [1 × nNeurons]

% Null distribution: permute labels and recompute F.
nullF = zeros(params.nBoot, nNeurons);
for b = 1:params.nBoot
    permLabels = levelLabelsVec(randperm(nTotalTrials));                       % shuffle labels across trials
    nullF(b,:) = computeFstat(DiffPooled, permLabels, levels);
end

pValOmnibus = mean(nullF >= F_obs, 1, 'omitnan');                              % [1 × nNeurons]

% =========================================================================
% Pairwise comparisons: two-sample permutation test for each level pair.
% Observed: mean(Diff_i) - mean(Diff_j).
% Null: randomly reassign trials between the two groups, recompute mean
% difference. Two-tailed: |observed| >= |null|.
% NaN-aware: each group's mean is computed per neuron using its own valid n.
% FDR correction across all pairs × neurons.
% =========================================================================
pairIdx = nchoosek(1:nLevels, 2);                                              % [nPairs × 2]
nPairs  = size(pairIdx, 1);

pValPairwise    = nan(nPairs, nNeurons);                                       % raw p-values
obsStatPairwise = nan(nPairs, nNeurons);                                       % observed mean differences

for pr = 1:nPairs
    i      = pairIdx(pr,1);
    j      = pairIdx(pr,2);
    Diff_i = allDiff{i};                                                       % [nTrials_i × nNeurons], possibly NaN-padded
    Diff_j = allDiff{j};
    nI     = size(Diff_i, 1);
    nJ     = size(Diff_j, 1);

    % Observed mean difference per neuron (NaN-aware).
    ObsDiff_ij            = mean(Diff_i, 1, 'omitnan') - mean(Diff_j, 1, 'omitnan');
    obsStatPairwise(pr,:) = ObsDiff_ij;

    % Pool trials and prepare zero-filled copy + per-row validity mask.
    DiffPair        = [Diff_i; Diff_j];                                        % [nI+nJ × nNeurons]
    isValidPair     = ~isnan(DiffPair);                                        % [nI+nJ × nNeurons]
    DiffPair0       = DiffPair;
    DiffPair0(~isValidPair) = 0;
    nTotal_pr       = nI + nJ;

    nullPair = zeros(params.nBoot, nNeurons);
    for b = 1:params.nBoot
        perm = randperm(nTotal_pr);
        idx1 = perm(1:nI);
        idx2 = perm(nI+1:end);

        % Per-neuron group sums and valid counts.
        s1 = sum(DiffPair0(idx1,:),    1);                                     % [1 × nNeurons]
        s2 = sum(DiffPair0(idx2,:),    1);
        n1 = sum(isValidPair(idx1,:),  1);
        n2 = sum(isValidPair(idx2,:),  1);

        m1 = s1 ./ max(n1, 1);
        m2 = s2 ./ max(n2, 1);
        m1(n1 == 0) = NaN;
        m2(n2 == 0) = NaN;

        nullPair(b,:) = m1 - m2;
    end

    % Two-tailed p-value.
    pValPairwise(pr,:) = mean(abs(nullPair) >= abs(ObsDiff_ij), 1, 'omitnan');
end

% FDR correction across all (pair × neuron) cells simultaneously.
if params.ApplyFDR && nPairs > 1
    pFlat     = pValPairwise(:);                                               % [nPairs*nNeurons × 1]
    validMask = ~isnan(pFlat);
    pAdj      = nan(size(pFlat));
    if any(validMask)
        [pAdj(validMask), ~, ~, ~] = fdr_BH(pFlat(validMask), 0.05, false);    % Benjamini-Hochberg, two-stage off
    end
    pValPairwiseAdj = reshape(pAdj, nPairs, nNeurons);                         % [nPairs × nNeurons]
else
    pValPairwiseAdj = pValPairwise;
end

% =========================================================================
% Store results in struct S and save.
% =========================================================================
S.categoryName      = params.compareCategory;                                  % name of compared category
S.categoryLevels    = levels;                                                  % actual level values
S.params            = params;
S.CategoryMaximized = params.CategoryMaximized;

% Per-level results — field named by category and level value.
for k = 1:nLevels
    % Build a valid MATLAB struct field name from category name + level value.
    fName_k = sprintf('%s_%g', lower(strtrim(params.compareCategory)), levels(k));
    fName_k = strrep(fName_k, '.', 'p');                                       % '0.5' -> '0p5' (valid field)
    fName_k = strrep(fName_k, '-', 'neg');                                     % '-1'  -> 'neg1'

    S.(fName_k).pvalsResponse = pValPerLevel(k,:);                             % [1 × nNeurons] p-values vs baseline
    S.(fName_k).ZScoreU       = zPerLevel(k,:);                                % [1 × nNeurons] bias-corrected z-score
    S.(fName_k).ObsStat       = obsStatLevels(k,:);                            % [1 × nNeurons] mean Diff (spikes/ms)
    S.(fName_k).nTrials       = size(allDiff{k}, 1);                           % padded row count for this level
    S.(fName_k).nValid        = sum(~isnan(allDiff{k}), 1);                    % [1 × nNeurons] valid-trial count per neuron
end

% Omnibus test.
S.omnibus.pVal  = pValOmnibus;                                                 % [1 × nNeurons] any difference across levels?
S.omnibus.F_obs = F_obs;                                                       % [1 × nNeurons] observed F-statistic

% Pairwise results.
for pr = 1:nPairs
    i        = pairIdx(pr,1);
    j        = pairIdx(pr,2);
    pairName = sprintf('%s_%g_vs_%g', ...
        lower(strtrim(params.compareCategory)), levels(i), levels(j));
    pairName = strrep(pairName, '.', 'p');
    pairName = strrep(pairName, '-', 'neg');

    S.pairwise.(pairName).pVal    = pValPairwise(pr,:);                        % [1 × nNeurons] raw
    S.pairwise.(pairName).pValAdj = pValPairwiseAdj(pr,:);                     % [1 × nNeurons] FDR-corrected
    S.pairwise.(pairName).obsDiff = obsStatPairwise(pr,:);                     % [1 × nNeurons] level_i minus level_j
end

fprintf('Saving results to %s\n', outputFile);
save(outputFile, '-struct', 'S');
results = S;

end % end main function


% =========================================================================
%% Local helper: cache-staleness guard (mirrors StatisticsPerNeuron's)
% =========================================================================
function tf = computationParamsMatch(pSaved, pNew)
% computationParamsMatch  True iff every computation-relevant param matches.
%   Generic (compares all fields) so new params (e.g. Lock) are covered
%   automatically; only non-computational flags are ignored. isequaln
%   handles NaN and the nested cell in Lock.
    ignore = {'overwrite'};                                       % flag that doesn't affect the result
    fn = union(fieldnames(pSaved), fieldnames(pNew));             % all fields on either struct
    tf = true;                                                    % assume match until proven otherwise
    for i = 1:numel(fn)
        f = fn{i};                                                % current field name
        if any(strcmp(f, ignore)), continue; end                 % skip ignored flags
        if ~isfield(pSaved, f) || ~isfield(pNew, f) || ...        % field missing on one side, or
                ~isequaln(pSaved.(f), pNew.(f))                   % values differ (NaN-safe, recursive)
            tf = false; return;                                   % any mismatch -> stale
        end
    end
end


% =========================================================================
%% Local helper: permutation-compatible one-way ANOVA F-statistic
% =========================================================================
function F = computeFstat(Diff, levelLabels, levels)
% computeFstat - One-way ANOVA F-statistic for permutation testing.
%
% Inputs:
%   Diff        : [nTrials × nNeurons] response-minus-baseline values
%                 (may contain NaN; treated as missing per neuron)
%   levelLabels : [nTrials × 1] integer group membership per trial
%   levels      : unique level values [nLevels × 1] (only its length is used)
%
% Output:
%   F           : [1 × nNeurons] F-statistic per neuron, NaN-safe
%
% The previous version pooled `nk` across neurons (`all(~isnan(...),2)`)
% which collapsed to zero whenever different neurons had different NaN
% patterns (the maximized branch). This rewrite computes nk, group means,
% the grand mean, df_between, and df_within per neuron.

    nLevels  = numel(levels);
    nNeurons = size(Diff, 2);

    % Validity mask and zero-filled copy: zero-filled values contribute
    % nothing to sums (so sums-of-squares stay correct as long as we use
    % the SAME zero-filled copy and the per-neuron group means).
    isValid          = ~isnan(Diff);                                           % [nTotal × nNeurons]
    nValidTotal      = sum(isValid, 1);                                        % [1 × nNeurons]
    Diff0            = Diff;
    Diff0(~isValid)  = 0;

    % Grand mean per neuron (using each neuron's own valid n).
    grandMean        = sum(Diff0, 1) ./ max(nValidTotal, 1);                   % [1 × nNeurons]

    SS_between   = zeros(1, nNeurons);
    SS_within    = zeros(1, nNeurons);
    nGroupsUsed  = zeros(1, nNeurons);                                         % count of non-empty groups per neuron (for df_between)

    for k = 1:nLevels
        mask     = levelLabels == k;                                           % trials assigned to group k under the current label vector
        Dg0      = Diff0(mask, :);                                             % zero-filled group block
        isVg     = isValid(mask, :);                                           % per-neuron validity within this group
        nk_per_n = sum(isVg, 1);                                               % [1 × nNeurons] valid trials in group k per neuron

        % Per-neuron group mean (use 1 in the denominator when nk=0 to
        % avoid division by zero; that group simply contributes nothing
        % because nk_per_n is 0 in the SS_between term).
        groupMean = sum(Dg0, 1) ./ max(nk_per_n, 1);                           % [1 × nNeurons]

        % SS_between = sum_k nk * (groupMean_k - grandMean)^2, per neuron.
        SS_between = SS_between + nk_per_n .* (groupMean - grandMean).^2;

        % SS_within = sum_k sum_i (x_ik - groupMean_k)^2, per neuron.
        % Subtract groupMean within the group, then zero out the contribution
        % of invalid (NaN-padded) rows so they don't add (groupMean)^2 noise.
        dev          = Dg0 - groupMean;                                        % broadcasts to [nk × nNeurons]
        dev(~isVg)   = 0;                                                      % invalid rows contribute 0
        SS_within    = SS_within + sum(dev.^2, 1);

        % Count this group for df purposes only if it contributed any
        % valid trials for that neuron.
        nGroupsUsed  = nGroupsUsed + (nk_per_n > 0);
    end

    df_between = max(nGroupsUsed - 1, 1);                                      % per neuron; floor at 1 to avoid 0/0
    df_within  = max(nValidTotal - nGroupsUsed, 1);                            % per neuron

    F                       = (SS_between ./ df_between) ./ (SS_within ./ df_within);
    F(SS_within == 0)       = 0;                                               % degenerate: all trials identical within groups
    F(nValidTotal == 0)     = NaN;                                             % truly empty neurons -> NaN
end