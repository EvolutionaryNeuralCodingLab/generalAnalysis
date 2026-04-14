function results = StatisticsPerNeuron(obj, params)
% StatisticsPerNeuron - Computes per-neuron response statistics vs baseline.
%
% For each neuron this function outputs:
%   pvalsResponse : p-value from a max-statistic sign-flip permutation test.
%                   Tests H0: no stimulus category drives a response above baseline.
%                   The max-statistic controls family-wise error rate across categories
%                   without requiring Bonferroni correction.
%
%   ZScoreU       : Data-driven z-score of neuronal response normalised by pooled
%                   baseline SD. Three modes controlled by MovingWindow and UseLOO:
%                   - MovingWindow=true : peak 300ms sliding window at preferred
%                     category (argmax of MW), baseline corrected.
%                   - MovingWindow=false, UseLOO=true : LOO cross-validated mean
%                     Diff at preferred category — unbiased across stimuli.
%                   - MovingWindow=false, UseLOO=false : direct mean Diff at
%                     preferred category — faster but subject to winner's curse.
%
%   ZScorePermutation : Permutation z-score — observed max-statistic normalised
%                   by the mean and SD of its own null distribution.
%                   Quantifies how many SDs above the null the observed response is.
%                   More comparable across stimuli than ZScoreU when stimulus
%                   durations or category counts differ substantially.
%                   Note: still partially affected by nCats and duration since
%                   nullSD scales with both. Use alongside ZScoreU, not instead.
%
%   prefCat       : Consensus preferred category index [1 × nNeurons].
%
%   validCats     : [nCats × nNeurons] logical mask. False where a category has
%                   >= EmptyTrialPerc fraction of zero-spike trials.
%
%   pValTTest     : p-value from one-sample t-test against zero, pooled across
%                   all valid categories per neuron.
%
%   tStat         : t-statistic corresponding to pValTTest [1 × nNeurons].
%
% Usage:
%   results = obj.StatisticsPerNeuron()
%   results = obj.StatisticsPerNeuron(nBoot=5000, UseLOO=false, overwrite=true)
%
% Reference for sign-flip permutation test:
%   Nichols & Holmes (2002) Human Brain Mapping 15:1-25

arguments (Input)
    obj
    params.nBoot                = 10000  % number of permutation iterations for null distribution
    params.EmptyTrialPerc       = 0.7    % exclude category if fraction of zero-spike trials >= this threshold
    params.FilterEmptyResponses = false  % whether to apply empty-trial category filtering
    params.overwrite            = false  % if true, recompute even if a saved file already exists
    params.randomSeed           = 42     % fixed seed for reproducible permutation results (required for publication)
    params.MovingWindowPval     = true   % if true: use per-trial sliding window max for pval
                                         % if false: use full stimulus duration mean
    params.UseLOO               = true   % if true: LOO cross-validated z-score (recommended)
                                         % if false: direct z-score at preferred category (faster, inflated)
                                         % ignored when MovingWindow=true (prefCat from argmax of MW)
    params.CapStimDuration      = true   % if true: cap stimulus duration at MaxStimDuration ms
                                         % before building response matrix. Ensures comparable
                                         % analysis windows across stimuli with different durations.
    params.MaxStimDuration      = 500    % maximum stimulus duration in ms when CapStimDuration=true.
                                         % Should be set to the duration of the shortest stimulus
                                         % (e.g. 500ms for rectGrid) for cross-stimulus comparability.
end

% -------------------------------------------------------------------------
% Load cached results if available
% -------------------------------------------------------------------------
if isfile(obj.getAnalysisFileName) && ~params.overwrite
    if nargout == 1
        fprintf('Loading saved results from file.\n');
        results = load(obj.getAnalysisFileName);  % return previously computed results
    else
        fprintf('Analysis already exists (use overwrite option to recalculate).\n');
    end
    return
end

% -------------------------------------------------------------------------
% Fix random seed for reproducibility
% Required for published code so permutation results are identical across runs
% -------------------------------------------------------------------------
rng(params.randomSeed);

% -------------------------------------------------------------------------
% Load spike-sorted units
% -------------------------------------------------------------------------
p     = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);  % load kilosort/phy output
label = string(p.label');                                             % unit quality labels as strings
goodU = p.ic(:, label == 'good');                                     % keep only somatic ('good') units
responseParams = obj.ResponseWindow;                                  % stimulus timing and category structure

% -------------------------------------------------------------------------
% Handle case with no somatic neurons — save empty struct and return
% -------------------------------------------------------------------------
if isempty(goodU)
    warning('%s has no somatic neurons, skipping experiment.\n', obj.dataObj.recordingName);
    S        = buildEmptyStruct(obj, responseParams);  % consistent empty output struct
    S.params = params;
    save(obj.getAnalysisFileName, '-struct', 'S');
    results = S;
    return
end

% -------------------------------------------------------------------------
% Sync diode triggers for stimulus alignment
% Wrapped in try/catch because trigger files may need to be regenerated
% on first run or after recording issues
% -------------------------------------------------------------------------
try
    obj.getSyncedDiodeTriggers;
catch
    obj.getSessionTime("overwrite", true);                                               % regenerate session time file
    obj.getDiodeTriggers("extractionMethod", 'digitalTriggerDiode', 'overwrite', true);  % re-extract diode triggers
    obj.getSyncedDiodeTriggers;                                                          % retry sync
end

% -------------------------------------------------------------------------
% Parse stimulus timing per condition
% Stimulus type determines loop structure:
%   linearlyMovingBall/Bar → one or two speed conditions (Speed1, Speed2)
%   StaticDriftingGrating  → Static and Moving phases
%   all others (rectGrid)  → single condition
% -------------------------------------------------------------------------
if isfield(responseParams, "Speed1")
    % BUG FIX: original code used length(obj.VST.speed) which returns total
    % number of trials — corrected to numel(unique(...)) for distinct speeds
    nSpeeds = numel(unique(obj.VST.speed));

    Times.Speed1      = responseParams.Speed1.C(:,1)';
    Durations.Speed1  = responseParams.Speed1.stimDur;
    trialsCats.Speed1 = round(numel(Times.Speed1) / size(responseParams.Speed1.NeuronVals, 2));
    MWs.Speed1        = responseParams.Speed1.NeuronVals(:,:,4)';  % [nCats × nNeurons]

    if nSpeeds > 1
        Times.Speed2      = responseParams.Speed2.C(:,1)';
        Durations.Speed2  = responseParams.Speed2.stimDur;
        trialsCats.Speed2 = round(numel(Times.Speed2) / size(responseParams.Speed2.NeuronVals, 2));
        MWs.Speed2        = responseParams.Speed2.NeuronVals(:,:,4)';
    end

    x = nSpeeds;

elseif isequal(obj.stimName, 'StaticDriftingGrating')
    Times.Moving      = responseParams.C(:,1)' + obj.VST.static_time * 1000;
    Durations.Moving  = responseParams.Moving.stimDur;
    trialsCats.Moving = round(numel(Times.Moving) / size(responseParams.Moving.NeuronVals, 2));
    MWs.Moving        = responseParams.Moving.NeuronVals(:,:,4)';

    Times.Static      = responseParams.C(:,1)';
    Durations.Static  = responseParams.Static.stimDur;
    trialsCats.Static = round(numel(Times.Static) / size(responseParams.Static.NeuronVals, 2));
    MWs.Static        = responseParams.Static.NeuronVals(:,:,4)';

    FieldNames = {'Static', 'Moving'};
    x = 2;

else
    directimesSorted = responseParams.C(:,1)';
    stimDur          = responseParams.stimDur;
    trialsCat        = round(numel(directimesSorted) / size(responseParams.NeuronVals, 2));
    MW               = responseParams.NeuronVals(:,:,4)';  % [nCats × nNeurons]
    x = 1;
end

% =========================================================================
% Main loop over stimulus conditions
% =========================================================================
for s = 1:x

    % --- Assign condition-specific variables ---
    if isfield(responseParams, "Speed1")
        fieldName        = sprintf('Speed%d', s);
        directimesSorted = Times.(fieldName);
        stimDur          = Durations.(fieldName);
        trialsCat        = trialsCats.(fieldName);
        MW               = MWs.(fieldName);
    end

    if isequal(obj.stimName, 'StaticDriftingGrating')
        fieldName        = FieldNames{s};
        directimesSorted = Times.(fieldName);
        stimDur          = Durations.(fieldName);
        trialsCat        = trialsCats.(fieldName);
        MW               = MWs.(fieldName);
    end

    % -------------------------------------------------------------------------
    % Cap stimulus duration if requested
    % Ensures the response matrix covers the same time span across stimuli,
    % preventing winner's curse inflation in moving window analyses caused by
    % longer stimuli providing more windows to search over.
    % For moving ball (2.3s) vs rectGrid (0.5s), capping at 500ms makes the
    % number of sliding window positions comparable.
    % Warning is issued when capping occurs so the user is aware.
    % -------------------------------------------------------------------------
    if params.CapStimDuration && stimDur > params.MaxStimDuration
        fprintf(['Warning: stimulus duration (%.0f ms) exceeds MaxStimDuration ' ...
                 '(%.0f ms) — capping response window for %s.\n'], ...
                 stimDur, params.MaxStimDuration, obj.stimName);
        effectiveStimDur = params.MaxStimDuration;  % capped duration used for Mr only
    else
        effectiveStimDur = stimDur;  % full duration — no capping needed
    end

    % --- Build spike count matrices ---
    % Mr: response window — capped at effectiveStimDur if CapStimDuration=true
    % Capping takes first MaxStimDuration ms of each trial, starting at stimulus onset
    Mr = BuildBurstMatrix(goodU, ...
        round(p.t), ...
        round(directimesSorted), ...
        round(effectiveStimDur));  % capped or full duration

    % Mb: baseline window — always uses 75% of inter-trial interval
    % Duration is independent of stimulus duration so no capping needed
    Mb = BuildBurstMatrix(goodU, ...
        round(p.t), ...
        round(directimesSorted - 0.75 * obj.VST.interTrialDelay * 1000), ...
        round(0.75 * obj.VST.interTrialDelay * 1000));

    % Always compute full-duration means — needed for empty-trial filtering
    % and for the non-moving-window z-score path regardless of MovingWindow flag
    responses = mean(Mr, 3);   % mean spikes/ms over response window: [nTrials × nNeurons]
    baselines = mean(Mb, 3);   % mean spikes/ms over baseline window: [nTrials × nNeurons]

    % -------------------------------------------------------------------------
    % Compute Diff — method depends on MovingWindow flag
    % MovingWindow=true : per-trial sliding window max, applied identically to
    %                     Mr and Mb so null distribution uses the same operation
    % MovingWindow=false: full duration mean, appropriate for sustained responses
    % -------------------------------------------------------------------------
    if params.MovingWindowPval
        winSize = responseParams.params.durationWindow;  % sliding window size in ms/bins

        % Guard: both response and baseline must be at least as long as winSize
        assert(size(Mr,3) >= winSize, ...
            ['Response window (%d ms) shorter than durationWindow (%d ms). ' ...
             'Reduce durationWindow or increase MaxStimDuration.'], ...
            size(Mr,3), winSize);
        assert(size(Mb,3) >= winSize, ...
            ['Baseline window (%d ms) shorter than durationWindow (%d ms). ' ...
             'Reduce durationWindow or increase interTrialDelay.'], ...
            size(Mb,3), winSize);

        % movmean along dim 3 — vectorised sliding window mean, no loop needed
        % 'Endpoints','discard' removes partial windows at array edges
        mrMov = movmean(Mr, winSize, 3, 'Endpoints', 'discard');  % [nTrials × nNeurons × nStepsR]
        mbMov = movmean(Mb, winSize, 3, 'Endpoints', 'discard');  % [nTrials × nNeurons × nStepsB]

        % Per-trial maximum over all valid window positions: [nTrials × nNeurons]
        responsesMW = max(mrMov, [], 3);
        baselinesMW = max(mbMov, [], 3);

        % Per-trial moving window difference — used for permutation test and z-score
        Diff = responsesMW - baselinesMW;  % [nTrials × nNeurons]

    else
        % Full duration mean — appropriate for sustained/spatially localised responses
        Diff = responses - baselines;  % [nTrials × nNeurons]
    end

    nNeurons = size(goodU, 2);                   % number of good units
    nCats    = round(size(Diff,1) / trialsCat);  % number of stimulus categories

    % Sanity check: total trials must equal nCats × trialsCat
    assert(size(Diff,1) == nCats * trialsCat, ...
        'Trial count (%d) not evenly divisible by trialsCat (%d). Check responseParams.', ...
        size(Diff,1), trialsCat);

    % -------------------------------------------------------------------------
    % Category-level empty-trial filtering
    % Uses full-duration responses (not moving window) for the zero-spike check
    % because MW inflates apparent spike counts for silent trials
    % -------------------------------------------------------------------------
    validCats = true(nCats, nNeurons);  % [nCats × nNeurons]; true = include

    if params.FilterEmptyResponses
        responsesReshaped = reshape(responses, trialsCat, nCats, nNeurons);  % [trialsCat × nCats × nNeurons]
        for c = 1:nCats
            for u = 1:nNeurons
                emptyTrials = responsesReshaped(:, c, u) == 0;  % zero-spike trials in this category
                perc        = sum(emptyTrials) / trialsCat;      % fraction of empty trials
                if perc >= params.EmptyTrialPerc
                    validCats(c, u) = false;  % exclude category c for neuron u
                end
            end
        end
    end

    % Neurons where ALL categories are invalid — statistics undefined
    noValidCat = all(~validCats, 1);  % [1 × nNeurons]

    % -------------------------------------------------------------------------
    % Observed max-statistic
    % Maximum mean Diff across valid categories per neuron [1 × nNeurons]
    % -------------------------------------------------------------------------
    DiffReshaped = reshape(Diff, trialsCat, nCats, nNeurons);         % [trialsCat × nCats × nNeurons]
    catMeans     = reshape(mean(DiffReshaped, 1), nCats, nNeurons);   % [nCats × nNeurons]

    catMeansMasked             = catMeans;
    catMeansMasked(~validCats) = -Inf;            % invalid categories cannot be preferred
    ObsStat                    = max(catMeansMasked, [], 1);  % [1 × nNeurons]

    % -------------------------------------------------------------------------
    % Max-statistic sign-flip permutation test
    %
    % H0: sign of Diff is random per trial (response = baseline distribution).
    % Sign-flipping simulates H0 without parametric assumptions.
    % Taking MAX across categories controls FWER (Nichols & Holmes 2002).
    % Fully vectorised via pagemtimes — no loop over nBoot.
    % -------------------------------------------------------------------------

    % Generate all sign vectors at once: [nTrials × nBoot], values ±1
    signs  = 2 * randi(2, size(Diff,1), params.nBoot) - 3;

    % Reshape signs to match category structure: [trialsCat × nCats × nBoot]
    signsR = reshape(signs, trialsCat, nCats, params.nBoot);

    % Permute for pagemtimes:
    %   DiffRp  : [nNeurons × trialsCat × nCats]
    %   signsRp : [trialsCat × nBoot    × nCats]
    DiffRp  = permute(DiffReshaped, [3 1 2]);
    signsRp = permute(signsR,       [1 3 2]);

    % Batched matrix multiply over category pages: [nNeurons × nBoot × nCats]
    catMeansAll = pagemtimes(DiffRp, signsRp) / trialsCat;

    % Permute to [nCats × nNeurons × nBoot] for masking and max
    catMeansAll = permute(catMeansAll, [3 1 2]);

    % Exclude invalid categories from null distribution
    validCats3D               = repmat(validCats, 1, 1, params.nBoot);  % [nCats × nNeurons × nBoot]
    catMeansAll(~validCats3D) = -Inf;

    % Null distribution: max across categories for each permutation [nBoot × nNeurons]
    % reshape instead of squeeze — safe when nNeurons=1 or nCats=1
    nullMax = reshape(max(catMeansAll, [], 1), params.nBoot, nNeurons);

    % p-value: proportion of null max-statistics >= observed (one-tailed, excitatory)
    pVal             = mean(nullMax >= ObsStat, 1);  % [1 × nNeurons]
    pVal(noValidCat) = NaN;

    % -------------------------------------------------------------------------
    % Permutation z-score
    % Observed stat normalised by the mean and SD of its own null distribution.
    % Answers: "how many SDs above the null is this neuron's observed response?"
    % Saved as a separate field from ZScoreU — the two metrics complement each
    % other and are appropriate for different comparisons.
    % Note: nullSD still partially scales with nCats and stimulus duration,
    % so this metric is not perfectly comparable across stimuli — see methods.
    % -------------------------------------------------------------------------
    nullMean          = mean(nullMax, 1);              % [1 × nNeurons] expected max under H0
    nullSD            = std(nullMax,  1);              % [1 × nNeurons] variability of null max
    zPerm             = (ObsStat - nullMean) ./ nullSD;% [1 × nNeurons] permutation z-score
    zPerm(nullSD==0)  = 0;    % degenerate null — set to 0
    zPerm(noValidCat) = NaN;  % undefined for fully invalid neurons

    % -------------------------------------------------------------------------
    % Data z-score (ZScoreU)
    % Three modes depending on MovingWindow and UseLOO flags:
    %
    % MovingWindow=true:
    %   prefCat = argmax(MW) — MW is [nCats × nNeurons] peak firing rate
    %   per category from sliding window already computed in ResponseWindow.
    %   z_mean = MW at prefCat minus mean baseline (both in spikes/ms).
    %   UseLOO is ignored in this mode.
    %
    % MovingWindow=false, UseLOO=true (recommended):
    %   LOO cross-validated mean Diff at preferred category.
    %   Preferred category identified on n-1 trials per fold.
    %   Prevents winner's curse inflation that scales with nCats.
    %
    % MovingWindow=false, UseLOO=false:
    %   Direct mean Diff at preferred category from all trials.
    %   Faster but inflated when nCats is large — exploration only.
    %
    % All modes normalised by pooled baseline SD across all trials,
    % more stable than per-category SD with few trials per category.
    % -------------------------------------------------------------------------
    sdBase = std(baselines, 0, 1);  % [1 × nNeurons] pooled baseline SD

    % if params.MovingWindowZS
    %     % Preferred category = argmax of MW per neuron
    %     % MW is [nCats × nNeurons] — already the best window mean per category
    %     catMWMasked             = MW;
    %     catMWMasked(~validCats) = -Inf;                     % exclude invalid categories
    %     [~, prefCat]            = max(catMWMasked, [], 1);  % [1 × nNeurons]
    % 
    %     % Extract MW at preferred category using linear indexing
    %     idx_pref  = prefCat + (0:nNeurons-1) * nCats;       % linear index into [nCats × nNeurons]
    %     mwPrefCat = MW(idx_pref);                            % [1 × nNeurons] MW at preferred cat
    %
    %     % Baseline correct: both MW and mean(baselines) are in spikes/ms
    %     z_mean = mwPrefCat - mean(baselines, 1);             % [1 × nNeurons]
    %
    % else
    if nCats == 1
        % Single category — preferred is trivially category 1
        % LOO and direct z-score are identical with one category
        prefCat = ones(1, nNeurons);   % only one category
        z_mean  = mean(Diff, 1);       % mean over all trials [1 × nNeurons]

    elseif params.UseLOO
        % --- LOO cross-validated z-score ---
        % Pre-compute per-category sums for efficient LOO mean: [nCats × nNeurons]
        totalSum = zeros(nCats, nNeurons);
        for c = 1:nCats
            rows          = (c-1)*trialsCat + 1 : c*trialsCat;  % trial rows for category c
            totalSum(c,:) = sum(Diff(rows,:), 1);                % sum over trials
        end

        z_loo_acc    = zeros(1, nNeurons);      % accumulates held-out Diff at preferred cat
        prefCatCount = zeros(nCats, nNeurons);  % tallies preferred category per fold

        for k = 1:trialsCat
            % LOO mean: subtract held-out trial k from total, divide by n-1
            % Indexing (0:nCats-1)*trialsCat+k gives the kth trial of each category
            looMean               = (totalSum - Diff((0:nCats-1)*trialsCat + k, :)) / (trialsCat - 1);
            looMeanMasked         = looMean;
            looMeanMasked(~validCats) = -Inf;               % exclude invalid categories

            [~, prefCatLOO] = max(looMeanMasked, [], 1);   % [1 × nNeurons] preferred cat this fold

            % Linear index into [nCats × nNeurons]
            idx               = prefCatLOO + (0:nNeurons-1) * nCats;
            prefCatCount(idx) = prefCatCount(idx) + 1;     % tally this fold's choice

            % Held-out trial contribution at preferred category
            testVals  = Diff((0:nCats-1)*trialsCat + k, :);  % [nCats × nNeurons]
            z_loo_acc = z_loo_acc + testVals(idx);            % accumulate held-out diff
        end

        z_mean       = z_loo_acc / trialsCat;         % mean held-out diff [1 × nNeurons]
        [~, prefCat] = max(prefCatCount, [], 1);      % consensus preferred category

    else
        % --- Direct z-score (no LOO) ---
        % Subject to winner's curse — use for exploration only
        catMeansDir             = reshape(mean(DiffReshaped, 1), nCats, nNeurons);
        catMeansDir(~validCats) = -Inf;
        [~, prefCat]            = max(catMeansDir, [], 1);         % [1 × nNeurons]
        idx                     = prefCat + (0:nNeurons-1) * nCats;
        z_mean                  = catMeansDir(idx);                % [1 × nNeurons]
    end
    % end

    % Normalise by pooled baseline SD
    z              = z_mean ./ sdBase;  % [1 × nNeurons]
    z(sdBase == 0) = 0;    % silent baseline — set to 0
    z(noValidCat)  = NaN;  % no valid categories — undefined

    % -------------------------------------------------------------------------
    % One-sample t-test pooled across all valid categories
    % H0: mean(Diff) = 0 across all valid trials.
    % Pooling maximises df and avoids cherry-picking the preferred category.
    % Permutation test is the primary criterion; t-test is a secondary check.
    % -------------------------------------------------------------------------
    pValTTest = zeros(1, nNeurons);
    tStat     = zeros(1, nNeurons);

    for u = 1:nNeurons
        if noValidCat(u)
            pValTTest(u) = NaN;
            tStat(u)     = NaN;
            continue
        end

        % Logical row mask: all trials belonging to valid categories for neuron u
        validRows = false(size(Diff, 1), 1);
        for c = 1:nCats
            if validCats(c, u)
                rows            = (c-1)*trialsCat + 1 : c*trialsCat;
                validRows(rows) = true;
            end
        end

        DiffValid            = Diff(validRows, u);              % valid trials for neuron u
        [~, pValTTest(u), ~, stats] = ttest(DiffValid);        % one-sample t-test vs zero
        tStat(u)             = stats.tstat;
    end

    pValTTest(noValidCat) = NaN;
    tStat(noValidCat)     = NaN;

    % -------------------------------------------------------------------------
    % Store results for this condition
    % -------------------------------------------------------------------------
    if isfield(responseParams, "Speed1") || isequal(obj.stimName, 'StaticDriftingGrating')
        S.(fieldName).pvalsResponse      = pVal;        % [1 × nNeurons] permutation p-values
        S.(fieldName).ZScoreU            = z;            % [1 × nNeurons] data z-score (LOO/direct/MW)
        S.(fieldName).ZScorePermutation  = zPerm;        % [1 × nNeurons] permutation z-score
        S.(fieldName).ObsDiff            = Diff;         % [nTrials × nNeurons] response minus baseline
        S.(fieldName).ObsResponse        = responses;    % [nTrials × nNeurons] full-duration response
        S.(fieldName).ObsBaseline        = baselines;    % [nTrials × nNeurons] baseline spike counts
        S.(fieldName).prefCat            = prefCat;      % [1 × nNeurons] preferred category index
        S.(fieldName).validCats          = validCats;    % [nCats × nNeurons] category validity mask
        S.(fieldName).MaxMovWinResponse  = max(MW,[],1); % [1 × nNeurons] peak MW response across cats
        S.(fieldName).pValTTest          = pValTTest;    % [1 × nNeurons] t-test p-values
        S.(fieldName).tStat              = tStat;        % [1 × nNeurons] t-statistics
    else
        S.pvalsResponse     = pVal;
        S.ZScoreU           = z;
        S.ZScorePermutation = zPerm;
        S.ObsDiff           = Diff;
        S.ObsResponse       = responses;
        S.ObsBaseline       = baselines;
        S.prefCat           = prefCat;
        S.validCats         = validCats;
        S.MaxMovWinResponse = max(MW,[],1);
        S.pValTTest         = pValTTest;
        S.tStat             = tStat;
    end

    S.params = params;  % store parameters alongside results for reproducibility

end % end condition loop

% --- Save and return ---
fprintf('Saving results to file.\n');
save(obj.getAnalysisFileName, '-struct', 'S');
results = S;

end % end main function


% =========================================================================
%% Local helper: build empty output struct when no neurons are found
% =========================================================================
function S = buildEmptyStruct(obj, responseParams)
% buildEmptyStruct - Returns empty results struct with correct field names.
% Ensures downstream code receives a consistent struct regardless of neuron count.

    emptyFields = {'pvalsResponse','ZScoreU','ZScorePermutation','ObsDiff', ...
                   'ObsResponse','ObsBaseline','prefCat','validCats', ...
                   'MaxMovWinResponse','pValTTest','tStat'};

    if isequal(obj.stimName, 'linearlyMovingBall') || isequal(obj.stimName, 'linearlyMovingBar')
        for f = emptyFields
            S.Speed1.(f{1}) = [];
        end
        if isfield(responseParams, "Speed2")
            for f = emptyFields
                S.Speed2.(f{1}) = [];
            end
        end

    elseif isequal(obj.stimName, 'StaticDriftingGrating')
        for cond = {'Static', 'Moving'}
            for f = emptyFields
                S.(cond{1}).(f{1}) = [];
            end
        end

    else
        for f = emptyFields
            S.(f{1}) = [];
        end
    end
end