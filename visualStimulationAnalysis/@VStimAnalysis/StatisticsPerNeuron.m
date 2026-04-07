function results = StatisticsPerNeuron(obj, params)
% StatisticsPerNeuron - Computes per-neuron response statistics vs baseline.
%
% For each neuron this function outputs:
%   pVal    : p-value from a max-statistic sign-flip permutation test.
%             Tests H0: no stimulus category drives a response above baseline.
%             The max-statistic controls family-wise error rate across categories
%             without requiring Bonferroni correction.
%
%   ZScoreU : Leave-one-out (LOO) cross-validated z-score at the preferred
%             stimulus category. On each LOO fold, the preferred category is
%             identified on all-but-one trials, and the held-out trial contributes
%             to the z-score estimate. This prevents winner's curse inflation that
%             would otherwise scale with nCats, making z-scores non-comparable
%             across stimuli with different category counts (e.g. rectGrid with
%             81 positions vs moving ball with 4 directions).
%
%   prefCat : Consensus preferred category — the category most frequently
%             selected as preferred across all LOO folds.
%
%   validCats: [nCats × nNeurons] logical mask. False where a category has
%             >= EmptyTrialPerc fraction of zero-spike trials for a given neuron.
%
% Usage:
%   results = obj.StatisticsPerNeuron()
%   results = obj.StatisticsPerNeuron(nBoot=5000, overwrite=true)
%
% Reference for sign-flip permutation test:
%   Nichols & Holmes (2002) Human Brain Mapping 15:1-25

arguments (Input)
    obj
    params.nBoot               = 10000  % number of permutation iterations for null distribution
    params.EmptyTrialPerc      = 0.7    % exclude category if fraction of zero-spike trials >= this threshold
    params.FilterEmptyResponses = true  % whether to apply empty-trial category filtering
    params.overwrite           = false  % if true, recompute even if a saved file already exists
    params.randomSeed          = 42     % fixed seed for reproducible permutation results (required for publication)
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
%   linearlyMovingBall   → one or two speed conditions (Speed1, Speed2)
%   StaticDriftingGrating → Static and Moving phases
%   all others (rectGrid) → single condition
% -------------------------------------------------------------------------
if isfield(responseParams, "Speed1")
    % BUG FIX: original code used length(obj.VST.speed) which returns total
    % number of trials, not unique speeds — caused loop to run hundreds of times
    nSpeeds = numel(unique(obj.VST.speed));  % number of distinct speed values

    Times.Speed1      = responseParams.Speed1.C(:,1)';                                    % trial onset times [1 × nTrials]
    Durations.Speed1  = responseParams.Speed1.stimDur;                                    % stimulus duration in ms
    trialsCats.Speed1 = round(numel(Times.Speed1) / size(responseParams.Speed1.NeuronVals, 2));  % trials per category

    if nSpeeds > 1
        Times.Speed2      = responseParams.Speed2.C(:,1)';
        Durations.Speed2  = responseParams.Speed2.stimDur;
        trialsCats.Speed2 = round(numel(Times.Speed2) / size(responseParams.Speed2.NeuronVals, 2));
    end

    x = nSpeeds;  % number of loop iterations

elseif isequal(obj.stimName, 'StaticDriftingGrating')
    % Moving phase onset is shifted by static_time relative to trial onset
    Times.Moving      = responseParams.C(:,1)' + obj.VST.static_time * 1000;
    Durations.Moving  = responseParams.Moving.stimDur;
    trialsCats.Moving = round(numel(Times.Moving) / size(responseParams.Moving.NeuronVals, 2));

    Times.Static      = responseParams.C(:,1)';
    Durations.Static  = responseParams.Static.stimDur;
    trialsCats.Static = round(numel(Times.Static) / size(responseParams.Static.NeuronVals, 2));

    FieldNames = {'Static', 'Moving'};  % loop will index these in order
    x = 2;

else
    % Single-condition stimuli (rectGrid, etc.)
    directimesSorted = responseParams.C(:,1)';
    stimDur          = responseParams.stimDur;
    trialsCat        = round(numel(directimesSorted) / size(responseParams.NeuronVals, 2));
    x = 1;
end

% =========================================================================
% Main loop over stimulus conditions
% =========================================================================
for s = 1:x

    % --- Assign condition-specific timing variables ---
    if isfield(responseParams, "Speed1")
        fieldName        = sprintf('Speed%d', s);   % 'Speed1' or 'Speed2'
        directimesSorted = Times.(fieldName);
        stimDur          = Durations.(fieldName);
        trialsCat        = trialsCats.(fieldName);
    end

    if isequal(obj.stimName, 'StaticDriftingGrating')
        fieldName        = FieldNames{s};           % 'Static' or 'Moving'
        directimesSorted = Times.(fieldName);
        stimDur          = Durations.(fieldName);
        trialsCat        = trialsCats.(fieldName);
    end

    % --- Build spike count matrices ---
    % Mr: spike counts in the response window (stimulus duration)
    Mr = BuildBurstMatrix(goodU, ...
        round(p.t), ...
        round(directimesSorted), ...
        round(stimDur));

    % Mb: spike counts in the baseline window (75% of inter-trial interval
    % before each trial onset — conservative buffer to avoid overlap with
    % the preceding stimulus)
    Mb = BuildBurstMatrix(goodU, ...
        round(p.t), ...
        round(directimesSorted - 0.75 * obj.VST.interTrialDelay * 1000), ...
        round(0.75 * obj.VST.interTrialDelay * 1000));

    responses = mean(Mr, 3);        % mean spike count over time bins: [nTrials × nNeurons]
    baselines = mean(Mb, 3);        % mean spike count over baseline bins: [nTrials × nNeurons]
    Diff      = responses - baselines;  % per-trial response minus baseline: [nTrials × nNeurons]

    nNeurons = size(goodU, 2);              % number of good units
    nCats    = round(size(Diff,1) / trialsCat);  % number of stimulus categories

    % Sanity check: total trials must equal nCats × trialsCat
    assert(size(Diff,1) == nCats * trialsCat, ...
        'Trial count (%d) is not evenly divisible by trialsCat (%d). Check responseParams.', ...
        size(Diff,1), trialsCat);

    % -------------------------------------------------------------------------
    % Category-level empty-trial filtering
    % Mark a category as invalid for a given neuron if the fraction of trials
    % with zero spikes meets or exceeds EmptyTrialPerc.
    % Invalid categories are excluded (not zeroed) to avoid biasing Diff toward 0.
    % -------------------------------------------------------------------------
    validCats = true(nCats, nNeurons);  % [nCats × nNeurons]; true = include in analysis

    if params.FilterEmptyResponses
        % Reshape responses to [trialsCat × nCats × nNeurons] for category indexing
        responsesReshaped = reshape(responses, trialsCat, nCats, nNeurons);

        for c = 1:nCats
            for u = 1:nNeurons
                emptyTrials = responsesReshaped(:, c, u) == 0;   % logical: zero-spike trials
                perc        = sum(emptyTrials) / trialsCat;       % fraction of empty trials
                if perc >= params.EmptyTrialPerc
                    validCats(c, u) = false;  % exclude category c for neuron u
                end
            end
        end
    end

    % Neurons where ALL categories are invalid — statistics are undefined
    noValidCat = all(~validCats, 1);  % [1 × nNeurons]

    % -------------------------------------------------------------------------
    % Reshape Diff into category structure
    % DiffReshaped: [trialsCat × nCats × nNeurons]
    % -------------------------------------------------------------------------
    DiffReshaped = reshape(Diff, trialsCat, nCats, nNeurons);

    % -------------------------------------------------------------------------
    % Observed max-statistic
    % For each neuron: maximum mean response-minus-baseline across valid categories.
    % Invalid categories are set to -Inf so they cannot contribute to the max.
    % -------------------------------------------------------------------------
    catMeans            = squeeze(mean(DiffReshaped, 1));   % [nCats × nNeurons]
    catMeansMasked      = catMeans;
    catMeansMasked(~validCats) = -Inf;                      % exclude invalid categories
    ObsStat             = max(catMeansMasked, [], 1);        % [1 × nNeurons]

    % -------------------------------------------------------------------------
    % Max-statistic sign-flip permutation test
    %
    % H0: for each trial the sign of (response - baseline) is random,
    %     meaning response and baseline are drawn from the same distribution.
    %
    % Under H0, randomly flipping the sign of each trial's difference is
    % equivalent to randomly reassigning which window is "response" and which
    % is "baseline". Repeating this nBoot times builds a null distribution of
    % the max category mean under H0.
    %
    % Taking the MAX across categories at each permutation automatically
    % controls family-wise error rate (FWER) across categories without
    % requiring Bonferroni correction (Nichols & Holmes 2002).
    %
    % Fully vectorised using pagemtimes for efficiency — no loop over nBoot.
    % -------------------------------------------------------------------------

    % Generate all sign vectors at once: [nTrials × nBoot], values ±1
    signs = 2 * randi(2, size(Diff,1), params.nBoot) - 3;

    % Reshape signs to match category structure: [trialsCat × nCats × nBoot]
    % This maps each trial's sign flip to its correct category slot
    signsR = reshape(signs, trialsCat, nCats, params.nBoot);

    % Permute for pagemtimes — pages correspond to categories:
    %   DiffRp  : [nNeurons × trialsCat × nCats]
    %   signsRp : [trialsCat × nBoot    × nCats]
    DiffRp  = permute(DiffReshaped, [3 1 2]);   % [nNeurons × trialsCat × nCats]
    signsRp = permute(signsR,       [1 3 2]);   % [trialsCat × nBoot    × nCats]

    % Batched matrix multiply over category pages:
    %   for each category: [nNeurons × trialsCat] × [trialsCat × nBoot] = [nNeurons × nBoot]
    %   result: [nNeurons × nBoot × nCats]
    catMeansAll = pagemtimes(DiffRp, signsRp) / trialsCat;

    % Permute to [nCats × nNeurons × nBoot] for masking and max operation
    catMeansAll = permute(catMeansAll, [3 1 2]);

    % Exclude invalid categories from null distribution (same mask as observed stat)
    validCats3D = repmat(validCats, 1, 1, params.nBoot);  % [nCats × nNeurons × nBoot]
    catMeansAll(~validCats3D) = -Inf;

    % Max across categories for each permutation and neuron: [nBoot × nNeurons]
    nullMax = squeeze(max(catMeansAll, [], 1))';   % squeeze: [nNeurons × nBoot], then transpose

    % p-value: proportion of null max-statistics >= observed max-statistic
    % One-tailed test for excitatory responses
    pVal             = mean(nullMax >= ObsStat, 1);  % [1 × nNeurons]
    pVal(noValidCat) = NaN;                          % undefined for neurons with no valid categories

    % -------------------------------------------------------------------------
    % Leave-one-out (LOO) cross-validated z-score
    %
    % With only ~10 trials per category, split-half would leave only 5 trials
    % for each half — too few for stable estimates. LOO maximises data usage:
    %   - Selection set: all (trialsCat - 1) trials → identify preferred category
    %   - Test set:      the single held-out trial  → contributes to z-score
    %
    % The preferred category is selected independently of the test trial on
    % every fold, preventing winner's curse inflation that would otherwise
    % differ across stimuli with different numbers of categories.
    %
    % Implementation is vectorised over categories and neurons.
    % The only loop runs trialsCat times (e.g., 10) — negligible cost.
    % -------------------------------------------------------------------------

    % Pre-compute sum across trials once: [1 × nCats × nNeurons]
    % Subtracting one trial from the total is cheaper than recomputing the mean
    totalSum = sum(DiffReshaped, 1);

    % Pre-allocate accumulators
    z_loo_acc    = zeros(1, nNeurons);      % accumulates held-out diff at preferred category
    prefCatCount = zeros(nCats, nNeurons);  % tallies how often each category is preferred per fold

    for k = 1:trialsCat
        % Category mean on all trials except trial k: [nCats × nNeurons]
        % Subtracting trial k from pre-computed total avoids recomputing the full mean
        looMean = squeeze((totalSum - DiffReshaped(k,:,:)) / (trialsCat - 1));

        % Exclude invalid categories from selection
        looMeanMasked             = looMean;
        looMeanMasked(~validCats) = -Inf;

        % Preferred category index on the selection data: [1 × nNeurons]
        [~, prefCatLOO] = max(looMeanMasked, [], 1);

        % Accumulate preferred category tally (used later for consensus prefCat)
        idx            = sub2ind([nCats, nNeurons], prefCatLOO, 1:nNeurons);
        prefCatCount(idx) = prefCatCount(idx) + 1;  % increment tally for this fold's choice

        % Held-out trial k: diff values at all categories: [nCats × nNeurons]
        testVals = squeeze(DiffReshaped(k,:,:));

        % Accumulate the held-out diff at the fold's preferred category
        z_loo_acc = z_loo_acc + testVals(idx);  % [1 × nNeurons]
    end

    % Average LOO diff across all held-out trials: [1 × nNeurons]
    z_loo_mean = z_loo_acc / trialsCat;

    % Consensus preferred category: most frequently selected across LOO folds
    % Used to extract baseline SD at a stable preferred location
    [~, prefCat] = max(prefCatCount, [], 1);  % [1 × nNeurons]

    % Baseline SD pooled across all trials and categories per neuron
    % More stable than per-category SD, especially with few trials per category.
    % Justified because baseline periods are pre-stimulus and should not vary
    % systematically across stimulus categories.
    sdBasePref = std(baselines, 0, 1);  % [1 × nNeurons] — std over all nTrials rows

    % Final z-score: LOO mean diff normalised by baseline SD
    z                  = z_loo_mean ./ sdBasePref;  % [1 × nNeurons]
    z(sdBasePref == 0) = 0;    % silent baseline (no variability): set to 0
    z(noValidCat)      = NaN;  % no valid categories: undefined

    % -------------------------------------------------------------------------
    % Store results for this condition
    % -------------------------------------------------------------------------
    if isfield(responseParams, "Speed1") || isequal(obj.stimName, 'StaticDriftingGrating')
        % Named sub-struct for multi-condition stimuli
        S.(fieldName).pvalsResponse = pVal;       % [1 × nNeurons] permutation p-values
        S.(fieldName).ZScoreU       = z;           % [1 × nNeurons] LOO cross-validated z-scores
        S.(fieldName).ObsDiff       = Diff;        % [nTrials × nNeurons] raw trial differences
        S.(fieldName).ObsResponse   = responses;   % [nTrials × nNeurons] response spike counts
        S.(fieldName).ObsBaseline   = baselines;   % [nTrials × nNeurons] baseline spike counts
        S.(fieldName).prefCat       = prefCat;     % [1 × nNeurons] consensus preferred category index
        S.(fieldName).validCats     = validCats;   % [nCats × nNeurons] category validity mask
    else
        % Flat struct for single-condition stimuli
        S.pvalsResponse = pVal;
        S.ZScoreU       = z;
        S.ObsDiff       = Diff;
        S.ObsResponse   = responses;
        S.ObsBaseline   = baselines;
        S.prefCat       = prefCat;
        S.validCats     = validCats;
    end

    S.params = params;  % store analysis parameters for reproducibility

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
% buildEmptyStruct - Returns an empty results struct with correct field names.
% Ensures downstream code receives a consistent struct regardless of neuron count.

    emptyFields = {'pvalsResponse','ZScoreU','ObsDiff','ObsResponse', ...
                   'ObsBaseline','prefCat','validCats'};

    if isequal(obj.stimName, 'linearlyMovingBall')
        for f = emptyFields
            S.Speed1.(f{1}) = [];  % empty Speed1 fields
        end
        if isfield(responseParams, "Speed2")
            for f = emptyFields
                S.Speed2.(f{1}) = [];  % empty Speed2 fields if second speed exists
            end
        end

    elseif isequal(obj.stimName, 'StaticDriftingGrating')
        for cond = {'Static', 'Moving'}
            for f = emptyFields
                S.(cond{1}).(f{1}) = [];  % empty fields for each grating condition
            end
        end

    else
        for f = emptyFields
            S.(f{1}) = [];  % flat empty struct for single-condition stimuli
        end
    end
end