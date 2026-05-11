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
% Output is saved to a separate file: analysisFileName_categoryname.mat
%
% Category names are matched case-insensitively to responseParams.colNames.
% For linearlyMovingBall, comparing 'speed' receives special handling since
% Speed1 and Speed2 are stored in separate struct fields.
%
% Usage:
%   results = obj.StatisticsPerNeuronPerCategory('compareCategory', 'size')
%   results = obj.StatisticsPerNeuronPerCategory('compareCategory', 'direction', ...
%               'nBoot', 5000, 'overwrite', true)
%
% Reference for permutation tests:
%   Nichols & Holmes (2002) Human Brain Mapping 15:1-25

arguments (Input)
    obj
    params.compareCategory  = ''      % category name to compare (case-insensitive)
    params.nBoot            = 10000   % permutation iterations
    params.BaseRespWindow   = 1000     % ms response window from stimulus onset
    params.BaselineBuffer   = 200     % ms buffer before stimulus onset for baseline
                                      % avoids contamination from off-responses of
                                      % preceding stimulus or anticipatory activity
    params.overwrite        = false   % recompute even if saved file exists
    params.randomSeed       = 42      % fixed seed for reproducibility
    params.ApplyFDR         = false    % Benjamini-Hochberg FDR correction for pairwise
    params.MovingWindowDuration = 200   % ms sliding window for moving ball per-trial peak response
                                     % Applied to full stimulus duration, response only (not baseline)
                                     % Only used when stimulus is linearlyMovingBall
    params.GratingType       = "moving"  %If the stimulus is grating, select it's mode.                                  
end

% -------------------------------------------------------------------------
% Validate input
% -------------------------------------------------------------------------
if isempty(strtrim(params.compareCategory))
    error('params.compareCategory must be specified (e.g. ''size'', ''direction'').');
end

% -------------------------------------------------------------------------
% Output file: append category name to base analysis filename
% -------------------------------------------------------------------------
outputFile = strrep(obj.getAnalysisFileName, '.mat', ...
    ['_' lower(strtrim(params.compareCategory)) '.mat']);

if isfile(outputFile) && ~params.overwrite
    if nargout == 1
        fprintf('Loading saved results from file.\n');
        results = load(outputFile);
    else
        fprintf('Analysis already exists (use overwrite option to recalculate).\n');
    end
    return
end

% -------------------------------------------------------------------------
% Fix random seed for reproducibility
% -------------------------------------------------------------------------
rng(params.randomSeed);

% -------------------------------------------------------------------------
% Load spike-sorted somatic units
% -------------------------------------------------------------------------
p        = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
label    = string(p.label');
goodU    = p.ic(:, label == 'good');
nNeurons = size(goodU, 2);

if isempty(goodU)
    warning('%s has no somatic neurons.', obj.dataObj.recordingName);
    results = [];
    return
end

responseParams = obj.ResponseWindow;

% -------------------------------------------------------------------------
% Sync diode triggers for stimulus alignment
% -------------------------------------------------------------------------
try
    obj.getSyncedDiodeTriggers;
catch
    obj.getSessionTime("overwrite", true);
    obj.getDiodeTriggers("extractionMethod", 'digitalTriggerDiode', 'overwrite', true);
    obj.getSyncedDiodeTriggers;
end

% -------------------------------------------------------------------------
% Identify stimulus type and set flags
% -------------------------------------------------------------------------
isMovingBall = isequal(obj.stimName, 'linearlyMovingBall') || ...
               isequal(obj.stimName, 'linearlyMovingBar');
isGratingMov    = isequal(obj.stimName, 'StaticDriftingGrating') && params.GratingType == "moving";
isGratingStat    = isequal(obj.stimName, 'StaticDriftingGrating') && params.GratingType == "static";
isSpeedComp  = isMovingBall && strcmpi(strtrim(params.compareCategory), 'speed');

% -------------------------------------------------------------------------
% Get C matrix, trial times, stimulus duration and category column names
% -------------------------------------------------------------------------
if isMovingBall
    nSpeeds = numel(unique(obj.VST.speed));

    if isSpeedComp
        % Speed comparison: each speed is a level, handled separately below
 

        if nSpeeds < 2
            fprintf(['Only one speed condition found in %s. ' ...
                'Cannot compare speeds.\n'], obj.stimName);
            results = [];
            return
        end

        nLevels = nSpeeds;
        levels  = (1:nSpeeds)';
        fprintf('Comparing %d speed conditions for %s.\n', nSpeeds, obj.stimName);
    else
        % colNames are the same regardless of speed — just need them for category matching
        colNames = responseParams.colNames{1}(5:end);
        % C will be overwritten with pooled version inside the response matrix block
        C = responseParams.Speed1.C;  % temporary — used only for catIdx/cCol detection
    end

elseif isGratingMov
    % Use Moving phase for grating
    C        = responseParams.C;
    C(:,1) = C(:,1) +obj.VST.static_time*1000;
    colNames = responseParams.colNames{1}(5:end);

else
    % All other stimuli (rectGrid, etc.)
    C        = responseParams.C;
    colNames = responseParams.colNames{1}(5:end);
end

% -------------------------------------------------------------------------
% Find category column in C and get unique levels
% colNames{k} corresponds to C(:, k+1) since C(:,1) is stimulus onset time
% -------------------------------------------------------------------------
if ~isSpeedComp
    catIdx = find(strcmpi(colNames, strtrim(params.compareCategory)));

    if isempty(catIdx)
        fprintf(['Category "%s" not found in this stimulus.\n' ...
                 'Available categories: %s\n'], ...
                 params.compareCategory, strjoin(colNames, ', '));
        results = [];
        return
    end

    cCol    = catIdx + 1;                 % column index in C (C(:,2) = first category)
    levels  = unique(C(:, cCol));         % unique level values [nLevels × 1]
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

if isSpeedComp
    allDiff      = cell(nSpeeds, 1);
    allBaselines = cell(nSpeeds, 1);

    for s = 1:nSpeeds
        fName_s      = sprintf('Speed%d', s);
        trialTimes_s = responseParams.(fName_s).C(:,1)';
        stimDur_s    = responseParams.(fName_s).stimDur;

        % Response: sliding window across full stimulus duration (moving ball always)
        Mr_s = BuildBurstMatrix(goodU, round(p.t), ...
            round(trialTimes_s), round(stimDur_s));

        assert(size(Mr_s,3) >= params.MovingWindowDuration, ...
            'Speed%d stimulus duration (%d ms) shorter than MovingWindowDuration (%d ms).', ...
            s, size(Mr_s,3), params.MovingWindowDuration);

        mrMov_s     = movmean(Mr_s, params.MovingWindowDuration, 3, ...
            'Endpoints', 'discard');
        responses_s = max(mrMov_s, [], 3);              % [nTrials_s × nNeurons] peak window

        % Baseline: fixed window — no moving window on baseline
        Mb_s = BuildBurstMatrix(goodU, round(p.t), ...
            round(trialTimes_s - min([obj.VST.interTrialDelay*1000 - params.BaselineBuffer, ...
            params.BaseRespWindow])), ...
            params.BaseRespWindow);

        baselines_s     = mean(Mb_s, 3);
        allDiff{s}      = responses_s - baselines_s;
        allBaselines{s} = baselines_s;

        fprintf('Speed%d: using %d ms sliding window over %d ms stimulus.\n', ...
            s, params.MovingWindowDuration, round(stimDur_s));
    end

    % Pooled baseline SD across all trials and speeds
    sdBase = std(vertcat(allBaselines{:}), 0, 1);   % [1 × nNeurons]

else
    % -------------------------------------------------------------------------
    % Standard comparison: build full matrices once, split by category level
    % For moving ball: sliding window across full stimulus duration to capture
    %   peak per-trial response regardless of when ball crosses the RF.
    %   Baseline remains fixed — no moving window on baseline to avoid false
    %   negative bias (supervisor recommendation).
    % For all other stimuli: fixed window from stimulus onset.
    % -------------------------------------------------------------------------

    if isMovingBall
        % Pool trials across ALL speeds instead of using max speed only
        % Each speed has different stimDur so build matrices per speed,
        % apply moving window to each, then concatenate
        nSpeeds = numel(unique(obj.VST.speed));

        % Concatenate C matrices from all speeds to get pooled category info
        C_all = [];
        for s = 1:nSpeeds
            fName_s = sprintf('Speed%d', s);
            C_all   = [C_all; responseParams.(fName_s).C];
        end
        C    = C_all;            % overwrite C with pooled version
        cCol = catIdx + 1;       % recalculate in case C structure changed

        % Rebuild levels from pooled C (should be same but ensures consistency)
        levels  = unique(C(:, cCol));
        nLevels = numel(levels);

        % Build response and baseline per speed, apply moving window, concatenate
        responsesList = cell(nSpeeds, 1);
        baselinesList = cell(nSpeeds, 1);

        for s = 1:nSpeeds
            fName_s      = sprintf('Speed%d', s);
            trialTimes_s = responseParams.(fName_s).C(:,1)';
            stimDur_s    = responseParams.(fName_s).stimDur;

            % Response: full stimulus duration with sliding window
            Mr_s = BuildBurstMatrix(goodU, round(p.t), ...
                round(trialTimes_s), round(stimDur_s));

            assert(size(Mr_s,3) >= params.MovingWindowDuration, ...
                'Speed%d: stimulus (%d ms) shorter than MovingWindowDuration (%d ms).', ...
                s, size(Mr_s,3), params.MovingWindowDuration);

            mrMov_s          = movmean(Mr_s, params.MovingWindowDuration, 3, ...
                'Endpoints', 'discard');
            responsesList{s} = max(mrMov_s, [], 3);   % [nTrials_s × nNeurons]

            % Baseline: fixed window — same approach for all speeds
            Mb_s = BuildBurstMatrix(goodU, round(p.t), ...
                round(trialTimes_s - min([obj.VST.interTrialDelay*1000 - params.BaselineBuffer, ...
                params.BaseRespWindow])), ...
                params.BaseRespWindow);
            baselinesList{s} = mean(Mb_s, 3);          % [nTrials_s × nNeurons]

            fprintf('Speed%d: %d ms sliding window over %d ms, %d trials pooled.\n', ...
                s, params.MovingWindowDuration, round(stimDur_s), size(Mr_s,1));
        end

        % Concatenate across speeds: [nTotalTrials × nNeurons]
        responsesFull = vertcat(responsesList{:});
        baselinesFull = vertcat(baselinesList{:});
        DiffFull      = responsesFull - baselinesFull;

        % Pooled baseline SD across all trials and speeds
        sdBase = std(baselinesFull, 0, 1);   % [1 × nNeurons]

        % Split Diff by category level using pooled C
        allDiff = cell(nLevels, 1);
        for k = 1:nLevels
            mask       = C(:, cCol) == levels(k);
            allDiff{k} = DiffFull(mask, :);
            fprintf('Level %g: %d trials (pooled across speeds)\n', levels(k), sum(mask));
        end

    else

        trialTimes = C(:,1)';

        % Fixed window for all other stimuli
        MrFull        = BuildBurstMatrix(goodU, round(p.t), ...
            round(trialTimes), params.BaseRespWindow);
        responsesFull = mean(MrFull, 3);                            % [nTrials × nNeurons]


        % Baseline: fixed window ending BaselineBuffer ms before onset
        % Same for all stimuli — no moving window on baseline
        MbFull        = BuildBurstMatrix(goodU, round(p.t), ...
            round(trialTimes - min([obj.VST.interTrialDelay*1000 - params.BaselineBuffer, ...
            params.BaseRespWindow])), ...
            params.BaseRespWindow);
        baselinesFull = mean(MbFull, 3);                                % [nTrials × nNeurons]

        DiffFull      = responsesFull - baselinesFull;                  % [nTrials × nNeurons]

        % Pooled baseline SD across all trials
        sdBase = std(baselinesFull, 0, 1);   % [1 × nNeurons]

        % Split Diff by category level — pool all other non-specified categories
        allDiff = cell(nLevels, 1);
        for k = 1:nLevels
            mask       = C(:, cCol) == levels(k);
            allDiff{k} = DiffFull(mask, :);
            fprintf('Level %g: %d trials\n', levels(k), sum(mask));
        end

    end
end

% =========================================================================
% Per-level responsiveness test
% Sign-flip permutation test: H0: mean(Diff) = 0 for this level
% Trials are pooled across all non-specified categories within each level
% One-tailed (excitatory): p = proportion of null >= observed mean
% Z-score: bias-corrected by subtracting null mean, normalised by sdBase
% =========================================================================
pValPerLevel  = nan(nLevels, nNeurons);    % [nLevels × nNeurons] p-values
zPerLevel     = nan(nLevels, nNeurons);    % [nLevels × nNeurons] z-scores
obsStatLevels = nan(nLevels, nNeurons);    % [nLevels × nNeurons] observed mean Diff

for k = 1:nLevels
    Diff_k    = allDiff{k};                 % [nTrials_k × nNeurons]
    nTrials_k = size(Diff_k, 1);

    ObsStat_k          = mean(Diff_k, 1);   % [1 × nNeurons]
    obsStatLevels(k,:) = ObsStat_k;

    % Vectorised sign-flip null distribution: [nBoot × nNeurons]
    signs_k    = 2 * randi(2, nTrials_k, params.nBoot) - 3;  % [nTrials_k × nBoot]
    nullDist_k = (signs_k' * Diff_k) / nTrials_k;             % [nBoot × nNeurons]

    % One-tailed p-value: proportion of null >= observed
    pValPerLevel(k,:) = mean(nullDist_k >= ObsStat_k, 1);

    % Bias-corrected z-score: (observed - null mean) / pooled baseline SD
    nullMean_k           = mean(nullDist_k, 1);                % [1 × nNeurons]
    z_k                  = (ObsStat_k - nullMean_k) ./ sdBase; % [1 × nNeurons]
    z_k(sdBase == 0)     = 0;
    zPerLevel(k,:)       = z_k;
end

% =========================================================================
% Omnibus test: permutation one-way ANOVA F-test
% H0: mean response is equal across all levels
% Pool all trials, permute level labels nBoot times
% More powerful than pairwise-only approach since it uses all data
% =========================================================================
DiffPooled      = vertcat(allDiff{:});   % [nTotalTrials × nNeurons]
nTotalTrials    = size(DiffPooled, 1);

% Level label vector: k repeated nTrials_k times per level
levelLabelsVec = cell2mat(arrayfun(@(k) ...
    k * ones(size(allDiff{k},1), 1), ...
    (1:nLevels)', 'UniformOutput', false));   % [nTotalTrials × 1]

% Observed F-statistic
F_obs = computeFstat(DiffPooled, levelLabelsVec, levels);   % [1 × nNeurons]

% Null distribution: permute level labels
nullF = zeros(params.nBoot, nNeurons);
for b = 1:params.nBoot
    permLabels = levelLabelsVec(randperm(nTotalTrials));
    nullF(b,:) = computeFstat(DiffPooled, permLabels, levels);
end

pValOmnibus = mean(nullF >= F_obs, 1);   % [1 × nNeurons]

% =========================================================================
% Pairwise comparisons: two-sample permutation test for each level pair
% Observed: mean(Diff_i) - mean(Diff_j)
% Null: randomly reassign trials between groups, recompute mean difference
% Two-tailed: |observed| >= |null|
% FDR correction across all pairs × neurons
% =========================================================================
pairIdx         = nchoosek(1:nLevels, 2);   % [nPairs × 2]
nPairs          = size(pairIdx, 1);

pValPairwise    = nan(nPairs, nNeurons);    % raw p-values
obsStatPairwise = nan(nPairs, nNeurons);    % observed mean differences

for pr = 1:nPairs
    i      = pairIdx(pr,1);
    j      = pairIdx(pr,2);
    Diff_i = allDiff{i};                    % [nTrials_i × nNeurons]
    Diff_j = allDiff{j};                    % [nTrials_j × nNeurons]
    nI     = size(Diff_i, 1);
    nJ     = size(Diff_j, 1);

    ObsDiff_ij            = mean(Diff_i,1) - mean(Diff_j,1);   % [1 × nNeurons]
    obsStatPairwise(pr,:) = ObsDiff_ij;

    % Pool trials and permute group assignment
    DiffPair   = [Diff_i; Diff_j];           % [nI+nJ × nNeurons]
    nTotal_pr  = nI + nJ;

    % Vectorised: weight matrix W [nTotal_pr × nBoot]
    % For each boot: first nI rows get weight +1/nI, rest get -1/nJ
    nullPair = zeros(params.nBoot, nNeurons);
    for b = 1:params.nBoot
        perm          = randperm(nTotal_pr);
        nullPair(b,:) = mean(DiffPair(perm(1:nI),:),    1) - ...
                        mean(DiffPair(perm(nI+1:end),:), 1);
    end

    % Two-tailed p-value
    pValPairwise(pr,:) = mean(abs(nullPair) >= abs(ObsDiff_ij), 1);
end

% FDR correction across all pairs × neurons simultaneously
if params.ApplyFDR && nPairs > 1
    pFlat     = pValPairwise(:);               % [nPairs*nNeurons × 1]
    validMask = ~isnan(pFlat);
    pAdj      = nan(size(pFlat));
    if any(validMask)
        [pAdj(validMask), ~, ~, ~] = fdr_BH(pFlat(validMask), 0.05, false);
    end
    pValPairwiseAdj = reshape(pAdj, nPairs, nNeurons);   % [nPairs × nNeurons]
else
    pValPairwiseAdj = pValPairwise;
end

% =========================================================================
% Store results
% =========================================================================
S.categoryName   = params.compareCategory;    % name of compared category
S.categoryLevels = levels;                    % actual level values
S.params         = params;

% Per-level results — field named by category and level value
for k = 1:nLevels
    % Build valid MATLAB field name from category + level value
    fName_k = sprintf('%s_%g', lower(strtrim(params.compareCategory)), levels(k));
    fName_k = strrep(fName_k, '.', 'p');    % replace decimal for valid field name
    fName_k = strrep(fName_k, '-', 'neg');  % replace negative sign

    S.(fName_k).pvalsResponse = pValPerLevel(k,:);    % [1 × nNeurons] p-values vs baseline
    S.(fName_k).ZScoreU       = zPerLevel(k,:);        % [1 × nNeurons] bias-corrected z-score
    S.(fName_k).ObsStat       = obsStatLevels(k,:);    % [1 × nNeurons] mean Diff (spikes/ms)
    S.(fName_k).nTrials       = size(allDiff{k}, 1);  % number of trials for this level
end

% Omnibus test
S.omnibus.pVal  = pValOmnibus;   % [1 × nNeurons] any difference across levels?
S.omnibus.F_obs = F_obs;          % [1 × nNeurons] observed F-statistic

% Pairwise results
for pr = 1:nPairs
    i        = pairIdx(pr,1);
    j        = pairIdx(pr,2);
    pairName = sprintf('%s_%g_vs_%g', ...
        lower(strtrim(params.compareCategory)), levels(i), levels(j));
    pairName = strrep(pairName, '.', 'p');
    pairName = strrep(pairName, '-', 'neg');

    S.pairwise.(pairName).pVal    = pValPairwise(pr,:);     % [1 × nNeurons] raw
    S.pairwise.(pairName).pValAdj = pValPairwiseAdj(pr,:);  % [1 × nNeurons] FDR corrected
    S.pairwise.(pairName).obsDiff = obsStatPairwise(pr,:);  % [1 × nNeurons] level_i minus level_j
end

fprintf('Saving results to %s\n', outputFile);
save(outputFile, '-struct', 'S');
results = S;

end % end main function


% =========================================================================
%% Local helper: permutation-compatible one-way ANOVA F-statistic
% =========================================================================
function F = computeFstat(Diff, levelLabels, levels)
% computeFstat - One-way ANOVA F-statistic for permutation testing.
%
% Inputs:
%   Diff        : [nTrials × nNeurons] response-minus-baseline values
%   levelLabels : [nTrials × 1] integer group membership per trial
%   levels      : unique level values [nLevels × 1]
%
% Output:
%   F           : [1 × nNeurons] F-statistic per neuron

    nLevels  = numel(levels);
    nTotal   = size(Diff, 1);
    nNeurons = size(Diff, 2);

    grandMean  = mean(Diff, 1);   % [1 × nNeurons]
    SS_between = zeros(1, nNeurons);
    SS_within  = zeros(1, nNeurons);

    for k = 1:nLevels
        mask      = levelLabels == levels(k);
        nk        = sum(mask);
        groupMean = mean(Diff(mask,:), 1);                          % [1 × nNeurons]
        SS_between = SS_between + nk * (groupMean - grandMean).^2;  % weighted group deviation
        SS_within  = SS_within + sum((Diff(mask,:) - groupMean).^2, 1);  % within-group variance
    end

    df_between = nLevels - 1;
    df_within  = nTotal - nLevels;

    F              = (SS_between / df_between) ./ (SS_within / df_within);
    F(SS_within==0) = 0;   % degenerate: all trials identical within groups
end