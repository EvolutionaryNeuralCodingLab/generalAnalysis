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
%       stimTypes       — which stimulus analyses to include (abbreviations)
%       splitBy         — category variable to split within each stim type
%                         (e.g. "size", "direction").  "" = no split.
%                         Experiments with <2 levels are automatically skipped.
%       splitLevels     — numeric vector of specific levels to use (e.g. [5 10 20]).
%                         Empty = use all available levels.  Experiments
%                         missing any of the requested levels are skipped.
%       binWidth        — PSTH bin width in ms
%       smooth          — Gaussian smoothing window in ms (0 = none)
%       TakeTopPercentTrials — fraction of trials to keep (1 = all trials)
%       byDepth         — stratify neurons by cortical depth
%
%   See the 'arguments' block below for the full parameter list and defaults.

% -------------------------------------------------------------------------
% Input validation via MATLAB arguments block
% -------------------------------------------------------------------------
arguments
    exList  double                                                         % vector of experiment IDs to include
    params.stimTypes  (1,:) string  = ["RG", "MB"]                         % stimulus types (abbreviations: MB, MBR, RG, SDGm, SDGs, NV, NI, FFF)
    params.splitBy    string        = ""                                   % category variable for within-stim split; "" = no split
    params.splitLevels double       = []                                   % specific levels to compare (e.g. [5 10 20]); empty = all available
    params.binWidth   double        = 10                                   % PSTH time-bin width in ms
    params.smooth     double        = 0                                    % Gaussian smoothing window in ms (0 = no smoothing)
    params.statType   string        = "maxPermuteTest"                     % which statistical test for p-values
    params.speed      string        = "max"                                % which speed condition for MB/MBR
    params.alpha      double        = 0.05                                 % significance threshold for neuron responsiveness
    params.shadeSTD   logical       = true                                 % shade ±SEM around the mean PSTH line
    params.postStim   double        = 500                                  % duration after stimulus onset to include (ms)
    params.preBase    double        = 200                                  % pre-stimulus baseline duration (ms)
    params.overwrite  logical       = false                                % if true, recompute even when a saved file exists
    params.overwriteResponse logical = false                               % if true, force recompute of ResponseWindow
    params.overwriteStats logical  = false                                 % if true, force recompute of statistics
    params.useCategoryPvals logical = false                                % if true and splitBy is active, use per-level p-values from StatisticsPerNeuronPerCategory (OR across levels) instead of general per-neuron p-values
    params.nBootCategory double    = 10000                                 % number of bootstrap iterations for StatisticsPerNeuronPerCategory
    params.TakeTopPercentTrials double = 1                                 % fraction of trials to keep (1 = all; see note below)
    params.zScore     logical       = false                                % z-score each neuron's PSTH to its own baseline
    params.PaperFig   logical       = false                                % export publication-quality figure via printFig
    params.byDepth    logical       = false                                % split neurons into 3 depth bins
    params.unionResponsive logical = false                                 % if true, include neurons responsive to ANY stimType (OR union across stim types)
end

% -------------------------------------------------------------------------
% NOTE ON TakeTopPercentTrials (default = 1 = all trials)
% -------------------------------------------------------------------------
% Selecting the top N% of trials by mean spike count inflates PSTH
% amplitudes and biases the response profile.  For a publication figure
% this is hard to justify unless there is a principled reason (e.g.
% attention gating in a behaving animal).  Default is 1 (all trials).
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Guard: splitBy and byDepth together create too many lines
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
    params.unionResponsive = false;
end

% -------------------------------------------------------------------------
% Guard: unionResponsive + useCategoryPvals is ambiguous
% -------------------------------------------------------------------------
if params.unionResponsive && params.useCategoryPvals
    error(['unionResponsive and useCategoryPvals cannot both be true. ', ...
           'The union pre-pass uses general per-neuron p-values (statType). ', ...
           'Use one mode at a time.']);
end

% -------------------------------------------------------------------------
% Load depth-bin info if byDepth is requested
% -------------------------------------------------------------------------
if params.byDepth
    depthFile = 'W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\NeuronDepths.mat';
    if ~exist(depthFile, 'file')
        error('NeuronDepths.mat not found. Run getNeuronDepths() first.');
    end
    D             = load(depthFile);
    depthTable    = D.depthTable;
    depthBinEdges = D.depthBinEdges;
    nDepthBins    = 3;
    fprintf('Depth bins loaded:\n');
    fprintf('  Bin 1 (shallow): %.0f – %.0f µm\n', depthBinEdges(1), depthBinEdges(2));
    fprintf('  Bin 2 (middle) : %.0f – %.0f µm\n', depthBinEdges(2), depthBinEdges(3));
    fprintf('  Bin 3 (deep)   : %.0f – %.0f µm\n', depthBinEdges(3), depthBinEdges(4));
else
    nDepthBins = 1;
end

% -------------------------------------------------------------------------
% Build save directory path
% -------------------------------------------------------------------------
NP_first  = loadNPclassFromTable(exList(1));
vs_first  = linearlyMovingBallAnalysis(NP_first);

basePath  = extractBefore(vs_first.getAnalysisFileName, 'lizards');
basePath  = [basePath 'lizards'];
saveDir   = fullfile(basePath, 'Combined_lizard_analysis');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% ---- Construct the filename ----
stimLabel    = strjoin(params.stimTypes, '-');
depthSuffix  = '';
if params.byDepth;  depthSuffix = '_byDepth'; end
splitSuffix  = '';
if params.splitBy ~= ""; splitSuffix = ['_by' char(params.splitBy)]; end
if ~isempty(params.splitLevels)
    lvlStr = strjoin(arrayfun(@(v) sprintf('%g',v), params.splitLevels, ...
        'UniformOutput', false), '_');
    splitSuffix = [splitSuffix '_lvl' lvlStr];
end
nameOfFile = sprintf('Ex_%d-%d_Combined_PSTHs_%s%s%s.mat', ...
    exList(1), exList(end), stimLabel, splitSuffix, depthSuffix);
fullSavePath = fullfile(saveDir, nameOfFile);

% -------------------------------------------------------------------------
% Decide: recompute or load from disk
% -------------------------------------------------------------------------
forloop = true;
if exist(fullSavePath, 'file') == 2 && ~params.overwrite
    S = load(fullSavePath);
    if isequal(S.expList, exList)
        fprintf('Loading saved PSTHs from:\n  %s\n', fullSavePath);
        forloop = false;
    else
        fprintf('Experiment list mismatch — recomputing.\n');
    end
end

% =========================================================================
%  MAIN COMPUTATION LOOP (skip if loaded from disk)
% =========================================================================
if forloop

    nStim = numel(params.stimTypes);
    nExp  = numel(exList);

    % =====================================================================
    %  DISCOVERY PASS — find valid sessions and category levels
    % =====================================================================
    % For each experiment × stim type, we try sessions [0, 1, 2] to find
    % one whose category column (splitBy) has ≥2 levels (or contains all
    % of the user-specified splitLevels).  Results are stored in sessionMap
    % so the main loop can skip invalid experiments without re-searching.
    %
    %   sessionMap(ei, s)  = session to use (-1 = skip this exp×stim)
    %   catLabelsAll{s}    = string array of category labels for stim s
    % =====================================================================

    sessionMap   = zeros(nExp, nStim);                                     % will hold session numbers; -1 = skip
    catLabelsAll = cell(nStim, 1);

    if params.splitBy ~= ""
        fprintf('\nDiscovering categories (splitBy = "%s") ...\n', params.splitBy);
    end

    for s = 1:nStim
        stimKey = params.stimTypes(s);

        if params.splitBy == ""
            % ---- No split: all experiments use default session (0) ------
            catLabelsAll{s}  = "all";
            sessionMap(:, s) = 0;
        else
            % ---- Split requested: scan experiments for valid sessions ---
            allLevelsFound = [];                                           % accumulate levels across experiments

            for ei = 1:nExp
                try
                    NP_tmp = loadNPclassFromTable(exList(ei));
                catch
                    sessionMap(ei, s) = -1;
                    continue
                end

                % Try sessions [0, 1, 2]; return first with ≥2 levels
                [~, sess, levels] = findValidSession( ...
                    NP_tmp, stimKey, params.speed, params.splitBy, ...
                    params.splitLevels, params.overwriteResponse);

                if sess < 0
                    sessionMap(ei, s) = -1;
                    fprintf('  Exp %d [%s]: no session with ≥2 levels of "%s" — will skip.\n', ...
                        exList(ei), stimKey, params.splitBy);
                else
                    sessionMap(ei, s) = sess;
                    allLevelsFound = [allLevelsFound; levels(:)]; %#ok<AGROW>
                    fprintf('  Exp %d [%s]: session %d has levels [%s]\n', ...
                        exList(ei), stimKey, sess, num2str(levels(:)', '%g '));
                end
            end

            % Determine the final global set of category levels
            uniqueVals = unique(allLevelsFound);

            % If the user requested specific levels, intersect
            if ~isempty(params.splitLevels)
                uniqueVals = intersect(uniqueVals, params.splitLevels(:));
            end

            if numel(uniqueVals) < 2
                fprintf('  [%s] splitBy="%s" has only %d global level — falling back to unsplit.\n', ...
                    stimKey, params.splitBy, numel(uniqueVals));
                catLabelsAll{s}  = "all";
                sessionMap(:, s) = 0;                                      % reset all to default session
            else
                catLabelsAll{s} = string(uniqueVals(:)');
                fprintf('  [%s] final category levels: %s\n', ...
                    stimKey, strjoin(catLabelsAll{s}, ', '));
            end
        end
    end

    % ----- Find the maximum number of categories across stim types -------
    maxCats = max(cellfun(@numel, catLabelsAll));

    % ----- Pre-allocate psthAll ------------------------------------------
    psthAll = cell(nStim, nDepthBins, maxCats);

    % ----- Time-axis parameters (locked on first successful experiment) --
    lockedPreBase = [];
    lockedNBins   = [];
    lockedEdges   = [];

    % =====================================================================
    %  MAIN EXPERIMENT LOOP
    % =====================================================================
    for ei = 1:nExp

        ex = exList(ei);
        fprintf('\n=== Experiment %d (%d/%d) ===\n', ex, ei, nExp);

        % ---- Load experiment --------------------------------------------
        try
            NP = loadNPclassFromTable(ex);
        catch ME
            warning('Could not load experiment %d: %s', ex, ME.message);
            appendNaNRow(psthAll, nStim, nDepthBins, catLabelsAll, lockedNBins);
            continue
        end

        % =============================================================
        %  Union-responsive pre-pass: for each stim type, identify
        %  responsive neurons using per-neuron p-values, then take the
        %  OR union.  The resulting index vector (unionENeurons) is
        %  applied to ALL stim types in the main loop below.
        %
        %  Rationale: good-sorted units are identical across stim types
        %  within a recording (spike sorting is stim-agnostic), so
        %  neuron indices are directly comparable across stim types.
        % =============================================================
        if params.unionResponsive
            eNeuronsPerStim = cell(1, nStim);                              % one index vector per stim type

            for su = 1:nStim
                stimTypeU = params.stimTypes(su);                          % current stim abbreviation
                sessU     = sessionMap(ei, su);                            % pre-computed session for this exp×stim

                if sessU < 0                                               % no valid session → no neurons from this stim
                    eNeuronsPerStim{su} = [];
                    continue
                end

                % ---- Build analysis object for this stim type -----------
                try
                    objU = buildStimObject(NP, stimTypeU, sessU);
                catch
                    eNeuronsPerStim{su} = [];
                    continue
                end
                if isempty(objU)
                    eNeuronsPerStim{su} = [];
                    continue
                end

                % ---- Check stimulus was actually presented --------------
                try
                    stimMissingU = isempty(objU.VST);
                catch
                    stimMissingU = true;
                end
                if stimMissingU
                    eNeuronsPerStim{su} = [];
                    continue
                end

                % ---- Ensure ResponseWindow is available -----------------
                try
                    objU.ResponseWindow('overwrite', params.overwriteResponse);
                catch
                    eNeuronsPerStim{su} = [];
                    continue
                end

                % ---- Load statistics ------------------------------------
                try
                    if params.statType == "BootstrapPerNeuron"
                        StatsU = objU.BootstrapPerNeuron;
                    elseif params.statType == "maxPermuteTest"
                        StatsU = objU.StatisticsPerNeuron;
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
                        pvalsU = StatsU.pvalsResponse;                     % flat struct fallback
                    catch
                        eNeuronsPerStim{su} = [];
                        continue
                    end
                end

                eNeuronsPerStim{su} = find(pvalsU < params.alpha);         % indices into goodU for this stim type
            end

            % ---- Take the union across all stim types -------------------
            unionENeurons = [];                                            % will hold sorted unique indices
            for su = 1:nStim
                unionENeurons = union(unionENeurons, eNeuronsPerStim{su}); % union returns sorted, unique
            end

            fprintf('  Union-responsive: %d neuron(s) responsive to ≥1 of [%s] in exp %d.\n', ...
                numel(unionENeurons), strjoin(params.stimTypes, ', '), ex);

            % ---- If the union is empty, skip the entire experiment ------
            if isempty(unionENeurons)
                fprintf('  No neurons responsive to any stim type — skipping exp %d.\n', ex);
                appendNaNRow(psthAll, nStim, nDepthBins, catLabelsAll, lockedNBins);
                continue                                                   % skip to next experiment
            end
        end

        % =================================================================
        %  Loop over stimulus types
        % =================================================================
        for s = 1:nStim

            stimType = params.stimTypes(s);

            % ---- Check session map: skip if no valid session found ------
            sess = sessionMap(ei, s);
            if sess < 0
                fprintf('  [%s] Skipping exp %d (no valid session).\n', stimType, ex);
                appendNaNRowForStim(psthAll, s, nDepthBins, catLabelsAll{s}, lockedNBins);
                continue
            end

            % ---- Construct the analysis object using the chosen session -
            try
                obj = buildStimObject(NP, stimType, sess);
            catch ME
                warning('Could not build %s (session %d) for exp %d: %s', ...
                    stimType, sess, ex, ME.message);
                appendNaNRowForStim(psthAll, s, nDepthBins, catLabelsAll{s}, lockedNBins);
                continue
            end

            if isempty(obj)
                appendNaNRowForStim(psthAll, s, nDepthBins, catLabelsAll{s}, lockedNBins);
                continue
            end

            % ---- Check that the stimulus was actually presented ---------
            % The constructor may succeed but VST is empty when the
            % stimulus protocol was not run in this experiment.
            try
                stimMissing = isempty(obj.VST);
            catch
                stimMissing = true;
            end
            if stimMissing
                fprintf('  [%s] Stimulus not present in exp %d — skipping.\n', stimType, ex);
                appendNaNRowForStim(psthAll, s, nDepthBins, catLabelsAll{s}, lockedNBins);
                continue
            end

            % ---- Ensure ResponseWindow is computed ----------------------
            try
                obj.ResponseWindow('overwrite', params.overwriteResponse);
            catch ME
                warning('  [%s] ResponseWindow failed for exp %d: %s — skipping.', ...
                    stimType, ex, ME.message);
                appendNaNRowForStim(psthAll, s, nDepthBins, catLabelsAll{s}, lockedNBins);
                continue
            end

            % ---- Load statistics (p-values per neuron) ------------------
            try
                if params.statType == "BootstrapPerNeuron"
                    Stats = obj.BootstrapPerNeuron;
                elseif params.statType == "maxPermuteTest"
                    Stats = obj.StatisticsPerNeuron;
                else
                    Stats = obj.ShufflingAnalysis;
                end
            catch ME
                warning('  [%s] Statistics failed for exp %d: %s — skipping.', ...
                    stimType, ex, ME.message);
                appendNaNRowForStim(psthAll, s, nDepthBins, catLabelsAll{s}, lockedNBins);
                continue
            end

            % ---- Determine field name and stim-onset offset -------------
            [fieldName, startStim] = getFieldAndOffset(obj, stimType, params.speed);

            % ---- Get sorted good units ----------------------------------
            p_sort = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
            label  = string(p_sort.label');
            goodU  = p_sort.ic(:, label == "good");

            % ---- Extract condition matrix C and stimulus onset times ----
            C = getConditionMatrix(obj, stimType, params.speed);
            directimesSorted = C(:, 1)' + startStim;

            % ---- Lock the time-axis on the first valid experiment -------
            preBase     = params.preBase;
            windowTotal = preBase + params.postStim;

            if isempty(lockedPreBase)
                lockedPreBase = preBase;
                lockedEdges   = 0 : params.binWidth : windowTotal;
                lockedNBins   = numel(lockedEdges) - 1;
                fprintf('  Locked window: preBase=%d ms, postStim=%d ms, nBins=%d\n', ...
                    lockedPreBase, params.postStim, lockedNBins);
            end

            % ---- Determine whether category split is active for this stim ---
            nCats         = numel(catLabelsAll{s});
            isSplitActive = params.splitBy ~= "" && ~isequal(catLabelsAll{s}, "all");

            % ---- Select responsive neurons ------------------------------
            % Two modes:
            %   (a) useCategoryPvals + active split: call
            %       StatisticsPerNeuronPerCategory and OR the per-level
            %       p-values.  A neuron is "responsive" if it passes the
            %       threshold for ANY of the compared levels.
            %   (b) Default: use the general per-neuron p-values from the
            %       overall statistics (statType).

            if isSplitActive && params.useCategoryPvals
                % ---- Category-level p-values ----------------------------
                try
                    catStats = obj.StatisticsPerNeuronPerCategory( ...
                        'compareCategory', char(params.splitBy), ...
                        'nBoot',           params.nBootCategory, ...
                        'overwrite',       params.overwriteStats);
                catch ME
                    warning('  [%s] StatisticsPerNeuronPerCategory failed for exp %d: %s — skipping.', ...
                        stimType, ex, ME.message);
                    appendNaNRowForStim(psthAll, s, nDepthBins, catLabelsAll{s}, lockedNBins);
                    continue
                end

                % OR across all category levels
                orMask = false;                                            % will broadcast to vector on first OR
                for ci = 1:nCats
                    levelVal = str2double(catLabelsAll{s}(ci));             % numeric level value
                    fName    = levelToFieldName(params.splitBy, levelVal);  % e.g. 'size_5'
                    if isfield(catStats, fName)
                        orMask = orMask | (catStats.(fName).pvalsResponse(:) < params.alpha);
                    end
                end
                eNeurons = find(orMask);

                fprintf('  [%s] Using per-category p-values: %d responsive neuron(s) in exp %d.\n', ...
                    stimType, numel(eNeurons), ex);
            else
                % ---- General per-neuron p-values ------------------------
                try
                    pvals = Stats.(fieldName).pvalsResponse;
                catch
                    pvals = Stats.pvalsResponse;
                end
                eNeurons = find(pvals < params.alpha);
            end

            % ---- Override with union set if requested -------------------
            if params.unionResponsive
                eNeurons = unionENeurons;                                  % replace per-stim selection with union across all stim types
                fprintf('  [%s] Using union-responsive set: %d neuron(s) in exp %d.\n', ...
                    stimType, numel(eNeurons), ex);
            end

            if isempty(eNeurons)
                fprintf('  [%s] No responsive neurons in exp %d.\n', stimType, ex);
                appendNaNRowForStim(psthAll, s, nDepthBins, catLabelsAll{s}, lockedNBins);
                continue
            end

            if ~(isSplitActive && params.useCategoryPvals)
                fprintf('  [%s] %d responsive neuron(s) in exp %d.\n', ...
                    stimType, numel(eNeurons), ex);
            end

            % ---- Extract per-trial category values (if splitting) -------
            catValues = [];
            if isSplitActive
                catCol    = getCategoryColumn(obj, stimType, params.speed, params.splitBy);
                catValues = catCol(:)';
            end

            % ==============================================================
            %  Loop over responsive neurons
            % ==============================================================
            psthRateNeurons = NaN(numel(eNeurons), lockedNBins, nCats);
            neuronBinIdx    = zeros(numel(eNeurons), 1);

            for ni = 1:numel(eNeurons)
                u = eNeurons(ni);

                % ---- Assign depth bin -----------------------------------
                if params.byDepth
                    depthRow = depthTable.Experiment == ex & depthTable.Unit == u;
                    if ~any(depthRow)
                        neuronBinIdx(ni) = 0;
                        continue
                    end
                    unitDepth = depthTable.Depth_um(depthRow);
                    if unitDepth <= depthBinEdges(2)
                        neuronBinIdx(ni) = 1;
                    elseif unitDepth <= depthBinEdges(3)
                        neuronBinIdx(ni) = 2;
                    else
                        neuronBinIdx(ni) = 3;
                    end
                else
                    neuronBinIdx(ni) = 1;
                end

                % ---- Build PSTH for each category -----------------------
                for ci = 1:nCats

                    if ~isSplitActive
                        trialMask = true(size(directimesSorted));
                    else
                        trialMask = abs(catValues - str2double(catLabelsAll{s}(ci))) < 1e-6;
                    end
                    catOnsets = directimesSorted(trialMask);

                    if isempty(catOnsets)
                        psthRateNeurons(ni, :, ci) = NaN(1, lockedNBins);
                        continue
                    end

                    % Build binary spike matrix: [nTrials × windowTotal_ms]
                    MRhist = BuildBurstMatrix( ...
                        goodU(:, u), ...
                        round(p_sort.t), ...
                        round(catOnsets - lockedPreBase), ...
                        round(windowTotal));
                    MRhist = squeeze(MRhist);

                    % ---- Optional: keep only top-N% of trials -----------
                    if ~isempty(params.TakeTopPercentTrials) && params.TakeTopPercentTrials < 1
                        MeanTrial  = mean(MRhist, 2);
                        [~, ind]   = sort(MeanTrial, 'descend');
                        nKeep      = max(1, round(numel(MeanTrial) * params.TakeTopPercentTrials));
                        MRhist     = MRhist(ind(1:nKeep), :);
                    end

                    nTrials = size(MRhist, 1);

                    % ---- Compute PSTH by direct bin-summation -----------
                    counts = zeros(1, lockedNBins);
                    for bi = 1:lockedNBins
                        msStart = lockedEdges(bi) + 1;
                        msEnd   = lockedEdges(bi + 1);
                        counts(bi) = sum(MRhist(:, msStart:msEnd), 'all');
                    end
                    psthRateNeurons(ni, :, ci) = (counts / nTrials) / (params.binWidth / 1000);

                end % category loop
            end % neuron loop

            % ==============================================================
            %  z-score each neuron individually (if requested)
            % ==============================================================
            tAxis = lockedEdges(1:end-1);
            if params.zScore
                baselineBins = tAxis < lockedPreBase;
                for ni = 1:size(psthRateNeurons, 1)
                    for ci = 1:nCats
                        trace = psthRateNeurons(ni, :, ci);
                        if all(isnan(trace)); continue; end
                        bMean = mean(trace(baselineBins), 'omitnan');
                        bStd  = std(trace(baselineBins), 0, 'omitnan');
                        if bStd > 0
                            psthRateNeurons(ni, :, ci) = (trace - bMean) / bStd;
                        else
                            psthRateNeurons(ni, :, ci) = NaN;
                        end
                    end
                end
            end

            % ==============================================================
            %  Average across neurons per depth bin × category, append
            % ==============================================================
            for b = 1:nDepthBins
                binNeurons = neuronBinIdx == b;

                if ~any(binNeurons)
                    for ci = 1:nCats
                        psthAll{s, b, ci} = appendOrInit(psthAll{s, b, ci}, NaN(1, lockedNBins));
                    end
                    continue
                end

                for ci = 1:nCats
                    catData = psthRateNeurons(binNeurons, :, ci);
                    if all(isnan(catData), 'all')
                        psthAll{s, b, ci} = appendOrInit(psthAll{s, b, ci}, NaN(1, lockedNBins));
                    else
                        psthExp = mean(catData, 1, 'omitnan');
                        psthAll{s, b, ci} = appendOrInit(psthAll{s, b, ci}, psthExp(:)');
                    end
                end

                fprintf('    [%s] Depth bin %d: %d neuron(s) in exp %d.\n', ...
                    stimType, b, sum(binNeurons), ex);
            end

        end % stim-type loop
    end % experiment loop

    % =====================================================================
    %  Save results to disk
    % =====================================================================
    S.expList       = exList;
    S.lockedEdges   = lockedEdges;
    S.lockedPreBase = lockedPreBase;
    S.params        = params;
    S.catLabelsAll  = catLabelsAll;

    for s = 1:numel(params.stimTypes)
        stimField = matlab.lang.makeValidName(params.stimTypes(s));
        for b = 1:nDepthBins
            for ci = 1:numel(catLabelsAll{s})
                fieldKey = sprintf('%s_bin%d_cat%d', stimField, b, ci);
                S.(fieldKey) = psthAll{s, b, ci};
            end
        end
    end

    save(fullSavePath, '-struct', 'S');
    fprintf('\nSaved PSTHs to:\n  %s\n', fullSavePath);

else
    % =================================================================
    %  Load psthAll from the saved struct S
    % =================================================================
    lockedEdges   = S.lockedEdges;
    lockedPreBase = S.lockedPreBase;
    catLabelsAll  = S.catLabelsAll;

    maxCats = max(cellfun(@numel, catLabelsAll));
    psthAll = cell(numel(params.stimTypes), nDepthBins, maxCats);

    for s = 1:numel(params.stimTypes)
        stimField = matlab.lang.makeValidName(params.stimTypes(s));
        for b = 1:nDepthBins
            for ci = 1:numel(catLabelsAll{s})
                fieldKey = sprintf('%s_bin%d_cat%d', stimField, b, ci);
                if isfield(S, fieldKey)
                    psthAll{s, b, ci} = S.(fieldKey);
                else
                    warning('Field "%s" not found in saved file.', fieldKey);
                    psthAll{s, b, ci} = [];
                end
            end
        end
    end
end

% =========================================================================
%  PLOTTING
% =========================================================================

tAxis     = lockedEdges(1:end-1);
tAxisPlot = tAxis - lockedPreBase;

% ---- Colour palette and label maps -------------------------------------
nStim      = numel(params.stimTypes);
baseColors = lines(nStim);

stimLegendMap = containers.Map( ...
    {'MB', 'MBR', 'RG', 'SDGm', 'SDGs', 'NV', 'NI', 'FFF'}, ...
    {'MB', 'MBR', 'RG', 'SDGm', 'SDGs', 'NV', 'NI', 'FFF'});

depthShades = [0.05, 0.45, 0.78];
binLabels   = {'shallow', 'middle', 'deep'};

% ---- First pass: smooth traces, compute mean & SEM, find global ylim ---
yMax = -Inf;
yMin =  Inf;

meanStore = cell(nStim, nDepthBins, max(cellfun(@numel, catLabelsAll)));
semStore  = cell(nStim, nDepthBins, max(cellfun(@numel, catLabelsAll)));
nExpStore = zeros(nStim, nDepthBins, max(cellfun(@numel, catLabelsAll)));

for s = 1:nStim
    for b = 1:nDepthBins
        for ci = 1:numel(catLabelsAll{s})
            data = psthAll{s, b, ci};
            if isempty(data); continue; end
            validRows = ~all(isnan(data), 2);
            data      = data(validRows, :);
            if isempty(data); continue; end

            nValid = size(data, 1);
            nExpStore(s, b, ci) = nValid;

            % Smooth each experiment's trace FIRST, then compute stats
            if params.smooth > 0
                smoothBins = round(params.smooth / params.binWidth);
                for ri = 1:nValid
                    data(ri, :) = smoothdata(data(ri, :), 'gaussian', smoothBins);
                end
            end

            meanTrace = mean(data, 1, 'omitnan');
            semTrace  = std(data, 0, 1, 'omitnan') / sqrt(nValid);

            meanStore{s, b, ci} = meanTrace;
            semStore{s, b, ci}  = semTrace;

            yMax = max(yMax, max(meanTrace + semTrace));
            yMin = min(yMin, min(meanTrace - semTrace));
        end
    end
end

yPad = (yMax - yMin) * 0.1;
if params.zScore
    yLims = [yMin - yPad, yMax + yPad];
else
    yLims = [max(0, yMin - yPad), yMax + yPad];
end

% ---- Create figure ------------------------------------------------------
fig = figure;
set(fig, 'Units', 'centimeters', 'Position', [5 5 10 7]);
ax  = axes(fig);
hold(ax, 'on');

legendHandles = [];
legendLabels  = {};

% ---- Category colour maps -----------------------------------------------
catColorMaps = cell(nStim, 1);
for s = 1:nStim
    nc = numel(catLabelsAll{s});
    if nc == 1
        catColorMaps{s} = baseColors(s, :);
    else
        cmap = zeros(nc, 3);
        for ci = 1:nc
            frac       = (ci - 1) / max(nc - 1, 1);
            cmap(ci,:) = baseColors(s,:) .* (1 - 0.5*frac) + [0.3 0.1 0]*frac;
            cmap(ci,:) = min(max(cmap(ci,:), 0), 1);
        end
        catColorMaps{s} = cmap;
    end
end

% ---- Plot each condition ------------------------------------------------
for s = 1:nStim

    stimKey = char(params.stimTypes(s));
    if isKey(stimLegendMap, stimKey)
        shortName = stimLegendMap(stimKey);
    else
        shortName = stimKey;
    end

    for b = 1:nDepthBins
        for ci = 1:numel(catLabelsAll{s})

            meanPSTH = meanStore{s, b, ci};
            semPSTH  = semStore{s, b, ci};
            if isempty(meanPSTH); continue; end

            nValid = nExpStore(s, b, ci);

            % ---- Determine line colour and legend label -----------------
            isSplitHere = params.splitBy ~= "" && ~isequal(catLabelsAll{s}, "all");
            if params.byDepth
                lineColor   = baseColors(s,:) * (1 - depthShades(b));
                legendLabel = sprintf('%s %s (%.0f–%.0f µm, n=%d)', ...
                    shortName, binLabels{b}, ...
                    depthBinEdges(b), depthBinEdges(b+1), nValid);
            elseif isSplitHere
                lineColor   = catColorMaps{s}(ci, :);
                legendLabel = sprintf('%s %s=%s (n=%d)', ...
                    shortName, params.splitBy, catLabelsAll{s}(ci), nValid);
            else
                lineColor   = baseColors(s,:);
                legendLabel = sprintf('%s (n=%d)', shortName, nValid);
            end

            % ---- SEM shading -------------------------------------------
            if params.shadeSTD && nValid > 1
                upper = meanPSTH + semPSTH;
                lower = meanPSTH - semPSTH;
                xFill = [tAxisPlot(:)', fliplr(tAxisPlot(:)')];
                yFill = [upper(:)',      fliplr(lower(:)')];
                fill(ax, xFill, yFill, lineColor, ...
                    'FaceAlpha', 0.08, 'EdgeColor', 'none');
            end

            % ---- Mean PSTH line ----------------------------------------
            h = plot(ax, tAxisPlot(:)', meanPSTH(:)', ...
                'Color', lineColor, 'LineWidth', 1.5);

            legendHandles(end+1) = h; %#ok<AGROW>
            legendLabels{end+1}  = legendLabel; %#ok<AGROW>

        end % category loop
    end % depth-bin loop
end % stim-type loop

% ---- Reference lines ----------------------------------------------------
xline(ax, 0,               'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(ax, params.postStim, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% ---- Axis labels and formatting -----------------------------------------
if params.zScore; yLabel = 'Z-score'; else; yLabel = 'Firing rate [spk/s]'; end

xlabel(ax, 'Time re. stim onset [ms]', 'FontName', 'Helvetica', 'FontSize', 8);
ylabel(ax, yLabel,                      'FontName', 'Helvetica', 'FontSize', 8);
xlim(ax, [tAxisPlot(1), tAxisPlot(end)]);
ylim(ax, yLims);

legend(legendHandles, legendLabels, 'Location', 'northwest', ...
    'FontName', 'Helvetica', 'FontSize', 7);

ax.FontName       = 'Helvetica';
ax.FontSize       = 8;
ax.YAxis.FontSize = 8;
ax.XAxis.FontSize = 8;
hold(ax, 'off');

title(ax, sprintf('N = %d experiments', numel(exList)), ...
    'FontName', 'Helvetica', 'FontSize', 10);



% ---- Export publication figure if requested -----------------------------
if params.PaperFig
    stimStr = strjoin(params.stimTypes, '-');
    vs_first.printFig(fig, sprintf('PSTH-%s%s%s', stimStr, splitSuffix, depthSuffix), ...
        PaperFig = params.PaperFig);
end

end  % end of main function


% #########################################################################
%                       LOCAL HELPER FUNCTIONS
% #########################################################################


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

obj          = [];
validSession = -1;
levels       = [];

for sess = [0, 1, 2]
    candidate = buildStimObject(NP, stimKey, sess);
    if isempty(candidate)
        continue
    end

    % Check that the stimulus was actually presented
    try
        if isempty(candidate.VST)
            continue
        end
    catch
        continue
    end

    % Ensure ResponseWindow is computed
    try
        candidate.ResponseWindow('overwrite', overwriteRW);
    catch
        continue
    end

    % Extract condition matrix and column names
    rw = candidate.ResponseWindow;
    [C, colNames] = getCmatrixLocal(rw, stimKey, speedParam);
    if isempty(C) || isempty(colNames)
        continue
    end

    % Find the category column by name
    catIdx = find(strcmpi(colNames, splitBy), 1);
    if isempty(catIdx)
        continue
    end

    % Column mapping: colNames(k) → C(:, k+1)  [C(:,1) = onset times]
    catColIdx   = catIdx + 1;
    rawCol      = C(:, catColIdx);
    rawCol      = rawCol(~isnan(rawCol));
    availLevels = uniquetol(rawCol, 1e-6);

    % If specific levels requested, verify they are all present
    if ~isempty(splitLevels)
        allPresent = true;
        for lv = splitLevels(:)'
            if ~any(abs(availLevels - lv) < 1e-6)
                allPresent = false;
                break
            end
        end
        if ~allPresent
            continue
        end
        useLevels = splitLevels(:);
    else
        useLevels = availLevels;
    end

    % Need ≥2 levels to be worth splitting
    if numel(useLevels) < 2
        continue
    end

    % Found a valid session
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
%   colNames = parameter names only (rw.colNames{1}(5:end)).
%   C(:,1) = onset times.  C(:, k+1) = colNames(k).

C        = [];
colNames = {};

% Try to read colNames (same regardless of speed field)
try
    allColNames = rw.colNames{1};
    colNames    = allColNames(5:end);
catch
    return
end

% Get C from the correct sub-field
switch stimKey
    case {"MB", "MBR"}
        if speedParam == "max"
            fld = 'Speed1';
        else
            fld = 'Speed2';
        end
        if isfield(rw, fld)
            C = rw.(fld).C;
        else
            speedFields = fieldnames(rw);
            speedFields = speedFields(startsWith(speedFields, 'Speed'));
            if ~isempty(speedFields)
                C = rw.(speedFields{end}).C;
            end
        end

    case "SDGm"
        if isfield(rw, 'C')
            C = rw.C;
        elseif isfield(rw, 'Moving') && isfield(rw.Moving, 'C')
            C = rw.Moving.C;
        end

    case "SDGs"
        if isfield(rw, 'C')
            C = rw.C;
        elseif isfield(rw, 'Static') && isfield(rw.Static, 'C')
            C = rw.Static.C;
        end

    otherwise
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
%   Returns [] if construction fails.

if nargin < 3; session = 0; end

obj = [];
try
    % SDGm and SDGs both use StaticDriftingGratingAnalysis
    switch stimKey
        case {"SDGm", "SDGs"},  ctorKey = "SDG";
        otherwise,              ctorKey = stimKey;
    end

    if session == 0
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
    obj = [];
end
end


function [fieldName, startStim] = getFieldAndOffset(obj, stimKey, speedParam)
% getFieldAndOffset  Return the response-table field name and the stimulus-
%   onset offset (ms) for the given stimulus abbreviation key.

startStim = 0;

switch stimKey
    case "MB"
        if speedParam == "max"; fieldName = 'Speed1'; else; fieldName = 'Speed2'; end
    case "MBR"
        if speedParam == "max"; fieldName = 'Speed1'; else; fieldName = 'Speed2'; end
    case "SDGs"
        fieldName = 'Static';
    case "SDGm"
        fieldName = 'Moving';
        startStim = obj.VST.static_time * 1000;
    otherwise
        fieldName = 'Speed1';
end
end


function C = getConditionMatrix(obj, stimType, speedParam)
% getConditionMatrix  Extract condition matrix C from ResponseWindow.

[fieldName, ~] = getFieldAndOffset(obj, stimType, speedParam);
NeuronResp     = obj.ResponseWindow;

try
    C = NeuronResp.(fieldName).C;
catch
    C = NeuronResp.C;
end
end


function catCol = getCategoryColumn(obj, stimType, speedParam, splitBy)
% getCategoryColumn  Extract the category column from C by matching
%   splitBy against column names in ResponseWindow.
%
%   Column layout:
%     colNames{1}(5:end) = stimulus-parameter names
%     C(:,1)   = onset times
%     C(:,2:)  = parameter columns
%     → paramNames(k) maps to C(:, k+1)

responseParams = obj.ResponseWindow;

allColNames = responseParams.colNames{1};
paramNames  = allColNames(5:end);

matchIdx = find(strcmpi(paramNames, splitBy), 1);
if isempty(matchIdx)
    error(['splitBy = "%s" does not match any column in colNames.\n' ...
           '  Available: %s'], splitBy, strjoin(string(paramNames), ', '));
end

colIdxInC = matchIdx + 1;                                                  % paramNames(1) → C(:,2)

C = getConditionMatrix(obj, stimType, speedParam);
catCol = C(:, colIdxInC);
end


function arr = appendOrInit(arr, newRow)
% appendOrInit  Append a row, or initialize if empty.
if isempty(arr)
    arr = newRow;
else
    arr = [arr; newRow];
end
end


function appendNaNRow(psthAll, nStim, nDepthBins, catLabelsAll, lockedNBins)
% appendNaNRow  Insert NaN placeholder for ALL stim × depth × cat conditions.
if isempty(lockedNBins); return; end
for s = 1:nStim
    for b = 1:nDepthBins
        for ci = 1:numel(catLabelsAll{s})
            psthAll{s, b, ci} = appendOrInit(psthAll{s, b, ci}, NaN(1, lockedNBins));
        end
    end
end
end


function appendNaNRowForStim(psthAll, s, nDepthBins, catLabels, lockedNBins)
% appendNaNRowForStim  Insert NaN placeholder for one stim type.
if isempty(lockedNBins); return; end
for b = 1:nDepthBins
    for ci = 1:numel(catLabels)
        psthAll{s, b, ci} = appendOrInit(psthAll{s, b, ci}, NaN(1, lockedNBins));
    end
end
end


function fName = levelToFieldName(catName, value)
% levelToFieldName  Build a valid field name matching the convention used by
%   StatisticsPerNeuronPerCategory.
%
%   Examples:
%     levelToFieldName("size", 5)     →  'size_5'
%     levelToFieldName("speed", 0.3)  →  'speed_0p3'
%     levelToFieldName("dir", -90)    →  'dir_neg90'

fName = sprintf('%s_%g', lower(strtrim(char(catName))), value);            % e.g. 'size_5' or 'speed_0.3'
fName = strrep(fName, '.', 'p');                                           % replace decimal point with 'p'
fName = strrep(fName, '-', 'neg');                                         % replace minus sign with 'neg'
end