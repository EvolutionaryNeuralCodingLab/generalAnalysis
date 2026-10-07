function [tempTable, summaryStats] = AllExpAnalysis(expList, params)
% AllExpAnalysis  Pool neural responses across Neuropixels recordings,
%   run pairwise statistical comparisons via hierarchical bootstrapping,
%   and generate publication-ready swarm and scatter plots.
%
%   Supports three modes:
%
%   MODE 1 — ACROSS-STIMULUS:
%     ComparePairs = {'SDGm','SDGs'} compares z-scores/spike rates between
%     different stimulus types.  Neurons significant for ANY stimulus in the
%     set are included (OR union).
%
%   MODE 2 — WITHIN-STIMULUS CATEGORY (all levels):
%     ComparePairs = {'MB'}, CompareCategory = "size" compares all levels of
%     a category within a single stimulus.  Uses StatisticsPerNeuronPerCategory
%     for per-level statistics.  Recordings without >=2 levels are skipped
%     (tries Session=1 then Session=2).
%
%   MODE 3 — SPECIFIC-LEVEL ACROSS-STIMULUS:
%     ComparePairs    = {'MB','SDGm'}
%     CompareCategory = {'direction','direction'}
%     CompareLevels   = {[0],[0]}
%     Compares specific level(s) of (possibly different) categories across
%     different stimuli.  Each stimulus has its own category and level list.
%     For a given stimulus, all requested levels must exist in a single
%     session of that stimulus, otherwise the experiment is skipped.
%
%   INPUTS
%     expList  (1,:) double   Row vector of experiment IDs from master Excel.
%     params   Name-value     See arguments block below.
%
%   OUTPUTS
%     tempTable  table   Fraction-responsive table (one row per insertion x
%                        item), filtered to insertions containing all
%                        compared items.
%     summaryStats struct  Central-tendency summaries behind every swarm panel.
%                        Fields: Zscore, SpkR, FracResponsive; each with
%                          .groups - one row per stimulus
%                          .diffs  - one row per pair (stimA, stimB, p)
%                        Columns: n, nAnimals, nInsertions, mean, median
%                        (plain pooled descriptives; `mean` is the plotted bar);
%                        bootMean, ciLo, ciHi (hierarchical bootstrap;
%                        ciLo/ciHi are the plotted uncertainty patch);
%                        effect, effectLo, effectHi (Cohen's d, d_z on the
%                        differences, with bootstrap CI); propPos, propPosLo,
%                        propPosHi (P(value > 0) with bootstrap CI - the
%                        common-language effect size on the differences).
%                        bootMean is diagnostic only - it is NOT an equal-
%                        weight-per-animal mean, and a large mean-vs-bootMean
%                        gap signals unbalanced neuron counts across animals.
%                        FracResponsive additionally carries .sensitivity when
%                        minSomaticN > 0: the same summary recomputed over
%                        insertions with at least that many somatic units, as a
%                        robustness check. The primary result never excludes.
%                        All blocks are printed to the console, and saved into
%                        the companion _analysisStruct.mat when PaperFig=true.
%                        See the plotSwarmBootstrapWithComparisons header.
%
%   EXAMPLES
%     % Mode 1
%     t = AllExpAnalysis([49:54 64:66], ComparePairs = {'SDGm','SDGs'});
%
%     % Mode 2
%     t = AllExpAnalysis([49:54 64:66], ...
%         ComparePairs = {'MB'}, CompareCategory = "size");
%
%     % Mode 3
%     t = AllExpAnalysis([49:54 64:66], ...
%         ComparePairs    = {'MB','SDGm'}, ...
%         CompareCategory = {'direction','direction'}, ...
%         CompareLevels   = {[0],[0]});
%
%   See also: hierBoot,  ,
%             StatisticsPerNeuronPerCategory

% =========================================================================
%                        ARGUMENTS BLOCK
% =========================================================================
arguments
    expList (1,:) double                            % Row vector of experiment IDs from master Excel

    % --- Mode selection ---
    params.ComparePairs     cell                    % Stimuli to compare (cell of char/string)
    params.CompareCategory                  = ""    % "" -> mode 1; string -> mode 2; cell -> mode 3
    params.CompareLevels    cell            = {}    % Cell of numeric vectors (non-empty -> mode 3)

    % --- Responsiveness filter ---
    params.useGeneralFilter logical         = false % Use general per-neuron p-values instead of per-level p
    params.threshold        double          = 0.05  % p-value cutoff for the responsiveness OR-mask

    % --- ResponseWindow parameters ---
    params.RespDurationWin  double          = 100   % Duration window (ms) for ResponseWindow computation
    params.overwriteResponse logical        = false % Force recomputation of ResponseWindow

    % --- StatisticsPerNeuron parameters (maxPermuteTest) ---
    params.BaseRespWindow   double          = 100   % Base response window (ms) for statistics computation
    params.SpatialGridMode  logical         = false % When true, forces BaseRespWindow = 200 ms
    params.maxCategory      logical         = false % Use max-category mode in StatisticsPerNeuron
    params.applyFDR         logical         = false % Apply FDR correction inside the statistics functions
    params.overwriteStats   logical         = false % Force recomputation of statistics
    params.CategoryMaximized                = ''    %Category to be maximized along levels of category to comapre
    params.durationWindow   double          = []   %Window length for moving window pval mode for moving ball
    % --- Bootstrap parameters ---
    params.nBoot            double          = 10000 % Iterations for pairwise hierarchical bootstrap
    params.nBootCategory    double          = 10000 % Iterations for per-category bootstrap
    params.nBootPval        double          = 1e5   % Iterations for the reported p-values only. Higher than
                                                    %   nBoot because a p near 0.001 is estimated from only
                                                    %   ~10 tail counts at 1e4, so Monte Carlo error is the
                                                    %   same order as the threshold. Cost scales with the
                                                    %   number of pairs: set to nBoot to restore old speed.
    params.pvalSeed         double          = 0     % RNG seed applied immediately before each p-value loop.
                                                    %   Without it the p-values inherit whatever global RNG
                                                    %   state the preceding plotting calls left behind, so
                                                    %   re-running the same analysis could change a reported p.
    params.minSomaticN      double          = 5     % SENSITIVITY analysis for the fraction panel ONLY.
                                                    %   An insertion with few somatic units yields a quantised
                                                    %   fraction (2 units -> only 0, 0.5 or 1), which inflates
                                                    %   the SD and attenuates the effect size, and the
                                                    %   bootstrap cannot down-weight it: the fraction analysis
                                                    %   hands hierBootMatchFreq exactly one value per insertion,
                                                    %   so every insertion carries equal weight regardless of
                                                    %   how many units it held.
                                                    %   The PRIMARY figure, p-value and annotation always use
                                                    %   EVERY insertion; this only adds a second, filtered
                                                    %   summary alongside them. 0 disables the extra work.
                                                    %   Never applied to the z-score/spike-rate panels: there
                                                    %   the observation is the neuron, and dropping sparse
                                                    %   insertions would remove marginally-significant units
                                                    %   and inflate the mean z-score.

    % --- Data extraction ---
    params.useZmean         logical         = true  % true = use z_mean field; false = peak spike rate
    params.useTtest         logical         = false % true = use ttest to obtain pvals

    % --- Post-hoc correction (applied AFTER extraction, distinct from applyFDR) ---
    params.useFDR           logical         = false % Benjamini-Hochberg FDR on extracted p-values

    params.markUnits       cell          = {}          % N×4 {NeurID,stimulus,animal,insertion}; row1->X, row2->+


    % --- Use specific categories to create minimum trial set ---
    params.CategoryLevels   cell    = {}    % per-stimulus factors kept SEPARATE (their crossing = levels).
                                            %   {} -> "all" (every factor separate) for every stimulus.
                                            %   1 entry -> broadcast to all; else one entry per ComparePairs stim.
                                            %   Each entry: string array of names, or "all". e.g. {["directions","offsets"],"all"}
    params.Lock             cell    = {}    % per-stimulus lock spec; each entry {name,value,...} or {} (no lock).
                                            %   Broadcast like CategoryLevels. e.g. {{'luminosities',255},{}}

    % --- Output ---
    params.overwrite        logical         = false % Force rerun of the entire per-experiment loop
    params.PaperFig         logical         = false % Save publication-quality figures via printFig


end

% =========================================================================
% SECTION 1 — DETECT MODE, VALIDATE, AND APPLY GLOBAL OVERRIDES
% =========================================================================

% --- SpatialGridMode override ---
% In SpatialGridMode the analysis window is fixed at 200 ms.  Apply the
% override here so that every downstream call sees the corrected value.
if params.SpatialGridMode
    params.BaseRespWindow = 300;    % override user-supplied value
end

% --- Detect operating mode from parameter combinations ---
if ~isempty(params.CompareLevels) || iscell(params.CompareCategory)
    % MODE 3: specific category levels compared across stimuli
    mode = 3;

    % CompareCategory must be a cell or string array with one entry per stimulus
    assert(iscell(params.CompareCategory) || isstring(params.CompareCategory), ...
        'Mode 3: CompareCategory must be a cell or string array, one per stimulus.');
    assert(numel(params.CompareCategory) == numel(params.ComparePairs), ...
        'Mode 3: CompareCategory must have same length as ComparePairs.');
    % Empty CompareLevels → auto-discover all levels per stimulus
    if isempty(params.CompareLevels)
        params.CompareLevels = repmat({[]}, 1, numel(params.ComparePairs));
    end
    assert(numel(params.CompareLevels) == numel(params.ComparePairs), ...
        'Mode 3: CompareLevels must have same length as ComparePairs.');

    

    % Normalise CompareCategory to a cell of char arrays for uniform handling
    catList = cell(1, numel(params.CompareCategory));
    for i = 1:numel(params.CompareCategory)
        if iscell(params.CompareCategory)
            catList{i} = char(strtrim(params.CompareCategory{i}));  % cell input
        else
            catList{i} = char(strtrim(params.CompareCategory(i)));  % string array input
        end
    end
    fprintf('=== Mode 3: specific-level cross-stimulus comparison ===\n');

elseif (ischar(params.CompareCategory) || isstring(params.CompareCategory)) && ...
       strtrim(string(params.CompareCategory)) ~= ""
    % MODE 2: all levels of one category within a single stimulus
    mode = 2;

    assert(numel(params.ComparePairs) == 1, ...
        'Mode 2: requires exactly one stimulus in ComparePairs.');

    stimName = params.ComparePairs{1};   % the single stimulus being decomposed
    catName  = char(strtrim(string(params.CompareCategory)));  % category column name
    fprintf('=== Mode 2: within-stimulus category "%s" in %s ===\n', catName, stimName);

else
    % MODE 1: direct across-stimulus comparison
    mode = 1;

    assert(numel(params.ComparePairs) >= 2, ...
        'Mode 1: requires >=2 stimuli in ComparePairs.');
end

% Boolean shortcuts used throughout the function
isCategoryMode      = (mode == 2);   % within-stimulus category comparison
isSpecificLevelMode = (mode == 3);   % specific (stim, cat, level) tuples

% Unique stimulus names that need to be loaded (deduplicated, preserving order)
stimsNeeded = unique(params.ComparePairs, 'stable');

% --- Focus specs: one entry per stimulus, aligned to stimsNeeded ---
nStim            = numel(stimsNeeded);                                    % number of stimuli
catLevelsPerStim = normalizeFocusSpec(params.CategoryLevels, nStim, "all");  % separate-factor lists
lockPerStim      = normalizeFocusSpec(params.Lock,           nStim, {});     % lock specs

% Focus is incompatible with the spatial-grid path (no catTable/t-test there).
for si = 1:nStim
    assert(~(params.SpatialGridMode && isFocusActive(catLevelsPerStim{si}, lockPerStim{si})), ...
        'CategoryLevels/Lock are not supported with SpatialGridMode=true.');
end

% Map each shared object key to its focus spec (first stimulus wins for shared SDG object).
specByKey = containers.Map();                                             % key -> struct(catLevels, lock)
for si = 1:nStim
    k = getObjKey(stimsNeeded{si});                                      % shared analysis-object key
    if ~specByKey.isKey(k)
        specByKey(k) = struct('catLevels', {catLevelsPerStim{si}}, 'lock', {lockPerStim{si}});
    end
end


% =========================================================================
% SECTION 2 — DIRECTORY SETUP AND CACHE MANAGEMENT
% =========================================================================

% Load the first experiment to extract directory paths for saving
NP0 = loadNPclassFromTable(expList(end));
vs0 = linearlyMovingBallAnalysis(NP0);

% Build the root path and combined-analysis save directory
rootPath = extractBefore(vs0.dataObj.recordingDir, 'lizards');  % path up to 'lizards'
rootPath = [rootPath 'lizards'];                                 % include 'lizards' folder
saveDir  = fullfile(rootPath, 'Combined_lizard_analysis');       % pooled data directory
if ~exist(saveDir, 'dir')
    mkdir(saveDir);  % create if it does not exist
end

% Construct a descriptive filename for the cached pooled data.
% Encoding all comparison parameters in the filename prevents accidental
% cache collisions when the same experiment list is analysed differently.
switch mode
    case 1
        % Mode 1: stimulus names joined by hyphens
        nameOfFile = sprintf('Ex_%d-%d_Combined_%s.mat', ...
            expList(1), expList(end), strjoin(stimsNeeded, '-'));
    case 2
        % Mode 2: stimulus + category name
        nameOfFile = sprintf('Ex_%d-%d_Combined_%s_%s.mat', ...
            expList(1), expList(end), stimName, lower(catName));
    case 3
        % Mode 3: each (stim, cat, levels) tuple encoded
        parts = cell(1, numel(stimsNeeded));
        for si = 1:numel(stimsNeeded)
            if isempty(catList{si})
                parts{si} = stimsNeeded{si};
            elseif isempty(params.CompareLevels{si})
                parts{si} = sprintf('%s-%s-all', stimsNeeded{si}, catList{si});
            else
                lvStr = strjoin(arrayfun(@(v) sprintf('%g', v), ...
                    params.CompareLevels{si}, 'UniformOutput', false), '_');
                parts{si} = sprintf('%s-%s-%s', stimsNeeded{si}, catList{si}, lvStr);
            end
        end
        nameOfFile = sprintf('Ex_%d-%d_SpecLvl_%s.mat', ...
            expList(1), expList(end), strjoin(parts, '__'));
end
savePath = fullfile(saveDir, nameOfFile);  % full path for the cache .mat

% Decide whether the per-experiment loop needs to run.
% Skip if the cache exists, matches the experiment list, and overwrite is off.
runLoop = true;
if exist(savePath, 'file') == 2 && ~params.overwrite
    S = load(savePath);                         % load cached struct
    if isfield(S, 'expList') && isequal(S.expList, expList)
        runLoop = false;                        % cache is valid — skip the loop
    end
end

% =========================================================================
% SECTION 3 — INITIALISE LONG-FORMAT TABLES
% =========================================================================

% TableStimComp: one row per (neuron x stimulus).  Holds z-scores and spike
% rates for every neuron that passes the responsiveness filter.
TableStimComp = table( ...
    categorical.empty(0,1), categorical.empty(0,1), categorical.empty(0,1), ...
    categorical.empty(0,1), categorical.empty(0,1), double.empty(0,1), double.empty(0,1), ...
   'VariableNames', {'animal','insertion','realInsertion','stimulus','NeurID','Z-score','SpkR'});

% TableRespNeurs: one row per (insertion x stimulus).  Counts of responsive
% neurons and total somatic neurons for fraction-responsive analysis.
TableRespNeurs = table( ...
    categorical.empty(0,1), categorical.empty(0,1), categorical.empty(0,1), ...
    double.empty(0,1), double.empty(0,1), ...
    'VariableNames', {'animal','insertion','stimulus','respNeur','totalSomaticN'});

% In mode 2, level labels are determined from the first valid recording.
% In mode 3, comparison labels are fixed by parameters from the start.
levelLabels     = {};    % mode 2: populated on first valid recording
fixedCompLabels = {};    % mode 3: canonical labels built from parameters

if isSpecificLevelMode
    allSpecified = all(cellfun(@(x) ~isempty(x), params.CompareLevels)) && ...
        all(cellfun(@(x) ~isempty(x), catList));
    if allSpecified
        [fixedCompLabels, ~] = buildMode3Items(stimsNeeded, catList, params.CompareLevels);
    end
end

% =========================================================================
% SECTION 4 — PER-EXPERIMENT LOOP
% =========================================================================

if runLoop

    % --- BUG FIX (was Bug #4): Map-based counters ---
    % Using containers.Map guarantees the same animal/insertion always
    % receives the same numeric index, regardless of expList ordering.
    % The previous sequential-counter approach silently assigned different
    % IDs if the same animal appeared non-contiguously in expList.
    animalMap    = containers.Map('KeyType','char','ValueType','double');
    insertionMap = containers.Map('KeyType','char','ValueType','double');
    nextAnimalIdx    = 0;   % running counter for unique animals
    nextInsertionIdx = 0;   % running counter for unique (animal, insertion) pairs

    for ex = expList

        % ---- 4a: Load recording and check stimulus availability ----

        NP = loadNPclassFromTable(ex);                       % load Neuropixels recording object
        fprintf('Processing recording: %s\n', NP.recordingName);

        % Load one analysis object per unique stimulus class
        [vsObjs, present] = loadStimulusObjects(NP, stimsNeeded);

        % Check that every required stimulus is present in this recording
        allPresent = true;
        for si = 1:numel(stimsNeeded)
            if ~present(stimsNeeded{si})
                allPresent = false;  break              % at least one missing — skip
            end
        end
        if ~allPresent
            fprintf('  -> Skipping: stimulus not present.\n');
            continue
        end

        % ---- 4b: Mode-specific session selection ----

        if isCategoryMode
            % Mode 2: find session of stimName with >=2 levels of catName
            key   = getObjKey(stimName);                % shared object key (e.g. 'SDG')
            vsObj = vsObjs(key);                        % retrieve the analysis object

            % Search sessions for >=2 category levels
            [levels, ~, vsObj] = findCategoryLevels(vsObj, NP, stimName, catName, params);

            if numel(levels) < 2
                fprintf('  -> Skipping: <2 levels for "%s".\n', catName);
                continue                                % not enough levels for comparison
            end

            vsObjs(key) = vsObj;                        % store back (may have changed session)

            % Convert numeric levels to field-name labels (e.g. 'size_5')
            currentLabels = arrayfun(@(v) levelToFieldName(catName, v), ...
                levels, 'UniformOutput', false);

            % Lock the level set on the first valid recording; skip any
            % subsequent recordings with a different level set.
            if isempty(levelLabels)
                levelLabels = currentLabels;
                fprintf('  Category levels locked: %s\n', strjoin(levelLabels, ', '));
            else
                if ~isequal(sort(currentLabels), sort(levelLabels))
                    fprintf('  -> Skipping: level mismatch (expected %s, got %s).\n', ...
                        strjoin(levelLabels,','), strjoin(currentLabels,','));
                    continue
                end
            end

       elseif isSpecificLevelMode
            sessionFound = true;
            for si = 1:numel(stimsNeeded)
                sn   = stimsNeeded{si};
                cat  = catList{si};
                lvls = params.CompareLevels{si};
                key  = getObjKey(sn);

                if isempty(cat)
                    % No category — general stats will be used in 4d
                    continue

                elseif isempty(lvls)
                    % Category specified but no levels → auto-discover all
                    vsObj = vsObjs(key);
                    [discoveredLvls, ~, vsObj] = findCategoryLevels( ...
                        vsObj, NP, sn, cat, params);
                    if numel(discoveredLvls) < 2
                        fprintf('  -> Skipping: <2 levels for "%s" in %s.\n', cat, sn);
                        sessionFound = false; break
                    end
                    params.CompareLevels{si} = discoveredLvls(:)';
                    vsObjs(key) = vsObj;

                else
                    % Specific levels → find session containing all of them
                    [vsObj, allFound] = findSessionWithLevels( ...
                        NP, sn, cat, lvls, params);
                    if ~allFound
                        fprintf('  -> Skipping: %s levels [%s] of "%s" not found.\n', ...
                            sn, num2str(lvls(:)','%g '), cat);
                        sessionFound = false; break
                    end
                    vsObjs(key) = vsObj;
                end
            end
            if ~sessionFound, continue; end

            % Build labels from (possibly auto-discovered) levels
            [currentLabels, ~] = buildMode3Items(stimsNeeded, catList, params.CompareLevels);
            if isempty(fixedCompLabels)
                fixedCompLabels = currentLabels;
                fprintf('  Labels locked: %s\n', strjoin(fixedCompLabels, ', '));
            end
        end

        % ---- 4c: Parse metadata and assign animal/insertion indices ----
        %
        % BUG FIX (was Bug #4): Uses map-based lookup instead of sequential
        % counters.  Safe for any ordering of expList.

        % Extract animal ID from recording name (e.g. 'PV123' or 'SA45')
        animalID = string(regexp(NP.recordingName, 'PV\d+', 'match', 'once'));
        if animalID == ""
            animalID = string(regexp(NP.recordingName, 'SA\d+', 'match', 'once'));
        end

        % Extract insertion number from the directory path
        insStr = regexp(NP.recordingDir, 'Insertion\d+', 'match', 'once');
        insNum = str2double(regexp(insStr, '\d+', 'match'));

        % Real (recording) insertion number = trailing token of the recording name,
        % e.g. 'PV103_Experiment_12_6_23_1' -> 1. This is the human-facing insertion
        % number used downstream to mark example units; robust to dropped insertions
        % because it is the actual recording label, not a sequential/offset index.
        recName  = NP.recordingName;                        % == vsObj.dataObj.recordingName
        nameTok  = strsplit(recName, '_');                  % {'PV103','Experiment','12','6','23','1'}
        realInsertionID = str2double(nameTok{end});         % trailing token is the insertion number
        assert(~isnan(realInsertionID), ...
            'recordingName does not end in a numeric insertion token: %s', recName);

        % Register animal in the map (assigns a new index only on first encounter)
        animalKey = char(animalID);
        if ~animalMap.isKey(animalKey)
            nextAnimalIdx = nextAnimalIdx + 1;
            animalMap(animalKey) = nextAnimalIdx;
        end

        % Register the (animal, insertion) pair in the map
        insKey = sprintf('%s__Ins%d', animalKey, insNum);
        if ~insertionMap.isKey(insKey)
            nextInsertionIdx = nextInsertionIdx + 1;
            insertionMap(insKey) = nextInsertionIdx;
        end
        insertionCount = insertionMap(insKey);          % stable numeric index for this insertion

        % ---- 4d: Run statistics and extract per-item data ----

        stimData       = struct();      % per-item z-scores, p-values, spike rates
        nUnits         = [];            % total somatic neuron count (set once)
        compLabels     = {};            % labels for items being compared
        generalPbyStim = struct();      % general per-neuron p-values (for useGeneralFilter)

        if isCategoryMode
            % ---- Mode 2: single stimulus, all levels of one category ----

            key   = getObjKey(stimName);
            vsObj = vsObjs(key);

            % BUG FIX: runStimStats requires 4 args (catLevels/lockSpec have no
            % default and are referenced unconditionally inside) — the previous
            % 2-arg call crashed with "Not enough input arguments" every time
            % Mode 2 ran. Fetch this stimulus's focus spec, same as Mode 1/3.
            sp = specByKey(key);                        % this stimulus's CategoryLevels/Lock focus
            % BUG FIX: capture and reuse this result directly. extractStimData
            % used to re-derive stats via a bare vsObj.StatisticsPerNeuron()
            % call (all defaults), which silently discarded CategoryLevels/Lock
            % (and any other non-default param) via a mismatched recompute.
            stats = runStimStats(vsObj, params, sp.catLevels, sp.lock);  % general per-neuron statistics, with focus
            vsObjs(key) = vsObj;                        % store back (handle class — redundant but safe)

            % Extract general p-values (for optional useGeneralFilter)
            [~, generalP, ~, ~] = extractStimData(stats, vsObj, stimName, params.useZmean, params.useTtest);
            generalPbyStim.(stimName) = generalP;

            % Auto-detect GratingType from stimulus abbreviation
            gratingType = detectGratingType(stimName);

            % Build name-value arguments for StatisticsPerNeuronPerCategory
            % Lock (sp.lock) threads this stimulus's focus spec into the
            % per-level trial sets, so a Lock used to select responsive
            % neurons is also honoured when building the compared levels.
            catStatsArgs = { ...
                'compareCategory', catName, ...
                'nBoot',           params.nBootCategory, ...
                'overwrite',       params.overwriteStats, ...
                'BaseRespWindow',  params.BaseRespWindow, ...
                'applyFDR',        params.applyFDR,...
                'Lock',            sp.lock, ...
                'CategoryMaximized'      params.CategoryMaximized};
            if ~isempty(gratingType)
                % Only pass GratingType for grating stimuli (SDGm/SDGs)
                catStatsArgs = [catStatsArgs, {'GratingType', gratingType}];
            end

            % Run per-category statistics
            catStats = vsObj.StatisticsPerNeuronPerCategory(catStatsArgs{:});

            % Extract per-level data from catStats
            for li = 1:numel(levelLabels)
                fName = levelLabels{li};                        % field name in catStats
                stimData.(fName).z    = catStats.(fName).ZScoreU(:);       % z-score
                if params.useTtest
                     stimData.(fName).p    = catStats.(fName).pValTTest(:); % p-value
                else
                    stimData.(fName).p    = catStats.(fName).pvalsResponse(:); % p-value
                end
                stimData.(fName).spkR = catStats.(fName).ObsStat(:);       % spike rate / z_mean
                if isempty(nUnits), nUnits = numel(stimData.(fName).z); end
            end
            compLabels = levelLabels;                   % items to compare = level labels

        elseif isSpecificLevelMode
            % ---- Mode 3: each stimulus contributes one or more (stim, cat, level) items ----

            for si = 1:numel(stimsNeeded)
                sn   = stimsNeeded{si};                 % stimulus name
                cat  = catList{si};                     % category for this stimulus
                lvls = params.CompareLevels{si};        % requested levels
                key  = getObjKey(sn);
                vsObj = vsObjs(key);

                % Run general per-neuron statistics, with focus
                % BUG FIX: capture and reuse this result directly, rather than
                % letting extractStimData re-derive stats via a bare
                % vsObj.StatisticsPerNeuron() call (all defaults), which
                % silently discarded CategoryLevels/Lock via a mismatched recompute.
                stats = runStimStats(vsObj, params, catLevelsPerStim{si}, lockPerStim{si});  % single call, with focus
                vsObjs(key) = vsObj;

                if isempty(cat)
                    % The call above already produced the correct,
                    % focus-respecting results for this stimulus — reuse them.
                    % BUG FIX: useTtest was previously omitted here, but
                    % extractStimData references useTTest unconditionally
                    % (line ~1547/1588) — this call would have thrown
                    % "Not enough input arguments" the moment this branch ran.
                    [z, p, spkR, ~] = extractStimData(stats, vsObj, sn, params.useZmean, params.useTtest);
                    generalPbyStim.(sn) = p;
                    stimData.(sn).z    = z(:);
                    stimData.(sn).p    = p(:);
                    stimData.(sn).spkR = spkR(:);
                    compLabels{end+1}  = sn;                    %#ok<AGROW>
                    if isempty(nUnits), nUnits = numel(z); end
                    continue
                end

                vsObjs(key) = vsObj;

                % Extract general p-values (for optional useGeneralFilter)
                [~, generalP, ~, ~] = extractStimData(stats, vsObj, sn, params.useZmean, params.useTtest);
                generalPbyStim.(sn) = generalP;

                % Auto-detect GratingType from stimulus abbreviation
                gratingType = detectGratingType(sn);

                % Build name-value arguments for StatisticsPerNeuronPerCategory
                % Lock (lockPerStim{si}) threads this stimulus's focus spec
                % into the per-level trial sets, matching the general filter.
                catStatsArgs = { ...
                    'compareCategory', cat, ...
                    'nBoot',           params.nBootCategory, ...
                    'overwrite',       params.overwriteStats, ...
                    'BaseRespWindow',  params.BaseRespWindow, ...
                    'applyFDR',        params.applyFDR, ...
                    'Lock',            lockPerStim{si}};
                if ~isempty(gratingType)
                    catStatsArgs = [catStatsArgs, {'GratingType', gratingType}];
                end

                % Run per-category statistics for this stimulus
                catStats = vsObj.StatisticsPerNeuronPerCategory(catStatsArgs{:});

                % Fields corresponding to levels
                allFields = fieldnames(catStats);

                % Keep only fields starting with category name
                levelFields = allFields(startsWith(lower(allFields), lower(cat) + "_"));

                % Numeric levels stored in struct
                storedLvls = catStats.categoryLevels(:);

                % Extract data for each requested level
                for lvi = 1:numel(lvls)
                    lv      = lvls(lvi);                        % numeric level value
                    % Find closest matching stored level
                    [~, idx] = min(abs(storedLvls - lv));

                    % Corresponding field name
                    fName = levelFields{idx};
                    cLabel  = makeCompLabel(sn, cat, lv);       % short composite label

                    stimData.(cLabel).z    = catStats.(fName).ZScoreU(:);
                    stimData.(cLabel).p    = catStats.(fName).pvalsResponse(:);
                    stimData.(cLabel).spkR = catStats.(fName).ObsStat(:);

                    compLabels{end+1} = cLabel;                 %#ok<AGROW>
                    if isempty(nUnits), nUnits = numel(stimData.(cLabel).z); end
                end
            end

            % Verify labels match (skip check if labels were just locked above)
            if ~isempty(fixedCompLabels) && ~isequal(sort(compLabels(:)), sort(fixedCompLabels(:)))
                fprintf('  -> Skipping: comparison label mismatch.\n');
                continue
            end

        else
            % ---- Mode 1: across-stimulus (one item per stimulus) ----

            % Run general statistics for every loaded stimulus object
            % BUG FIX: capture and store each object's results directly (keyed
            % the same way as vsObjs) instead of letting extractStimData
            % re-derive stats via a bare vsObj.StatisticsPerNeuron() call (all
            % defaults), which silently discarded CategoryLevels/Lock via a
            % mismatched recompute.
            objKeys = keys(vsObjs);
            statsByKey = containers.Map();
            for k = 1:numel(objKeys)
                key   = objKeys{k};
                vsObj = vsObjs(key);
                sp    = specByKey(key);
                statsByKey(key) = runStimStats(vsObj, params, sp.catLevels, sp.lock);
                vsObjs(key) = vsObj;
            end

            % Extract z-scores, p-values, spike rates for each stimulus
            for si = 1:numel(stimsNeeded)
                sn  = stimsNeeded{si};
                key = getObjKey(sn);
                [z, p, spkR, ~] = extractStimData(statsByKey(key), vsObjs(key), sn, params.useZmean, params.useTtest);
                stimData.(sn).z    = z(:);
                stimData.(sn).p    = p(:);
                stimData.(sn).spkR = spkR(:);
                generalPbyStim.(sn) = p(:);             % general p = per-stimulus p in mode 1
                if isempty(nUnits), nUnits = numel(z); end
            end
            compLabels = stimsNeeded;                   % items to compare = stimulus names
        end

        % ---- 4e: Optional post-hoc FDR correction ----
        % NOTE: This is the AllExpAnalysis-level FDR (Benjamini-Hochberg on
        % the extracted p-values).  It is DISTINCT from params.applyFDR,
        % which is applied inside the statistics functions themselves.
        if params.useFDR
            for ci = 1:numel(compLabels)
                cl = compLabels{ci};
                stimData.(cl).p = bhFDR(stimData.(cl).p);  % correct per-item p-values
            end
        end

        % ---- 4f: Significance mask ----
        % Build a logical OR mask: a neuron passes if it is significant for
        % at least one comparison item.

        if params.useGeneralFilter && (isCategoryMode || isSpecificLevelMode)
            % Use general per-stimulus p-values (StatisticsPerNeuron), OR'd
            % across stimuli.  This avoids filtering on the same per-level
            % p-values used for comparison.
            orMask = false(nUnits, 1);
            stimNames_ = fieldnames(generalPbyStim);
            for si = 1:numel(stimNames_)
                gp = generalPbyStim.(stimNames_{si});
                if params.useFDR, gp = bhFDR(gp); end  % apply FDR if requested
                orMask = orMask | (gp < params.threshold);
            end
        else
            % Default: OR across all per-item p-values
            orMask = false(nUnits, 1);
            for ci = 1:numel(compLabels)
                cl = compLabels{ci};
                orMask = orMask | (stimData.(cl).p < params.threshold);
            end
        end

        unitIDs = find(orMask);     % indices of neurons passing the filter
        nSig    = numel(unitIDs);   % count of significant neurons

        % ---- 4g: Append to TableStimComp ----
        % Add one block of rows per comparison item (only significant neurons).
        if nSig > 0
            for ci = 1:numel(compLabels)
                cl = compLabels{ci};
                newRows = table( ...
                    repmat(categorical(cellstr(animalID)), nSig, 1), ...   % animal column
                    repmat(categorical(insertionCount),    nSig, 1), ...   % insertion column
                    repmat(categorical(realInsertionID),   nSig, 1), ...   %Real Insertion ID
                    repmat(categorical(cellstr(cl)),       nSig, 1), ...   % stimulus/item column
                    categorical(unitIDs), ...                                          % neuron ID column
                    stimData.(cl).z(orMask), ...                            % z-score column
                    stimData.(cl).spkR(orMask), ...                         % spike rate column
                    'VariableNames', {'animal','insertion','realInsertion','stimulus','NeurID','Z-score','SpkR'});
                TableStimComp = [TableStimComp; newRows];                   %#ok<AGROW>
            end
        end

        % ---- 4h: Append to TableRespNeurs ----
        % One summary row per (insertion x comparison item): responsive
        % neuron count and total somatic neuron count.
        for ci = 1:numel(compLabels)
            cl    = compLabels{ci};
            nResp = sum(stimData.(cl).p < params.threshold);  % neurons below threshold
            newRow = table( ...
                categorical(cellstr(animalID)), ...             % animal
                categorical(insertionCount), ...                % insertion
                categorical(cellstr(cl)), ...                   % stimulus/item
                nResp, nUnits, ...                              % responsive count, total count
                'VariableNames', TableRespNeurs.Properties.VariableNames);
            TableRespNeurs = [TableRespNeurs; newRow];          %#ok<AGROW>
        end

        fprintf('  -> %d / %d units pass filter.\n', nSig, nUnits);

    end  % end for ex

    % =========================================================================
    % SECTION 5 — SAVE POOLED DATA
    % =========================================================================

    % Describe the analysis mode in plain text (for figure annotation)
    switch mode
        case 1, modeDesc = 'across-stimulus';
        case 2, modeDesc = 'within-stimulus-category';
        case 3, modeDesc = 'specific-level-cross-stim';
    end

    % Build a metadata sub-struct capturing every parameter needed to
    % reproduce the analysis.  This travels with the pooled data AND with
    % saved figures, so any figure can be traced to its configuration.
    analysisMetadata = struct( ...
        'mode',             mode, ...                   % numeric mode (1/2/3)
        'modeDescription',  modeDesc, ...               % human-readable label
        'ComparePairs',     {params.ComparePairs}, ...  % stimuli being compared
        'RespDurationWin',  params.RespDurationWin, ... % ResponseWindow duration (ms)
        'BaseRespWindow',   params.BaseRespWindow, ...  % statistics base window (ms)
        'SpatialGridMode',  params.SpatialGridMode, ... % whether grid mode is active
        'maxCategory',      params.maxCategory, ...     % max-category flag for StatisticsPerNeuron
        'applyFDR',         params.applyFDR, ...        % FDR inside statistics functions
        'useFDR',           params.useFDR, ...          % post-hoc BH FDR in AllExpAnalysis
        'threshold',        params.threshold, ...       % responsiveness p-value cutoff
        'useZmean',         params.useZmean, ...        % z_mean vs peak spike rate
        'useGeneralFilter', params.useGeneralFilter, ...% general vs per-item filter
        'nBoot',            params.nBoot, ...           % bootstrap iterations (pairwise)
        'nBootCategory',    params.nBootCategory,...    % bootstrap iterations (category)
        'durationWindow',   params.durationWindow,...     % Window lenght for moving window
        'CategoryLevels',   {params.CategoryLevels}, ...  % focus: separate factors per stimulus
        'Lock',             {params.Lock});            % focus: locked factor values per stimulus

    % Add mode-specific fields
    if isCategoryMode
        analysisMetadata.stimName    = stimName;        % the decomposed stimulus
        analysisMetadata.catName     = catName;         % category column name
        analysisMetadata.levelLabels = levelLabels;     % resolved level labels
    end
    if isSpecificLevelMode
        analysisMetadata.catList         = catList;
        analysisMetadata.CompareLevels   = params.CompareLevels;
        analysisMetadata.fixedCompLabels = fixedCompLabels;
    end

    % Pack into save struct
    S.expList           = expList;
    S.TableStimComp     = TableStimComp;
    S.TableRespNeurs    = TableRespNeurs;
    S.params            = params;               % full params (superset of metadata)
    S.mode              = mode;
    S.analysisMetadata  = analysisMetadata;      % curated subset for figure annotation
    if isCategoryMode,      S.levelLabels     = levelLabels;     end
    if isSpecificLevelMode, S.fixedCompLabels = fixedCompLabels; end

    % Write to disk
    save(savePath, '-struct', 'S');
    fprintf('Saved pooled data to %s\n', savePath);
end

% =========================================================================
% SECTION 5b — RESTORE VARIABLES FROM CACHE
% =========================================================================
% BUG FIX (was Bug #6): When runLoop is false, S was loaded from disk but
% mode-dependent local variables were never set, causing downstream errors.

if ~runLoop
    mode = S.mode;                                      % restore numeric mode
    isCategoryMode      = (mode == 2);
    isSpecificLevelMode = (mode == 3);

    if isCategoryMode && isfield(S, 'levelLabels')
        levelLabels = S.levelLabels;                    % restore level labels
    end
    if isSpecificLevelMode && isfield(S, 'fixedCompLabels')
        fixedCompLabels = S.fixedCompLabels;            % restore comparison labels
    end

    fprintf('Loaded pooled data from cache: %s\n', savePath);
end

% =========================================================================
% SECTION 6 — GUARD: ABORT IF NO DATA
% =========================================================================

if isempty(S.TableStimComp) || height(S.TableStimComp) == 0
    warning('AllExpAnalysis:noUnits', 'No significant units found. Returning empty.');
    tempTable    = table();
    summaryStats = struct('Zscore', [], 'SpkR', [], 'FracResponsive', []);  % keep the output defined on every path
    return
end

% Replace NaN z-scores and spike rates with zero (prevents plotting issues)
S.TableStimComp.('Z-score')(isnan(S.TableStimComp.('Z-score'))) = 0;
S.TableStimComp.SpkR(isnan(S.TableStimComp.SpkR))               = 0;

% Defensive: ensure animal/insertion/stimulus are string-based categoricals.
% Prevents numeric-named categoricals from being silently converted to
% double inside plotSwarmBootstrapWithComparisons.
S.TableStimComp.animal    = categorical(cellstr(string(S.TableStimComp.animal)));
S.TableStimComp.insertion = categorical(cellstr(string(S.TableStimComp.insertion)));
S.TableStimComp.stimulus  = categorical(cellstr(string(S.TableStimComp.stimulus)));

% =========================================================================
% SECTION 7 — SHARED PLOTTING SETUP
% =========================================================================

% Create a fresh analysis object for the first experiment (used for printFig)
NP = loadNPclassFromTable(expList(end));
vs = linearlyMovingBallAnalysis(NP, 'MultipleOffsets', false, 'Multiplesizes', false);

[~, figTag, ~] = fileparts(nameOfFile);   % e.g. 'Ex_49-97_SpecLvl_SDGm-angles-0_1p57...'
% Build a consistent colour map: one colour per animal, shared across all plots
animalOrder  = categories(S.TableStimComp.animal);      % sorted unique animal names
nAnimals     = numel(animalOrder);                       % number of animals
sharedCmap   = lines(nAnimals);                          % Nx3 colour matrix
animalIdxAll = double(S.TableStimComp.animal);           % per-row animal index (for scatter colouring)

if (mode == 3 || mode == 2) && ...
        (any(contains(string(S.TableStimComp.stimulus), "SDG")) && any(contains(string(S.TableStimComp.stimulus), "ang")))

    stimStr = string(S.TableStimComp.stimulus);

    % Rows containing "SDG"
    idxSDG = contains(stimStr, "SDG");

    % Replace angle values only in SDG rows
    stimStr(idxSDG) = replace(stimStr(idxSDG), ...
        ["0",  "90", "180", "270"], ...
        ["1.57","0", "4.71","3.14"]);

    % Convert back to categorical
    S.TableStimComp.stimulus = categorical(stimStr);

end


% Smallest p-value the bootstrap can resolve: the two-tailed statistic doubles
% the smaller tail, so with zero counts on one side the minimum non-zero value
% is 2/nBootPval. Passed to the plotting and printing helpers so a p of exactly
% 0 is reported as '< that resolution' rather than against a hard-coded floor.
pValRes = 2 / params.nBootPval;

% All pairwise combinations of comparison items
compLabels = cellstr(categories(S.TableStimComp.stimulus));  % unique sorted item labels
pairsAll   = nchoosek(compLabels, 2);                        % Kx2 cell of pairs

% Display-name substitutions for axis labels
labelMap = {'RG','SB'; 'SDGs','SG'; 'SDGm','MG'};

% =========================================================================
% SECTION 8 — Z-SCORE PAIRWISE COMPARISON
% =========================================================================

% Compute a hierarchical-bootstrap p-value for each pair (z-score metric).
% Seed first so the reported p-values are reproducible and independent of the
% global RNG state left by any earlier plotting call.
rng(params.pvalSeed);
pValsZ = zeros(1, size(pairsAll, 1));
for pi = 1:size(pairsAll, 1)
    pValsZ(pi) = bootstrapPairDifference( ...
        S.TableStimComp, pairsAll(pi,:), params.nBootPval, 'Z-score');
end

% Upper y-limit: ceiling of max z-score plus headroom for significance brackets
ZscoreYlim = ceil(max(S.TableStimComp.('Z-score'))) +0.1*ceil(max(S.TableStimComp.('Z-score')));

% Generate swarm + bootstrap plot for z-scores
[fig,~,figAllDiffs,summaryStats.Zscore] = plotSwarmBootstrapWithComparisons( ...
    S.TableStimComp, pairsAll, pValsZ, {'Z-score'}, ...
    yLegend = 'Z-score', yMaxVis = ZscoreYlim, pValRes = pValRes, ...
    diff = true, plotMeanSem = true, Alpha = 0.7,   markUnits = params.markUnits)  ;        % N×4 {NeurID,stimulus,animal,insertion}; row1->X, row2->+);

% Echo the numbers behind the panel: mean/median, the plotted 95% CI, and the
% effect sizes with their bootstrap CIs, per stimulus and per pair.
printStatsSummary('Z-score', summaryStats.Zscore, pValRes);

% Apply consistent font formatting to all axes in the figure
ax = findobj(fig, 'Type', 'axes');
for i = 1:numel(ax)
    formatAxes(ax(i), 8, 'helvetica');
end

% Set figure size and colour map
figure(fig);
set(fig, 'Units', 'centimeters', 'Position', [20 20 10 6]);
colormap(fig, sharedCmap);

%set(gca,'Clipping','off')
pause(2)
% Save figure if PaperFig mode is active
if params.PaperFig
    vs.printFig(fig, sprintf('Zscore-Swarm-%s', figTag), ...
        PaperFig = params.PaperFig);
end

ax = findobj(figAllDiffs, 'Type', 'axes');
for i = 1:numel(ax)
    formatAxes(ax(i), 8, 'helvetica');
end

if ~isempty(figAllDiffs)
    figure(figAllDiffs);
    pause(2)
    % Save figure if PaperFig mode is active
    if params.PaperFig
        vs.printFig(figAllDiffs, sprintf('AllDiff-Zscore-Swarm-%s', figTag), ...
            PaperFig = params.PaperFig);
    end
end

groupcounts(S.TableStimComp, 'stimulus')                       % pooled rows per stimulus: are MB and MBR equal?
g = groupcounts(S.TableStimComp, {'insertion','stimulus'});    % rows per (insertion, stimulus)
unstack(g, 'GroupCount', 'stimulus')                           % one row per insertion; spot where MB ~= MBR
% For exactly two items, also produce a paired scatter plot
if numel(compLabels) == 2
    fig = plotPairScatter(S.TableStimComp, compLabels, ...
        'Z-score', sharedCmap, animalIdxAll, labelMap);
    title('Z-score');
    if params.PaperFig
        vs.printFig(fig, sprintf('Zscore-Scatter-%s', figTag), ...
            PaperFig = params.PaperFig);
    end
end

% =========================================================================
% SECTION 9 — SPIKE-RATE PAIRWISE COMPARISON
% =========================================================================

% Compute a hierarchical-bootstrap p-value for each pair (spike rate metric).
% Re-seed: this loop runs after the z-score panels were drawn, so without it the
% p-values would depend on how many random draws that plotting consumed.
rng(params.pvalSeed);
pValsSpk = zeros(1, size(pairsAll, 1));
for pi = 1:size(pairsAll, 1)
    pValsSpk(pi) = bootstrapPairDifference( ...
        S.TableStimComp, pairsAll(pi,:), params.nBootPval, 'SpkR');
end

% Upper y-limit for spike rate
spkMax = max(S.TableStimComp.SpkR) +0.1*max(S.TableStimComp.SpkR);

% Generate swarm + bootstrap plot for spike rates
[fig,~,figAllDiffs,summaryStats.SpkR] = plotSwarmBootstrapWithComparisons( ...
    S.TableStimComp, pairsAll, pValsSpk, {'SpkR'}, ...
    yLegend = 'SpkR', yMaxVis = spkMax, pValRes = pValRes, ...
    diff = true, plotMeanSem = true, Alpha = 0.7, markUnits = params.markUnits);

% Echo the spike-rate summary (same layout as the z-score one)
printStatsSummary('Spike rate', summaryStats.SpkR, pValRes);

% Apply consistent font formatting
ax = findobj(fig, 'Type', 'axes');
for i = 1:numel(ax)
    formatAxes(ax(i), 8, 'helvetica');
end

% Set figure size and colour map
figure(fig);
set(fig, 'Units', 'centimeters', 'Position', [20 20 10 6]);
colormap(fig, sharedCmap);
pause(2)
% Save figure if PaperFig mode is active
if params.PaperFig
    vs.printFig(fig, sprintf('SpkRate-Swarm-%s', figTag), ...
        PaperFig = params.PaperFig);
end

ax = findobj(figAllDiffs, 'Type', 'axes');
for i = 1:numel(ax)
    formatAxes(ax(i), 8, 'helvetica');
end

if ~isempty(figAllDiffs)
    figure(figAllDiffs);
    pause(2)

    % Save figure if PaperFig mode is active
    if params.PaperFig
        vs.printFig(figAllDiffs, sprintf('AllDiff-SpkRate-Swarm-%s', figTag), ...
            PaperFig = params.PaperFig);
    end
end

% For exactly two items, produce a paired scatter plot
if numel(compLabels) == 2
    fig = plotPairScatter(S.TableStimComp, compLabels, ...
        'SpkR', sharedCmap, animalIdxAll, labelMap);
    title('Spk. rate');
    if params.PaperFig
        vs.printFig(fig, sprintf('SpkRate-Scatter-%s', figTag), ...
            PaperFig = params.PaperFig);
    end
end

% =========================================================================
% SECTION 10 — FRACTION-RESPONSIVE ANALYSIS
% =========================================================================

compLabels = cellstr(categories(S.TableRespNeurs.stimulus));  % unique sorted item labels

% Find groups by insertion, then check which insertions contain ALL items
[G, ~] = findgroups(S.TableRespNeurs.insertion);
hasAll  = splitapply( ...
    @(s) all(ismember(categorical(compLabels), s)), ...
    S.TableRespNeurs.stimulus, G);

% Filter to only insertions that have data for every comparison item
tempTable = S.TableRespNeurs( ...
    hasAll(G) & ismember(S.TableRespNeurs.stimulus, categorical(compLabels)), :);

% --- BUG FIX (was Bug #3): Hierarchical bootstrap for fraction-responsive ---
% The previous flat bootstrp(@mean, diffs) ignored the nesting of insertions
% within animals.  Using hierBoot is consistent with the mixed model.

pairsAll   = nchoosek(compLabels, 2);

% Re-seed: this loop runs last, after both sets of swarm panels were drawn.
% Without it the fraction p-values inherit the RNG state those plots left, so
% a p sitting near a threshold could change on re-run with identical data.
rng(params.pvalSeed);
pValsFrac = zeros(1, size(pairsAll, 1));
for pi = 1:size(pairsAll, 1)
    % Per-insertion paired fraction difference + hierarchical bootstrap.
    % Extracted to its own function so the primary result, the sensitivity
    % analysis below, and any downstream summary table all use one
    % implementation and cannot drift apart.
    pValsFrac(pi) = bootstrapFractionDifference( ...
        S.TableRespNeurs, pairsAll(pi,:), params.nBootPval);
end

% Compute total responsive neurons per insertion (across stimuli)
[G, ~]                  = findgroups(tempTable.insertion);
totals                  = splitapply(@sum, tempTable.respNeur, G);
tempTable.TotalRespNeur = totals(G);

% Generate swarm plot for fraction responsive
[fig,~,figAllDiffs,summaryStats.FracResponsive] = plotSwarmBootstrapWithComparisons( ...
    tempTable, pairsAll, pValsFrac, ...
    {'respNeur','totalSomaticN'}, ...
    fraction = true, showBothAndDiff = false, pValRes = pValRes, ...
    yLegend = 'Responsive/total units', ...
    diff = false, filled = false, Xjitter = 'none', ...
    Alpha = 0.6, drawLines = true);

% Echo the fraction-responsive summary (insertion-level, so n = # insertions)
printStatsSummary('Fraction responsive', summaryStats.FracResponsive, pValRes);

% =========================================================================
% SECTION 10b — FRACTION-RESPONSIVE SENSITIVITY ANALYSIS
% =========================================================================
% Insertions with very few somatic units give a quantised fraction estimate
% (2 units -> only 0, 0.5 or 1). They cannot be down-weighted by the
% bootstrap, because this panel supplies exactly one value per insertion, so
% every insertion counts equally whether it held 2 units or 113. That
% quantisation noise inflates the SD and therefore ATTENUATES Cohen's d.
%
% The primary result above is left untouched (no exclusions to justify in
% Methods). This recomputes the same comparison over the higher-yield
% insertions only, so both numbers come out of a single run and the paper can
% state that the conclusion does not depend on the low-yield sessions.
%
% NOTE: this is a guard on the DENOMINATOR (units recorded), which is set by
% insertion quality and is independent of the visual response. Never filter on
% the NUMERATOR (units responsive) — that selects on the dependent variable and
% removes precisely the insertions where MB and SB agreed.
if params.minSomaticN > 0

    % Insertions passing the yield threshold (evaluated per row; totalSomaticN
    % is constant across the stimulus rows of one insertion)
    keepRows = tempTable.totalSomaticN >= params.minSomaticN;
    tempTableSens = tempTable(keepRows, :);

    nInsAll  = numel(unique(tempTable.insertion));           % insertions in the primary analysis
    nInsKept = numel(unique(tempTableSens.insertion));       % insertions surviving the threshold

    if nInsKept < 3
        % Too little left to bootstrap; say so rather than emitting a stat
        warning('AllExpAnalysis:sensitivityTooSmall', ...
            ['minSomaticN=%g leaves only %d insertion(s); skipping the ', ...
             'fraction-responsive sensitivity analysis.'], ...
            params.minSomaticN, nInsKept);
    else
        % ---- Filtered p-value: same loop and same seed as the primary ----
        rng(params.pvalSeed);
        pValsFracSens = zeros(1, size(pairsAll, 1));
        for pi = 1:size(pairsAll, 1)
            % Same routine as the primary p-value above, on the filtered table
            pValsFracSens(pi) = bootstrapFractionDifference( ...
                tempTableSens, pairsAll(pi,:), params.nBootPval);
        end

        % ---- Filtered CI and effect sizes ----
        % Reuse the plotting routine rather than duplicating the statistics:
        % computeSwarmStats is a local function of plotSwarmBootstrapWithComparisons
        % and cannot be called directly, and going through the same entry point
        % guarantees these numbers are produced identically to the primary ones.
        % The figure it creates is closed immediately — only the stats are kept.
        [figSens, ~, figSensDiffs, sensStats] = plotSwarmBootstrapWithComparisons( ...
            tempTableSens, pairsAll, pValsFracSens, ...
            {'respNeur','totalSomaticN'}, ...
            fraction = true, showBothAndDiff = false, pValRes = pValRes, ...
            yLegend = 'Responsive/total units', ...
            diff = false, filled = false, Xjitter = 'none', ...
            Alpha = 0.6, drawLines = true);

        close(figSens);                                    % discard the transient figure
        if ~isempty(figSensDiffs) && isgraphics(figSensDiffs)
            close(figSensDiffs);                           % and the all-pairs one, if any
        end

        % Record alongside the primary result
        sensStats.minSomaticN = params.minSomaticN;        % threshold applied
        sensStats.nDropped    = nInsAll - nInsKept;        % insertions removed
        summaryStats.FracResponsive.sensitivity = sensStats;

        fprintf(['\n### SENSITIVITY ANALYSIS: insertions with >= %g somatic units ', ...
                 '(%d of %d kept, %d dropped).\n', ...
                 '### The primary result above uses ALL insertions; this is a ', ...
                 'robustness check, not a replacement.\n'], ...
            params.minSomaticN, nInsKept, nInsAll, nInsAll - nInsKept);
        printStatsSummary(sprintf('Fraction responsive [minSomaticN>=%g]', ...
            params.minSomaticN), sensStats, pValRes);
    end
end

set(findobj(fig, 'Type', 'axes'), 'Clipping', 'off');   % target the primary panel explicitly: gca is unreliable here now that the sensitivity block above creates and closes figures

% Count total unique responsive neurons across all (animal, insertion) pairs
totalResp = 0;
animals   = unique(S.TableStimComp.animal);
for a = 1:numel(animals)
    inser = unique(S.TableStimComp.insertion(S.TableStimComp.animal == animals(a)));
    for in = 1:numel(inser)
        totalResp = totalResp + length(unique( ...
            S.TableStimComp.NeurID( ...
                S.TableStimComp.insertion == inser(in) & ...
                S.TableStimComp.animal    == animals(a))));
    end
end


% --- Compute total somatic neurons across all unique insertions ---
totalSomatic   = 0;                                                     % initialize accumulator
allInsertions  = unique(S.TableRespNeurs.insertion);                     % unique insertion IDs across all animals
for in = 1:numel(allInsertions)                                         % iterate over every insertion (BUG FIX: was "= numel(inser)", only hit last index)
    insRows = S.TableRespNeurs.totalSomaticN( ...                       % extract totalSomaticN for all rows matching this insertion
        S.TableRespNeurs.insertion == allInsertions(in));
    totalSomatic = totalSomatic + insRows(1);                           % take first value (identical across stimulus rows within an insertion)
end

% Build annotation string with per-item responsive neuron counts
perItemN  = arrayfun(@(s) sum(tempTable.respNeur(tempTable.stimulus == s)), ...
    categorical(compLabels));
annotParts = arrayfun(@(i) sprintf('%s = %d', compLabels{i}, perItemN(i)), ...
    1:numel(compLabels), 'UniformOutput', false);
annotStr = ['TS=' num2str(totalSomatic) '-' 'TR=' num2str(totalResp) '-' strjoin(annotParts, '-')];

% Apply consistent font formatting
ax = findobj(fig, 'Type', 'axes');
for i = 1:numel(ax)
    formatAxes(ax(i), 8, 'helvetica');
end

% Set figure size and labels
figure(fig);
set(fig, 'Units', 'centimeters', 'Position', [20 20 5 6]);
ylabel('Responsive / Total responsive');
title('');

% % Shift axes up slightly to make room for the annotation
% pos    = get(gca, 'Position');
% pos(2) = pos(2) + 2;
% set(gca, 'Position', pos);

% Add bottom annotation with neuron counts
annotation('textbox', [0.1, 0.01, 0.8, 0.04], ...
    'String', annotStr, 'EdgeColor', 'none', ...
    'FontSize', 5, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FitBoxToText', false);

% Save figure if PaperFig mode is active
if params.PaperFig
    vs.printFig(fig, sprintf('ResponsiveUnits-%s', figTag), ...
        PaperFig = params.PaperFig);
end

if ~isempty(figAllDiffs)
    figure(figAllDiffs);
    % Save figure if PaperFig mode is active
    if params.PaperFig
        vs.printFig(figAllDiffs, sprintf('AllDiff-RespUnits-%s', figTag), ...
            PaperFig = params.PaperFig);
    end
end

% =========================================================================
% SECTION 11 — SAVE ANALYSIS STRUCT TO FIGURE DIRECTORY
% =========================================================================
% When PaperFig is true, save a companion .mat alongside the figures so
% that every figure folder is self-contained: figures + the exact analysis
% configuration that produced them.

if params.PaperFig
    % Retrieve the figure save directory from the analysis object
    % NOTE: adjust this path if vs.printFig uses a different convention
    rootPath = extractBefore(vs0.dataObj.recordingDir, 'lizards');  % 'W:\Large_scale_mapping_NP\'

    figSaveDir = [rootPath 'Paper_figs'];
    % Use the same base name as the analysis cache, with a suffix
    [~, cacheName, ~] = fileparts(nameOfFile);
    figStructPath = fullfile(figSaveDir, [cacheName '_analysisStruct.mat']);

    % Attach the summary statistics so the figure folder also carries the exact
    % numbers plotted (means, CIs, effect sizes) and not just the raw tables.
    % Safe to mutate S here: the analysis cache was already written in SECTION 5,
    % and this is the last use of S in the function.
    S.summaryStats = summaryStats;

    % Save the full struct (duplicates ~tens of KB; guarantees self-containment)
    save(figStructPath, '-struct', 'S');
    fprintf('Saved analysis struct to figure directory: %s\n', figStructPath);
end

end  % end function AllExpAnalysis


% #########################################################################
%                       LOCAL HELPER FUNCTIONS
% #########################################################################

function printStatsSummary(label, stats, pValRes)
% printStatsSummary  Console report of one metric's swarm-panel statistics.
%   Prints the per-stimulus block and the pairwise-difference block in a
%   paper-ready layout: observed central tendency with its bootstrap CI, then
%   the effect sizes with theirs. The effect size is labelled d for single
%   groups (one-sample Cohen's d vs zero) and d_z for differences (the paired
%   Cohen's d, since those rows hold per-unit differences).

    fprintf('\n=== %s: per-stimulus summary ===\n', label);
    printStatsBlock(stats.groups, 'd', pValRes);         % single-group effect size

    fprintf('\n=== %s: pairwise differences ===\n', label);
    printStatsBlock(stats.diffs,  'd_z', pValRes);       % paired effect size

end % printStatsSummary


function printStatsBlock(T, effectName, pValRes)
% printStatsBlock  Print one summary table, two lines per row.

    % Nothing computed for this metric (e.g. no pairs requested)
    if isempty(T)
        fprintf('  (no rows)\n');
        return
    end

    % The p column exists only on the differences table
    hasP = ismember('p', T.Properties.VariableNames);

    for i = 1:height(T)

        % Header line: label, sample sizes, and the pair p-value when present.
        % formatPValue keeps this identical to the figure label (plain decimal,
        % never scientific notation) so console and figure cannot disagree.
        if hasP
            [pNum, pOp] = formatPValue(T.p(i), pValRes);
            pStr = sprintf('   p %s %s', pOp, pNum);
        else
            pStr = '';
        end
        fprintf('  %-16s n=%-5d (%d animals, %d insertions)%s\n', ...
            T.group(i), T.n(i), T.nAnimals(i), T.nInsertions(i), pStr);

        % Central tendency: observed mean (the plotted bar), median, and the
        % hierarchical-bootstrap interval (the plotted patch).
        % bootMean is the mean of the bootstrap distribution and is NOT what is
        % drawn -- the bar is the observed mean. It is printed as a diagnostic:
        % a large gap signals unbalanced per-animal counts, because
        % hierBootMatchFreq pools every resampled observation flat, making the
        % resampled statistic a ratio with a random denominator whose
        % expectation drifts away from the observed mean.
        fprintf(['     mean %9.4f   median %9.4f   95%% CI [%9.4f, %9.4f]' ...
                 '   bootMean %9.4f (%+.4f)\n'], ...
            T.mean(i), T.median(i), T.ciLo(i), T.ciHi(i), ...
            T.bootMean(i), T.bootMean(i) - T.mean(i));

        % Effect sizes with their bootstrap intervals
        fprintf('     %-4s %8.3f  95%% CI [%7.3f, %7.3f]     P(>0) %.3f  95%% CI [%.3f, %.3f]\n', ...
            effectName, T.effect(i), T.effectLo(i), T.effectHi(i), ...
            T.propPos(i), T.propPosLo(i), T.propPosHi(i));
    end

end % printStatsBlock


function specOut = normalizeFocusSpec(specIn, nStim, defaultVal)
% normalizeFocusSpec  Broadcast/validate a per-stimulus focus spec cell.
    if isempty(specIn)
        specOut = repmat({defaultVal}, 1, nStim);        % default for every stimulus
    elseif numel(specIn) == 1
        specOut = repmat(specIn(1),    1, nStim);         % single entry -> broadcast
    else
        assert(numel(specIn) == nStim, ...                % else must be one per stimulus
            'Focus spec has %d entries but %d stimuli.', numel(specIn), nStim);
        specOut = specIn;
    end
end

function tf = isFocusActive(catLevels, lockSpec)
% isFocusActive  True if this stimulus uses non-default focus.
    isAll = isempty(catLevels) || (isscalar(string(catLevels)) && string(catLevels) == "all");
    tf    = ~isempty(lockSpec) || ~isAll;                 % lock and/or a real subset of factors
end


function gt = detectGratingType(stimName)
% detectGratingType  Auto-detect the GratingType parameter from stimulus
%   abbreviation.  Returns 'moving' for SDGm, 'static' for SDGs, or ''
%   for non-grating stimuli (in which case GratingType is not passed).
    switch stimName
        case 'SDGm', gt = 'moving';    % moving grating
        case 'SDGs', gt = 'static';    % static grating
        otherwise,   gt = '';           % not a grating stimulus
    end
end


function [labels, items] = buildMode3Items(stimsNeeded, catList, levelsCell)
% buildMode3Items  Build canonical comparison labels and (stim, cat, lvl) tuples.
%   labels{k} = 'MB_dir_0' etc.
%   items(k)  = struct('stim','MB', 'cat','direction', 'lv',0).

labels = {};
items  = struct('stim', {}, 'cat', {}, 'lv', {});

for si = 1:numel(stimsNeeded)
    sn   = stimsNeeded{si};
    cat  = catList{si};
    lvls = levelsCell{si};

    if isempty(cat)
        % No category → single item labeled with plain stimulus name
        labels{end+1} = sn;                                             %#ok<AGROW>
        items(end+1)  = struct('stim', sn, 'cat', '', 'lv', NaN);      %#ok<AGROW>
        continue
    end

    for lvi = 1:numel(lvls)
        lv = lvls(lvi);             % single level value
        labels{end+1} = makeCompLabel(sn, cat, lv);                     %#ok<AGROW>
        items(end+1)  = struct('stim', sn, 'cat', cat, 'lv', lv);      %#ok<AGROW>
    end
end
end


function lbl = makeCompLabel(stimName, catName, levelValue)
% makeCompLabel  Short composite label: '<stim>_<cat3>_<value>'.
%   Category truncated to 3 chars; decimals -> 'p'; negative -> 'neg'.

    catAbbr = lower(catName);                           % lowercase category
    if strlength(catAbbr) > 3
        catAbbr = extractBetween(catAbbr, 1, 3);       % truncate to 3 characters
        catAbbr = char(catAbbr);
    end
    lbl = sprintf('%s_%s_%g', stimName, catAbbr, levelValue);   % e.g. 'MB_dir_0'
    lbl = strrep(lbl, '.', 'p');                        % 0.3 -> 0p3
    lbl = strrep(lbl, '-', 'neg');                      % -1  -> neg1
end


function [vsObj, allFound] = findSessionWithLevels(NP, stimName, catName, requestedLevels, params)
% findSessionWithLevels  Find a session of stimName whose category column
%   contains ALL requested levels.  Tries Session=1 then Session=2.
%
%   ResponseWindow is recomputed with params.overwriteResponse before
%   reading colNames/C, to ensure stale cached column names are refreshed.

    vsObj    = [];
    allFound = false;

    for session = [1, 2]
        candidate = createStimulusObject(NP, stimName, session);    % try this session
        if isempty(candidate) || isempty(candidate.VST)
            continue                                                % session not available
        end

        % Recompute ResponseWindow (fixes stale column names if overwrite on)
        candidate.ResponseWindow( ...
            'overwrite',      params.overwriteResponse, ...
            'durationWindow', params.RespDurationWin);
        rw = candidate.ResponseWindow;

        % Extract the condition matrix and its column names
        [C, colNames] = getCmatrix(rw, stimName);
        if isempty(C) || isempty(colNames)
            continue
        end

        % Find the requested category column (case-insensitive)
        catIdx = find(strcmpi(colNames, catName));
        if isempty(catIdx)
            continue                                    % category not in this stimulus
        end

        catColIdx   = catIdx + 1;                       % +1 because colNames excludes first 4 cols
        availLevels = uniquetol(C(~isnan(C(:, catColIdx)), catColIdx), 1e-6);

        % Check that every requested level is present
        ok = true;
        for lv = requestedLevels(:)'
            if ~any(abs(availLevels - lv) < 1e-2)
                ok = false;  break
            end
        end

        if ok
            vsObj    = candidate;
            allFound = true;
            return
        end
    end
end


function [levels, catColIdx, vsObj] = findCategoryLevels(vsObj, NP, stimName, catName, params)
% findCategoryLevels  Find unique category levels in a recording (mode 2).
%   Tries Session=1, then Session=2.

    levels    = [];
    catColIdx = 0;

    for session = [1, 2]
        fprintf('  Trying Session=%d...\n', session);
        vsObj = createStimulusObject(NP, stimName, session);
        if isempty(vsObj) || isempty(vsObj.VST)
            continue
        end

        % Recompute ResponseWindow
        vsObj.ResponseWindow( ...
            'overwrite',      params.overwriteResponse, ...
            'durationWindow', params.RespDurationWin);
        rw = vsObj.ResponseWindow;

        % Extract condition matrix
        [C, colNames] = getCmatrix(rw, stimName);
        if isempty(C) || isempty(colNames)
            continue
        end

        % Look for the category column (case-insensitive)
        catIdx = find(strcmpi(colNames, catName));
        if isempty(catIdx)
            fprintf('  Category "%s" not found. Available: %s\n', ...
                catName, strjoin(colNames, ', '));
            return
        end

        catColIdx = catIdx + 1;                         % offset for first 4 metadata cols
        rawCol    = C(:, catColIdx);                    % raw values
        rawCol    = rawCol(~isnan(rawCol));             % remove NaNs
        levels    = uniquetol(rawCol, 1e-6);            % unique levels with tolerance

        if numel(levels) >= 2
            fprintf('  Found %d levels of "%s" (session %d): [%s]\n', ...
                numel(levels), catName, session, num2str(levels', '%.4g '));
            return                                      % success — exit early
        else
            fprintf('  Only %d level of "%s" in session %d.\n', ...
                numel(levels), catName, session);
        end
    end
end


function [C, colNames] = getCmatrix(rw, stimName)
% getCmatrix  Extract the condition matrix C and column names from a
%   ResponseWindow struct.  Returns empty if not available.

    C        = [];
    colNames = {};

    switch stimName
        case {'MB', 'MBR'}
            % MB/MBR store per-speed sub-structs; use the last (fastest) speed
            speedFields = fieldnames(rw);
            speedFields = speedFields(startsWith(speedFields, 'Speed'));
            if isempty(speedFields), return; end
            maxField = speedFields{end};
            C        = rw.(maxField).C;
            colNames = rw.colNames{1}(5:end);           % skip first 4 metadata columns

        case 'SDGm'
            if isfield(rw, 'C')
                C        = rw.C;
                colNames = rw.colNames{1}(5:end);
            end

        case 'SDGs'
            if isfield(rw, 'C')
                C        = rw.C;
                colNames = rw.colNames{1}(5:end);
            end

        otherwise
            % Generic fallback
            if isfield(rw, 'C')
                C        = rw.C;
                colNames = rw.colNames{1}(5:end);
            end
    end
end


function vsObj = createStimulusObject(NP, stimName, session)
% createStimulusObject  Create an analysis object, optionally with Session.
%   Returns [] on failure.

    vsObj = [];
    try
        key = getObjKey(stimName);      % shared key (e.g. 'SDG' for both SDGm/SDGs)
        if session == 0
            % No session specified — use default
            switch key
                case 'MB',  vsObj = linearlyMovingBallAnalysis(NP);
                case 'RG',  vsObj = rectGridAnalysis(NP);
                case 'MBR', vsObj = linearlyMovingBarAnalysis(NP);
                case 'SDG', vsObj = StaticDriftingGratingAnalysis(NP);
                case 'NI',  vsObj = imageAnalysis(NP);
                case 'NV',  vsObj = movieAnalysis(NP);
                case 'FFF', vsObj = fullFieldFlashAnalysis(NP);
            end
        else
            % Explicit session number
            switch key
                case 'MB',  vsObj = linearlyMovingBallAnalysis(NP, 'Session', session);
                case 'RG',  vsObj = rectGridAnalysis(NP, 'Session', session);
                case 'MBR', vsObj = linearlyMovingBarAnalysis(NP, 'Session', session);
                case 'SDG', vsObj = StaticDriftingGratingAnalysis(NP, 'Session', session);
                case 'NI',  vsObj = imageAnalysis(NP, 'Session', session);
                case 'NV',  vsObj = movieAnalysis(NP, 'Session', session);
                case 'FFF', vsObj = fullFieldFlashAnalysis(NP, 'Session', session);
            end
        end
    catch ME
        fprintf('  Could not create %s Session=%d: %s\n', stimName, session, ME.message);
        vsObj = [];
    end
end


function [vsObjs, present] = loadStimulusObjects(NP, stimsNeeded)
% loadStimulusObjects  Load one analysis object per unique stimulus class.
%   Returns a containers.Map of objects and a Map of presence flags.

    vsObjs  = containers.Map();     % key -> analysis object
    present = containers.Map();     % stimName -> true/false

    for si = 1:numel(stimsNeeded)
        sn  = stimsNeeded{si};      % stimulus name
        key = getObjKey(sn);        % shared object key

        if ~vsObjs.isKey(key)
            % First encounter of this key — create the object
            obj = createStimulusObject(NP, sn, 0);

            if isempty(obj) || isempty(obj.VST)
                fprintf('  %s: stimulus not found.\n', key);
                present(sn) = false;
            else
                present(sn) = true;
            end

            if ~isempty(obj)
                vsObjs(key) = obj;
            end
        else
            % Object already loaded for this key — check presence
            if ~present.isKey(sn)
                present(sn) = vsObjs.isKey(key) && ~isempty(vsObjs(key).VST);
            end
        end
    end
end


function key = getObjKey(stimName)
% getObjKey  Map stimulus abbreviation to shared analysis-object key.
%   SDGm and SDGs share a single StaticDriftingGratingAnalysis object.
    switch stimName
        case {'SDGm','SDGs'},  key = 'SDG';
        otherwise,             key = stimName;
    end
end


function fName = levelToFieldName(catName, value)
% levelToFieldName  Build a valid MATLAB field name matching
%   StatisticsPerNeuronPerCategory's convention.
%   e.g. ('size', 5) -> 'size_5', ('speed', 0.3) -> 'speed_0p3'.

    fName = sprintf('%s_%g', lower(strtrim(catName)), value);
    fName = strrep(fName, '.', 'p');    % decimal point -> 'p'
    fName = strrep(fName, '-', 'neg');  % minus sign    -> 'neg'
end


function stats = runStimStats(vsObj, params, catLevels, lockSpec)   % + two focus args
    assert(isa(vsObj, 'handle'), 'runStimStats:notHandle', ...
        'Analysis object must be a handle class for in-place mutation.');

    vsObj.ResponseWindow( ...                                % Step 1: response traces + condition matrix
        'overwrite',      params.overwriteResponse, ...
        'durationWindow', params.RespDurationWin);

    % BUG FIX: durationWindow was previously never passed here, so
    % StatisticsPerNeuron silently kept its own default (500 ms) regardless
    % of the value requested at the AllExpAnalysis level.
    %
    % BUG FIX: the result is now captured and returned so callers can pass it
    % straight into extractStimData. Previously extractStimData re-derived
    % stats via a bare vsObj.StatisticsPerNeuron() call (all arguments at
    % their defaults) — since StatisticsPerNeuron's cache guard compares ALL
    % computation params, that parameterless call would detect a mismatch
    % against the focused/tuned computation just performed here and silently
    % RECOMPUTE with defaults (CategoryLevels="all", Lock={}, durationWindow=
    % 500, ...), overwriting the correct cache and discarding the focus spec.
    stats = vsObj.StatisticsPerNeuron( ...                   % Step 2: max-permutation + t-test stats
        'overwrite',       params.overwriteStats, ...
        'BaseRespWindow',  params.BaseRespWindow, ...
        'SpatialGridMode', params.SpatialGridMode, ...
        'maxCategory',     params.maxCategory, ...
        'applyFDR',        params.applyFDR, ...
        'durationWindow',  params.durationWindow, ...       % sliding-window width for MovingWindowPVal (moving ball only)
        'CategoryLevels',  catLevels, ...                    % focus: factors kept separate
        'Lock',            lockSpec);                        % focus: locked factor values
end


function [z, p, spkR, spkDiff] = extractStimData(stats, vsObj, stimName, useZmean, useTTest)
% extractStimData  Pull z-scores, p-values, and spike rate from the
%   StatisticsPerNeuron results.
%
%   INPUTS
%     stats     struct  Result of runStimStats's vsObj.StatisticsPerNeuron(...)
%                       call. BUG FIX: this used to be re-derived here via a
%                       bare vsObj.StatisticsPerNeuron() (all defaults), which
%                       silently discarded CategoryLevels/Lock/durationWindow/
%                       etc. via a mismatched cache recompute. Callers must
%                       now pass through the exact struct runStimStats returned.
%     vsObj     Analysis object (used only for ResponseWindow, whose cache
%               has no params-match guard — reloading it here is safe)
%     stimName  char  Stimulus abbreviation ('MB','SDGm','SDGs', etc.)
%     useZmean  logical  true = use z_mean; false = use peak spike rate
%
%   OUTPUTS
%     z        (N,1) double   Z-scores per neuron
%     p        (N,1) double   P-values per neuron
%     spkR     (N,1) double   Spike rate (or z_mean) per neuron
%     spkDiff  (N,1) double   Response minus baseline difference

    rw = vsObj.ResponseWindow;   % safe: no params-match cache logic, just reloads the file runStimStats already wrote

    switch stimName

        case 'MB'
            % MB has multiple speeds — find best speed per neuron (lowest p)
            speedFields = fieldnames(stats);
            speedFields = speedFields(contains(speedFields, 'Speed'));
            nSpeeds = numel(speedFields);

            allP  = [];     % (nNeurons x nSpeeds) p-value matrix
            allZ  = [];     % (nNeurons x nSpeeds) z-score matrix
            allR  = [];     % (nNeurons x nSpeeds) spike rate matrix
            allDf = [];     % (nNeurons x nSpeeds) difference matrix

            for iS = 1:nSpeeds
                sName  = speedFields{iS};           % e.g. 'Speed1', 'Speed2'
                subTmp = stats.(sName);              % statistics sub-struct
                rwTmp  = rw.(sName);                 % response window sub-struct
                
                if useTTest
                    allP(:,iS) = subTmp.pValTTest(:);
                    
                else
                    allP(:,iS) = subTmp.pvalsResponse(:);
                end 
                
                
                allZ(:,iS) = subTmp.ZScoreU(:);                             %#ok<AGROW>

                if useZmean && isfield(subTmp, 'z_mean')
                    allR(:,iS) = subTmp.z_mean(:);                          %#ok<AGROW>
                else
                    allR(:,iS) = max(rwTmp.NeuronVals(:,:,4), [], 2);      %#ok<AGROW>
                end
                allDf(:,iS) = max(rwTmp.NeuronVals(:,:,5), [], 2);         %#ok<AGROW>
            end

            % Select the speed with the lowest p-value for each neuron
            [p, bestIdx] = min(allP, [], 2);
            nNeurons = size(allP, 1);
            linIdx   = sub2ind(size(allP), (1:nNeurons)', bestIdx);
            z        = allZ(linIdx);
            spkR     = allR(linIdx);
            spkDiff  = allDf(linIdx);
           
            return

        case 'MBR'
            sub = stats.Speed1;  rwSub = rw.Speed1;    % bar: single speed
        case 'SDGm'
            sub = stats.Moving;  rwSub = rw.Moving;    % moving grating sub-struct
        case 'SDGs'
            sub = stats.Static;  rwSub = rw.Static;    % static grating sub-struct
        otherwise
            sub = stats;         rwSub = rw;            % generic: top-level struct
    end

    % Extract per-neuron values from the selected sub-struct
    z = sub.ZScoreU(:);                                 % z-score

    if useTTest
        p = sub.pValTTest(:);                           % p-value
    else
        p = sub.pvalsResponse(:);                           % p-value
    end
   

    if useZmean && isfield(sub, 'z_mean')
        spkR = sub.z_mean(:);                           % z-scored mean response
    else
        spkR = max(rwSub.NeuronVals(:,:,4), [], 2);    % peak spike rate
    end

    spkDiff = max(rwSub.NeuronVals(:,:,5), [], 2);     % peak response - baseline
end




function fig = plotPairScatter(tbl, compLabels, metric, cmap, animalIdx, labelMap)
% plotPairScatter  Scatter plot: first comparison item (x) vs second (y).
%   Points coloured by animal identity.

    fig = figure;

    % Separate data for each stimulus
    mask1 = tbl.stimulus == compLabels{1};
    mask2 = tbl.stimulus == compLabels{2};
    v1    = tbl.(metric)(mask1);            % x-axis values
    v2    = tbl.(metric)(mask2);            % y-axis values
    cIdx  = animalIdx(mask1);               % colour index per point

    % Scatter with animal-coloured markers
    scatter(v1, v2, 7, cmap(cIdx,:), 'filled', 'MarkerFaceAlpha', 0.3);
    hold on;  axis equal;

    % Unity line (equal-response reference)
    lims = [min(tbl.(metric)), max(tbl.(metric))];
    plot(lims, lims, 'k--', 'LineWidth', 1.5);
    xlim(lims);  ylim(lims);

    % Apply display-name substitutions for axis labels
    xLab = compLabels{1};  yLab = compLabels{2};
    for li = 1:size(labelMap, 1)
        xLab = strrep(xLab, labelMap{li,1}, labelMap{li,2});
        yLab = strrep(yLab, labelMap{li,1}, labelMap{li,2});
    end
    xlabel(xLab);  ylabel(yLab);
    colormap(fig, cmap);

    % Consistent font formatting and figure size
    formatAxes(gca, 8, 'helvetica');
    set(fig, 'Units', 'centimeters', 'Position', [20 20 5 5]);
end


function formatAxes(ax, fontSize, fontName)
% formatAxes  Apply consistent font styling to an axes object.
    ax.YAxis.FontSize = fontSize;  ax.YAxis.FontName = fontName;
    ax.XAxis.FontSize = fontSize;  ax.XAxis.FontName = fontName;
end


function pAdj = bhFDR(pVals)
% bhFDR  Benjamini-Hochberg FDR correction.
%   Adjusts a vector of p-values to control the false discovery rate.

    n = numel(pVals);                           % number of tests
    [pSorted, sortIdx] = sort(pVals(:));        % sort ascending
    ranks = (1:n)';                             % rank of each sorted p-value

    pAdj = pSorted .* n ./ ranks;               % BH adjustment: p * n / rank
    pAdj = min(pAdj, 1);                        % cap at 1

    % Enforce monotonicity: each adjusted p must be <= the one above it
    for k = n-1:-1:1
        pAdj(k) = min(pAdj(k), pAdj(k+1));
    end

    pAdj(sortIdx) = pAdj;                       % restore original order
end