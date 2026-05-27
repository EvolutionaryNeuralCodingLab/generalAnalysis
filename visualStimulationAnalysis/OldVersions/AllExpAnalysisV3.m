function [tempTable] = AllExpAnalysis(expList, params)
% AllExpAnalysis  Pool neural responses across Neuropixels recordings,
%   run pairwise statistical comparisons via hierarchical bootstrapping,
%   and generate publication-ready swarm and scatter plots.
%
%   This function:
%     1. Iterates over a list of experiments, loading pre-computed per-neuron
%        statistics (z-scores, p-values, spike rates) for each stimulus.
%     2. For each recording, identifies neurons responsive to ANY stimulus
%        in ComparePairs (OR union) and adds them to a long-format table.
%     3. Computes pairwise hierarchical bootstrap tests between stimuli.
%     4. Plots swarm charts and scatter plots for z-scores and spike rates.
%     5. Computes a fraction-responsive comparison across insertions.
%
%   INPUTS
%     expList  (1,:) double   Row vector of experiment IDs from master Excel.
%     params   Name-value     See arguments block below.
%
%   OUTPUTS
%     tempTable  table   Fraction-responsive table (one row per insertion ×
%                        stimulus), filtered to insertions containing all
%                        compared stimuli.
%
%   EXAMPLE
%     tempTable = AllExpAnalysis([49:54 64:66], ...
%         ComparePairs = {'SDGm','SDGs'}, ...
%         StatMethod   = 'maxPermuteTest', ...
%         PaperFig     = true);
%
%   See also: hierBoot, plotSwarmBootstrapWithComparisons

% =========================================================================
%                        ARGUMENTS BLOCK
% =========================================================================
arguments
    expList (1,:) double                        % Experiment IDs from master Excel
    params.ComparePairs cell                    % Stimuli to compare, e.g. {'SDGm','SDGs'}
                                                %   Neurons significant for ANY of these
                                                %   are included.  Statistics are pairwise.
    params.threshold        double  = 0.05      % p-value cutoff for responsiveness
    params.StatMethod       string  = 'ObsWindow'  % 'ObsWindow' | 'bootsrapRespBase' | 'maxPermuteTest'
    params.overwrite        logical = false      % Recompute and overwrite saved pooled file
    params.overwriteResponse logical = false     % Force re-run of ResponseWindow
    params.overwriteStats   logical = false      % Force re-run of per-neuron statistics
    params.RespDurationWin  double  = 100        % Response window duration (ms)
    params.shuffles         double  = 2000       % Shuffles / bootstrap iterations for per-neuron stats
    params.useZmean         logical = true       % Use z_mean (response−baseline normalised by null)
                                                 %   instead of raw peak spike rate
    params.useFDR           logical = false       % Apply Benjamini-Hochberg FDR correction per recording
    params.PaperFig         logical = false       % Save figures via vs.printFig
    params.nBoot            double  = 10000       % Bootstrap iterations for group-level tests
end

% =========================================================================
% SECTION 1 — SETUP: DETERMINE STIMULI, PATHS, AND CACHING
% =========================================================================

% Unique stimulus names that need to be loaded and analysed
stimsNeeded = unique(params.ComparePairs, 'stable');  % e.g. {'SDGm','SDGs'}

% Number of experiments to process
nExp = numel(expList);

% Load the first experiment to extract directory paths
NP0 = loadNPclassFromTable(expList(1));            % Neuropixels recording object
vs0 = linearlyMovingBallAnalysis(NP0);             % analysis object (used for path only)

% Build the output directory: <root>/lizards/Combined_lizard_analysis/
rootPath = extractBefore(vs0.getAnalysisFileName, 'lizards');  % path up to 'lizards'
rootPath = [rootPath 'lizards'];                                % append 'lizards'
saveDir  = fullfile(rootPath, 'Combined_lizard_analysis');      % subdirectory for pooled results
if ~exist(saveDir, 'dir')                                       % create if absent
    mkdir(saveDir);
end

% Construct a descriptive filename for the cached pooled data
nameOfFile = sprintf('Ex_%d-%d_Combined_%s.mat', ...
    expList(1), expList(end), strjoin(stimsNeeded, '-'));

savePath = fullfile(saveDir, nameOfFile);  % full path to cached .mat file

% Decide whether the per-experiment loop needs to run:
%   Skip if a cached file exists with the same experiment list AND overwrite=false
runLoop = true;                                          % default: run the loop
if exist(savePath, 'file') == 2 && ~params.overwrite     % cached file found
    S = load(savePath);                                  % load cached struct
    if isfield(S, 'expList') && isequal(S.expList, expList)
        runLoop = false;                                 % cache is valid → skip loop
    end
end

% =========================================================================
% SECTION 2 — INITIALISE LONG-FORMAT TABLES
% =========================================================================

% TableStimComp: one row per neuron × stimulus.
%   Contains z-scores and spike rates for neurons responsive to ANY stimulus
%   in ComparePairs (OR union across all stimuli).
TableStimComp = table( ...
    categorical.empty(0,1), ...   % animal   — animal ID (e.g. 'PV97')
    categorical.empty(0,1), ...   % insertion — insertion counter (unique per probe track)
    categorical.empty(0,1), ...   % stimulus  — stimulus abbreviation (e.g. 'SDGm')
    categorical.empty(0,1), ...   % NeurID    — unit index within the recording
    double.empty(0,1), ...        % Z-score   — z-score for this neuron × stimulus
    double.empty(0,1), ...        % SpkR      — spike rate (or z_mean) for this neuron × stimulus
    'VariableNames', {'animal','insertion','stimulus','NeurID','Z-score','SpkR'});

% TableRespNeurs: one row per insertion × stimulus.
%   Counts how many neurons are responsive to each stimulus (self-significant),
%   and the total number of somatic units in that recording.
TableRespNeurs = table( ...
    categorical.empty(0,1), ...   % animal
    categorical.empty(0,1), ...   % insertion
    categorical.empty(0,1), ...   % stimulus
    double.empty(0,1), ...        % respNeur       — count of responsive neurons
    double.empty(0,1), ...        % totalSomaticN  — total sorted units in recording
    'VariableNames', {'animal','insertion','stimulus','respNeur','totalSomaticN'});

% =========================================================================
% SECTION 3 — PER-EXPERIMENT LOOP
% =========================================================================

if runLoop

    % Counters for unique animals and insertions across the experiment list
    animalCount    = 0;          % running count of distinct animals
    insertionCount = 0;          % running count of distinct probe insertions
    prevAnimal     = "";         % animal ID from the previous iteration
    prevInsertion  = 0;          % insertion number from the previous iteration

    j = 1;   % 1-based experiment counter (indexes cell arrays if needed later)

    for ex = expList  % ---- iterate over each experiment ID ----

        % ------------------------------------------------------------------
        % 3a — Load the recording and check stimulus availability
        % ------------------------------------------------------------------

        NP = loadNPclassFromTable(ex);                % load Neuropixels recording object
        fprintf('Processing recording: %s\n', NP.recordingName);  % status message

        % Load analysis objects and check which stimuli are present
        %   vsObjs  — containers.Map: objKey → analysis object
        %   present — containers.Map: stimName → logical
        [vsObjs, present] = loadStimulusObjects(NP, stimsNeeded);

        % Skip this experiment if ANY needed stimulus is absent
        allPresent = true;                              % assume all present
        for si = 1:numel(stimsNeeded)
            if ~present(stimsNeeded{si})
                allPresent = false;                     % at least one missing
                break
            end
        end
        if ~allPresent
            fprintf('  → Skipping: not all stimuli present.\n');
            continue                                    % skip to next experiment
        end

        % ------------------------------------------------------------------
        % 3b — Parse metadata and update animal / insertion counters
        %      (only reached if all stimuli are present)
        % ------------------------------------------------------------------

        % Extract animal ID from the recording name (expects 'PV##' or 'SA##')
        animalID = string(regexp(NP.recordingName, 'PV\d+', 'match', 'once'));
        if animalID == ""                                   % fallback naming convention
            animalID = string(regexp(NP.recordingName, 'SA\d+', 'match', 'once'));
        end

        % Extract insertion number from filename (e.g. 'Insertion2' → 2)
        insStr = regexp( ...
            linearlyMovingBallAnalysis(NP).getAnalysisFileName, ...
            'Insertion\d+', 'match', 'once');              % match 'Insertion#'
        insNum = str2double(regexp(insStr, '\d+', 'match'));% parse the digit(s)

        % Update animal and insertion counters
        animalChanged = (animalID ~= prevAnimal);           % new animal?
        if animalChanged
            animalCount = animalCount + 1;                  % increment animal counter
            prevAnimal  = animalID;                         % update tracker
        end
        if insNum ~= prevInsertion || animalChanged         % new insertion?
            insertionCount = insertionCount + 1;            % increment insertion counter
            prevInsertion  = insNum;                        % update tracker
        end

        % ------------------------------------------------------------------
        % 3c — Run ResponseWindow + statistics for each stimulus
        % ------------------------------------------------------------------

        objKeys = keys(vsObjs);                % unique analysis-object keys
        for k = 1:numel(objKeys)
            key   = objKeys{k};                % e.g. 'SDG', 'MB', 'RG'
            vsObj = vsObjs(key);               % the analysis object
            runStimStats(vsObj, params);        % ResponseWindow + chosen stat method
            vsObjs(key) = vsObj;               % store back (handle class, but explicit)
        end

        % ------------------------------------------------------------------
        % 3d — Extract z-scores, p-values, spike rates for each stimulus
        %      All stimuli are guaranteed present at this point.
        % ------------------------------------------------------------------

        stimData = struct();                            % one sub-struct per stimulus
        nUnits   = [];                                  % total unit count (set from first stim)

        for si = 1:numel(stimsNeeded)
            sn  = stimsNeeded{si};                      % stimulus name (e.g. 'SDGm')
            key = getObjKey(sn);                        % shared-object key (e.g. 'SDG')

            [z, p, spkR, spkDiff] = extractStimData( ...
                vsObjs(key), sn, params.StatMethod, params.useZmean);
            stimData.(sn).z       = z(:);               % force column vector
            stimData.(sn).p       = p(:);
            stimData.(sn).spkR    = spkR(:);
            stimData.(sn).spkDiff = spkDiff(:);

            if isempty(nUnits)
                nUnits = numel(z);                      % set unit count from first stimulus
            end
        end

        % ------------------------------------------------------------------
        % 3e — Optional: Benjamini-Hochberg FDR correction per recording
        % ------------------------------------------------------------------

        if params.useFDR
            for si = 1:numel(stimsNeeded)
                sn = stimsNeeded{si};
                stimData.(sn).p = bhFDR(stimData.(sn).p);  % adjust p-values for FDR
            end
        end

        % ------------------------------------------------------------------
        % 3f — Build OR significance mask across all compared stimuli
        %      A neuron is included if it is significant for ANY stimulus
        %      in ComparePairs.
        % ------------------------------------------------------------------

        orMask = false(nUnits, 1);                      % initialise all-false mask
        for si = 1:numel(stimsNeeded)
            sn = stimsNeeded{si};
            orMask = orMask | (stimData.(sn).p < params.threshold);  % OR with each stimulus
        end

        unitIDs = find(orMask);                          % indices of neurons passing the OR filter
        nSig    = numel(unitIDs);                        % count of significant neurons

        % ------------------------------------------------------------------
        % 3g — Append rows to TableStimComp (neuron-level pairwise table)
        %      Each significant neuron gets one row PER stimulus.
        % ------------------------------------------------------------------

        if nSig > 0
            for si = 1:numel(stimsNeeded)
                sn = stimsNeeded{si};

                % Build a mini-table for this stimulus × this recording
                newRows = table( ...
                    repmat(categorical(cellstr(animalID)), nSig, 1), ...  % animal column
                    repmat(categorical(insertionCount),    nSig, 1), ...  % insertion column
                    repmat(categorical(cellstr(sn)),       nSig, 1), ...  % stimulus column
                    categorical(unitIDs), ...                             % neuron ID column
                    stimData.(sn).z(orMask), ...                          % z-score column
                    stimData.(sn).spkR(orMask), ...                       % spike rate column
                    'VariableNames', {'animal','insertion','stimulus','NeurID','Z-score','SpkR'});

                TableStimComp = [TableStimComp; newRows];   % append to pooled table
            end
        end

        % ------------------------------------------------------------------
        % 3h — Append rows to TableRespNeurs (insertion-level counts)
        %      One row per stimulus: how many neurons respond to THIS stimulus
        %      (self-significant, not OR union), and total unit count.
        % ------------------------------------------------------------------

        for si = 1:numel(stimsNeeded)
            sn      = stimsNeeded{si};
            nResp   = sum(stimData.(sn).p < params.threshold);  % self-responsive count
            newRow = table( ...
                categorical(cellstr(animalID)), ...              % animal
                categorical(insertionCount), ...                 % insertion
                categorical(cellstr(sn)), ...                    % stimulus
                nResp, ...                                       % respNeur
                nUnits, ...                                      % totalSomaticN
                'VariableNames', TableRespNeurs.Properties.VariableNames);
            TableRespNeurs = [TableRespNeurs; newRow];           % append row
        end

        fprintf('  → %d / %d units pass OR filter.\n', nSig, nUnits);
        j = j + 1;   % advance experiment counter

    end  % ---- end for ex = expList ----

    % =====================================================================
    % SECTION 4 — SAVE POOLED DATA
    % =====================================================================

    S.expList        = expList;         % experiment IDs that were processed
    S.TableStimComp  = TableStimComp;   % neuron-level pairwise table
    S.TableRespNeurs = TableRespNeurs;  % insertion-level responsive counts
    S.params         = params;          % parameter snapshot for reproducibility

    save(savePath, '-struct', 'S');     % save struct fields as top-level vars
    fprintf('Saved pooled data to %s\n', savePath);

end  % end if runLoop

% =========================================================================
% SECTION 5 — GUARD: ABORT EARLY IF NO SIGNIFICANT NEURONS WERE FOUND
% =========================================================================

if isempty(S.TableStimComp) || height(S.TableStimComp) == 0
    warning('AllExpAnalysis:noUnits', ...
        'No significant units found for comparison of %s. Returning empty.', ...
        strjoin(stimsNeeded, ' vs '));
    tempTable = table();              % return empty table
    return
end

% Replace any residual NaN z-scores or spike rates with 0
%   (conservative: treat NaN as "no response" for bootstrap)
S.TableStimComp.('Z-score')(isnan(S.TableStimComp.('Z-score'))) = 0;
S.TableStimComp.SpkR(isnan(S.TableStimComp.SpkR))               = 0;

% =========================================================================
% SECTION 6 — SHARED PLOTTING SETUP
% =========================================================================

% Reload an analysis object for figure-saving paths
NP = loadNPclassFromTable(expList(1));
vs = linearlyMovingBallAnalysis(NP);

% Build a shared colormap so every animal gets the same colour across all panels
animalOrder = categories(S.TableStimComp.animal);   % canonical alphabetical ordering
nAnimals    = numel(animalOrder);                   % number of distinct animals
sharedCmap  = lines(nAnimals);                      % nAnimals × 3 RGB matrix

% Numeric animal index for each row (used for colour lookup)
animalIdxAll = double(S.TableStimComp.animal);

% Generate all pairwise combinations of stimuli for statistical testing
%   e.g. {'SDGm','SDGs'} → one pair; {'MB','RG','MBR'} → three pairs
pairsAll = nchoosek(stimsNeeded, 2);        % nPairs × 2 cell array of pairs

% Label-replacement map for display (internal abbreviation → paper label)
labelMap = {'RG','SB'; 'SDGs','SG'; 'SDGm','MG'};

% =========================================================================
% SECTION 7 — Z-SCORE PAIRWISE COMPARISON
% =========================================================================

% --- 7a: Hierarchical bootstrap for each pair ---

pValsZ = zeros(1, size(pairsAll, 1));       % one p-value per pair

for pi = 1:size(pairsAll, 1)               % iterate over stimulus pairs
    % Compute per-neuron Z-score differences and run hierarchical bootstrap
    pValsZ(pi) = bootstrapPairDifference( ...
        S.TableStimComp, pairsAll(pi,:), params.nBoot, 'Z-score');
end

% --- 7b: Swarm plot of Z-scores ---

ZscoreYlim = ceil(max(S.TableStimComp.('Z-score'))) + 4;  % y-axis ceiling

fig = plotSwarmBootstrapWithComparisons( ...
    S.TableStimComp, pairsAll, pValsZ, {'Z-score'}, ...
    yLegend    = 'Z-score', ...
    yMaxVis    = ZscoreYlim, ...
    diff       = true, ...
    plotMeanSem = true, ...
    Alpha      = 0.7);

formatAxes(gca, 8, 'helvetica');                          % consistent font styling
set(fig, 'Units', 'centimeters', 'Position', [20 20 10 6]);
colormap(fig, sharedCmap);                                 % enforce shared colour scheme

if params.PaperFig
    vs.printFig(fig, sprintf('Zscore-Swarm-%s', strjoin(stimsNeeded,'-')), ...
        PaperFig = params.PaperFig);
end

% --- 7c: Scatter plot — first stimulus vs second stimulus (Z-score) ---

if numel(stimsNeeded) == 2
    fig = plotPairScatter(S.TableStimComp, stimsNeeded, ...
        'Z-score', sharedCmap, animalIdxAll, labelMap);
    title('Z-score');

    if params.PaperFig
        vs.printFig(fig, sprintf('Zscore-Scatter-%s', strjoin(stimsNeeded,'-')), ...
            PaperFig = params.PaperFig);
    end
end

% =========================================================================
% SECTION 8 — SPIKE-RATE PAIRWISE COMPARISON
% =========================================================================

% --- 8a: Hierarchical bootstrap for each pair ---

pValsSpk = zeros(1, size(pairsAll, 1));

for pi = 1:size(pairsAll, 1)
    pValsSpk(pi) = bootstrapPairDifference( ...
        S.TableStimComp, pairsAll(pi,:), params.nBoot, 'SpkR');
end

% --- 8b: Swarm plot of spike rates ---

spkMax = max(S.TableStimComp.SpkR);                        % y-axis ceiling

fig = plotSwarmBootstrapWithComparisons( ...
    S.TableStimComp, pairsAll, pValsSpk, {'SpkR'}, ...
    yLegend    = 'SpkR', ...
    yMaxVis    = spkMax, ...
    diff       = true, ...
    plotMeanSem = true, ...
    Alpha      = 0.7);

formatAxes(gca, 8, 'helvetica');
colormap(fig, sharedCmap);
set(fig, 'Units', 'centimeters', 'Position', [20 20 10 6]);

if params.PaperFig
    vs.printFig(fig, sprintf('SpkRate-Swarm-%s', strjoin(stimsNeeded,'-')), ...
        PaperFig = params.PaperFig);
end

% --- 8c: Scatter plot — first stimulus vs second stimulus (spike rate) ---

if numel(stimsNeeded) == 2
    fig = plotPairScatter(S.TableStimComp, stimsNeeded, ...
        'SpkR', sharedCmap, animalIdxAll, labelMap);
    title('Spk. rate');

    if params.PaperFig
        vs.printFig(fig, sprintf('SpkRate-Scatter-%s', strjoin(stimsNeeded,'-')), ...
            PaperFig = params.PaperFig);
    end
end

% =========================================================================
% SECTION 9 — FRACTION-RESPONSIVE ANALYSIS
%   Compares the proportion of responsive neurons between stimuli,
%   bootstrapping at the insertion level (no hierarchy needed because
%   there is one data point per insertion).
% =========================================================================

% Find insertions that contain ALL compared stimuli
[G, ~] = findgroups(S.TableRespNeurs.insertion);       % group by insertion
hasAll  = splitapply( ...                               % check each group
    @(s) all(ismember(categorical(stimsNeeded), s)), ...%   does it contain every stimulus?
    S.TableRespNeurs.stimulus, G);

% Restrict to complete insertions and relevant stimuli only
tempTable = S.TableRespNeurs( ...
    hasAll(G) & ismember(S.TableRespNeurs.stimulus, categorical(stimsNeeded)), :);

% Bootstrap the difference in responsive fraction for each pair
pValsFrac = zeros(1, size(pairsAll, 1));

for pi = 1:size(pairsAll, 1)

    diffs = [];   % will hold one fraction-difference per insertion

    for ins = unique(S.TableRespNeurs.insertion)'        % iterate over insertions

        % Find rows for this insertion × each stimulus in the pair
        idx1 = S.TableRespNeurs.insertion == categorical(ins) & ...
               S.TableRespNeurs.stimulus  == pairsAll{pi,1};
        idx2 = S.TableRespNeurs.insertion == categorical(ins) & ...
               S.TableRespNeurs.stimulus  == pairsAll{pi,2};

        if any(idx1) && any(idx2)
            total = S.TableRespNeurs.totalSomaticN(idx1);       % shared denominator
            f1    = S.TableRespNeurs.respNeur(idx1) / total;    % fraction responsive stim1
            f2    = S.TableRespNeurs.respNeur(idx2) / total;    % fraction responsive stim2
            diffs(end+1, 1) = f1 - f2;                          % per-insertion difference
        end
    end

    % Simple bootstrap of the mean difference (one value per insertion → flat)
    bootDiff      = bootstrp(params.nBoot, @mean, diffs);  % nBoot × 1 bootstrap means
    pValsFrac(pi) = mean(bootDiff <= 0);                   % p-value: prop ≤ 0
end

% Add a total-responsive column (sum across stimuli within each insertion)
[G, ~]                  = findgroups(tempTable.insertion);
totals                  = splitapply(@sum, tempTable.respNeur, G);
tempTable.TotalRespNeur = totals(G);

% Plot fraction-responsive with significance annotation
fig = plotSwarmBootstrapWithComparisons( ...
    tempTable, pairsAll, pValsFrac, ...
    {'respNeur','totalSomaticN'}, ...
    fraction       = true, ...
    showBothAndDiff = false, ...
    yLegend        = 'Responsive/total units', ...
    diff           = false, ...
    filled         = false, ...
    Xjitter        = 'none', ...
    Alpha          = 0.6, ...
    drawLines      = true);

% Compute summary counts for the annotation
totalResp  = sum(tempTable.respNeur);                                     % all stims combined
perStimN   = arrayfun(@(s) sum(tempTable.respNeur(tempTable.stimulus == s)), ...
                      categorical(stimsNeeded));                          % per-stimulus counts

% Build annotation string: 'TR = 45 - SDGm = 28 - SDGs = 17'
annotParts = arrayfun(@(i) sprintf('%s = %d', stimsNeeded{i}, perStimN(i)), ...
                      1:numel(stimsNeeded), 'UniformOutput', false);
annotStr   = ['TR = ' num2str(totalResp) ' - ' strjoin(annotParts, ' - ')];

formatAxes(gca, 8, 'helvetica');
set(fig, 'Units', 'centimeters', 'Position', [20 20 5 6]);
ylabel('Responsive / Total responsive');
title('');

% Shift axes up slightly to make room for bottom annotation
pos    = get(gca, 'Position');       % [left bottom width height]
pos(2) = pos(2) + 0.05;             % push bottom edge up
set(gca, 'Position', pos);

% Place annotation at the bottom of the figure
annotation('textbox', [0.1, 0.01, 0.8, 0.04], ...
    'String',              annotStr, ...
    'EdgeColor',           'none', ...
    'FontSize',            9, ...
    'FontWeight',          'bold', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment',   'middle', ...
    'FitBoxToText',        false);

if params.PaperFig
    vs.printFig(fig, sprintf('ResponsiveUnits-%s', strjoin(stimsNeeded,'-')), ...
        PaperFig = params.PaperFig);
end

end  % end function AllExpAnalysis


% #########################################################################
%                       LOCAL HELPER FUNCTIONS
% #########################################################################

function [vsObjs, present] = loadStimulusObjects(NP, stimsNeeded)
% loadStimulusObjects  Load one analysis object per unique stimulus class.
%
%   Several stimuli (e.g. SDGm and SDGs) share the same analysis object.
%   This function loads each object at most once, and checks whether each
%   stimulus was actually recorded by inspecting the VST property.
%
%   INPUTS
%     NP          Neuropixels recording object
%     stimsNeeded cell array of stimulus abbreviations
%
%   OUTPUTS
%     vsObjs   containers.Map  objKey → analysis object
%     present  containers.Map  stimName → logical (true if recorded)

    vsObjs  = containers.Map();          % cache of loaded analysis objects
    present = containers.Map();          % presence flag per stimulus

    for si = 1:numel(stimsNeeded)
        sn  = stimsNeeded{si};           % stimulus name
        key = getObjKey(sn);             % shared object key

        % Load the analysis object if not already cached
        if ~vsObjs.isKey(key)
            try
                switch key
                    case 'MB',  obj = linearlyMovingBallAnalysis(NP);
                    case 'RG',  obj = rectGridAnalysis(NP);
                    case 'MBR', obj = linearlyMovingBarAnalysis(NP);
                    case 'SDG', obj = StaticDriftingGratingAnalysis(NP);
                    case 'NI',  obj = imageAnalysis(NP);
                    case 'NV',  obj = movieAnalysis(NP);
                    case 'FFF', obj = fullFieldFlashAnalysis(NP);
                end

                % Check if the stimulus was actually presented
                if isempty(obj.VST)
                    fprintf('  %s: stimulus not found in recording.\n', key);
                    present(sn) = false;    % VST empty → not recorded
                else
                    present(sn) = true;     % VST populated → was recorded
                end

                vsObjs(key) = obj;          % cache the object

            catch ME
                fprintf('  %s: could not load (%s).\n', key, ME.message);
                present(sn) = false;        % constructor failed → not present
            end
        else
            % Object already loaded; still need to set presence for this stimulus name
            % (e.g. SDGm present ≠ SDGs present → both use the same object,
            %  but both are present if the object loaded successfully)
            if ~present.isKey(sn)
                % If the shared object loaded with non-empty VST, mark present
                if isKey(vsObjs, key) && ~isempty(vsObjs(key).VST)
                    present(sn) = true;
                else
                    present(sn) = false;
                end
            end
        end
    end
end


function key = getObjKey(stimName)
% getObjKey  Map a stimulus abbreviation to its analysis-object key.
%   SDGm and SDGs both map to 'SDG' because they share one object.

    switch stimName
        case {'SDGm','SDGs'},  key = 'SDG';
        otherwise,             key = stimName;    % MB, RG, MBR, NI, NV, FFF
    end
end


function runStimStats(vsObj, params)
% runStimStats  Run ResponseWindow and the chosen statistical method.
%
%   Dispatches to ShufflingAnalysis, BootstrapPerNeuron, or
%   StatisticsPerNeuron depending on params.StatMethod.

    % Compute or load the response window
    vsObj.ResponseWindow( ...
        'overwrite',       params.overwriteResponse, ...
        'durationWindow',  params.RespDurationWin);

    % Run the chosen statistical method
    switch params.StatMethod
        case 'ObsWindow'
            vsObj.ShufflingAnalysis( ...
                'overwrite',    params.overwriteStats, ...
                'N_bootstrap',  params.shuffles);

        case 'bootsrapRespBase'
            vsObj.BootstrapPerNeuron('overwrite', params.overwriteStats);

        case 'maxPermuteTest'
            vsObj.StatisticsPerNeuron('overwrite', params.overwriteStats);

        otherwise
            error('AllExpAnalysis:badMethod', ...
                'Unknown StatMethod "%s".', params.StatMethod);
    end
end


function [z, p, spkR, spkDiff] = extractStimData(vsObj, stimName, statMethod, useZmean)
% extractStimData  Pull z-scores, p-values, spike rate, and spike-rate
%   difference from a stats struct, navigating the stimulus-specific nesting.
%
%   The struct layout varies by stimulus type:
%     Flat:   stats.ZScoreU           (RG, FFF, NI, NV)
%     Speed:  stats.Speed1.ZScoreU    (MB prefers Speed2; MBR uses Speed1)
%     Moving/Static: stats.Moving.*   (SDGm) or stats.Static.* (SDGs)

    % --- Retrieve the stats struct (dispatch on statistical method) ---
    switch statMethod
        case 'ObsWindow',       stats = vsObj.ShufflingAnalysis;
        case 'bootsrapRespBase', stats = vsObj.BootstrapPerNeuron;
        case 'maxPermuteTest',  stats = vsObj.StatisticsPerNeuron;
    end

    rw = vsObj.ResponseWindow;   % response-window struct (spike rates stored here)

    % --- Navigate to the correct sub-struct for this stimulus ---
    switch stimName
        case 'MB'
            % MB has Speed1 and optionally Speed2 (faster, more salient).
            % Prefer Speed2 for z-scores/p-values; spike rate from Speed1.
            sub   = stats.Speed1;
            rwSub = rw.Speed1;
            if isfield(stats, 'Speed2')
                sub   = stats.Speed2;           % z-scores/p from faster speed
                % NOTE: spike rate intentionally comes from Speed1 (original convention)
                rwSub = rw.Speed1;
            end

        case 'MBR'
            sub   = stats.Speed1;               % moving bar: Speed1 only
            rwSub = rw.Speed1;

        case 'SDGm'
            sub   = stats.Moving;               % drifting gratings: moving condition
            rwSub = rw.Moving;

        case 'SDGs'
            sub   = stats.Static;               % drifting gratings: static condition
            rwSub = rw.Static;

        otherwise   % RG, FFF, NI, NV — flat struct
            sub   = stats;
            rwSub = rw;
    end

    % --- Extract z-scores and p-values ---
    z = sub.ZScoreU(:);                         % force column vector
    p = sub.pvalsResponse(:);

    % --- Extract spike rate ---
    if useZmean && isfield(sub, 'z_mean')
        spkR = sub.z_mean(:);                   % normalised response (z_mean)
    else
        spkR = max(rwSub.NeuronVals(:,:,4), [], 2);  % peak rate across directions
    end

    % --- Extract spike-rate difference (response – baseline) ---
    spkDiff = max(rwSub.NeuronVals(:,:,5), [], 2);

    % --- Override spike rate for bootstrap method (uses observed responses) ---
    if strcmp(statMethod, 'bootsrapRespBase') && isfield(sub, 'ObsResponse')
        spkR = mean(sub.ObsResponse, 1)';       % mean across repeats
    end
end


function pVal = bootstrapPairDifference(tbl, pair, nBoot, metric)
% bootstrapPairDifference  Hierarchical bootstrap test for a single pair.
%
%   Computes per-neuron differences (stim1 − stim2), then resamples at the
%   animal → insertion → neuron hierarchy.
%
%   INPUTS
%     tbl     table      Long-format table with columns: insertion, stimulus,
%                        animal, and the metric column.
%     pair    {1×2} cell Stimulus pair, e.g. {'SDGm','SDGs'}.
%     nBoot   double     Number of bootstrap iterations.
%     metric  char       Column name to compare ('Z-score' or 'SpkR').
%
%   OUTPUT
%     pVal    double     Proportion of bootstrap means ≤ 0.

    diffs   = [];           % per-neuron differences pooled across insertions
    insers  = [];           % insertion label for each difference
    animals = [];           % animal label for each difference

    for ins = unique(tbl.insertion)'                     % iterate over insertions

        % Select rows: this insertion × each stimulus in the pair
        idx1 = tbl.insertion == categorical(ins) & tbl.stimulus == pair{1};
        idx2 = tbl.insertion == categorical(ins) & tbl.stimulus == pair{2};

        V1 = tbl.(metric)(idx1);                         % metric values for stim 1
        V2 = tbl.(metric)(idx2);                         % metric values for stim 2

        if isempty(V1) || isempty(V2)
            continue                                      % skip incomplete insertions
        end

        animal = unique(tbl.animal(idx1));                % animal for this insertion

        diffs   = [diffs;   V1 - V2];                    %#ok<AGROW> append differences
        insers  = [insers;  double(repmat(ins,    numel(V1), 1))]; %#ok<AGROW>
        animals = [animals; double(repmat(animal, numel(V1), 1))]; %#ok<AGROW>
    end

    % Run hierarchical bootstrap (resample animals → insertions → neurons)
    bootMeans = hierBoot(diffs, nBoot, insers, animals);
    pVal      = mean(bootMeans <= 0);                     % one-sided p-value
end


function fig = plotPairScatter(tbl, stimsNeeded, metric, cmap, animalIdx, labelMap)
% plotPairScatter  Scatter the first vs second stimulus for a given metric.
%
%   Each dot is one neuron; colour = animal identity.

    fig = figure;

    % Extract data for each stimulus
    mask1 = tbl.stimulus == stimsNeeded{1};               % rows for stimulus 1
    mask2 = tbl.stimulus == stimsNeeded{2};               % rows for stimulus 2
    v1    = tbl.(metric)(mask1);                          % metric values for stim 1
    v2    = tbl.(metric)(mask2);                          % metric values for stim 2
    cIdx  = animalIdx(mask1);                             % animal colour index (aligned with mask1)

    % Scatter with animal-coded colour
    scatter(v1, v2, 7, cmap(cIdx,:), 'filled', 'MarkerFaceAlpha', 0.3);
    hold on;
    axis equal;

    % Identity line (y = x)
    lims = [min(tbl.(metric)), max(tbl.(metric))];        % data range
    plot(lims, lims, 'k--', 'LineWidth', 1.5);
    xlim(lims);  ylim(lims);

    % Axis labels — apply display-name substitutions
    xLab = stimsNeeded{1};
    yLab = stimsNeeded{2};
    for li = 1:size(labelMap, 1)
        xLab = strrep(xLab, labelMap{li,1}, labelMap{li,2});
        yLab = strrep(yLab, labelMap{li,1}, labelMap{li,2});
    end
    xlabel(xLab);   ylabel(yLab);
    colormap(fig, cmap);

    formatAxes(gca, 8, 'helvetica');
    set(fig, 'Units', 'centimeters', 'Position', [20 20 5 5]);
end


function formatAxes(ax, fontSize, fontName)
% formatAxes  Apply consistent font styling to an axes object.
    ax.YAxis.FontSize = fontSize;   ax.YAxis.FontName = fontName;
    ax.XAxis.FontSize = fontSize;   ax.XAxis.FontName = fontName;
end


function pAdj = bhFDR(pVals)
% bhFDR  Benjamini-Hochberg FDR correction.
%   Adjusts a vector of p-values to control the false discovery rate.

    n          = numel(pVals);                  % number of tests
    [pSorted, sortIdx] = sort(pVals(:));        % sort ascending
    ranks      = (1:n)';                        % integer ranks

    % BH adjustment: p_adj(k) = min( p(k)*n/k , 1 ), enforced monotone
    pAdj       = pSorted .* n ./ ranks;         % raw BH adjustment
    pAdj       = min(pAdj, 1);                  % cap at 1
    for k = n-1:-1:1                            % enforce monotonicity from bottom up
        pAdj(k) = min(pAdj(k), pAdj(k+1));
    end

    % Unsort back to original order
    pAdj(sortIdx) = pAdj;
end