function [tblOut, figs] = AllExpDirectionTuning(expList, params)
% AllExpDirectionTuning  Pool OSI and DSI across multiple Neuropixels
%   recordings, run hierarchical bootstrap statistics for pairwise
%   stimulus comparisons, and generate publication-ready swarm plots.
%
%   This function follows the same design pattern as AllExpAnalysis:
%     1. Loop over experiments, load each stimulus object
%     2. Call DirectionTuning to compute per-neuron OSI/DSI
%     3. Pool results into a long table (one row per neuron × stimulus × metric)
%     4. Run hierBoot (Saravanan et al. 2020) for pairwise comparisons
%     5. Plot with plotSwarmBootstrapWithComparisons
%
%   INPUTS
%     expList  (1,:) double   Row vector of experiment IDs from master Excel.
%     params   Name-value     See arguments block below.
%
%   OUTPUTS
%     tblOut   table   Long table with columns:
%                        animal, insertion, stimulus, NeurID, metric, value
%     figs     struct  Figure handles (.OSI, .DSI, .prefDir)
%
%   EXAMPLES
%     % Compare OSI between moving gratings and moving ball:
%     [tbl, figs] = AllExpDirectionTuning([49:54 64:66], ...
%         stimuli = {'SDGm', 'MB'}, plot = true, PaperFig = true);
%
%     % Single stimulus — just pool and plot, no comparison:
%     [tbl, figs] = AllExpDirectionTuning([49:54], ...
%         stimuli = {'SDGm'}, plot = true);
%
%   See also: DirectionTuning, hierBoot, plotSwarmBootstrapWithComparisons,
%             AllExpAnalysis

% =========================================================================
%                         ARGUMENTS BLOCK
% =========================================================================
arguments
    expList (1,:) double                                    % row vector of experiment IDs

    % --- Stimuli to compare ---
    params.stimuli          cell    = {'SDGm'}              % cell of stimulus abbreviations to pool
                                                            % (e.g. {'SDGm','MB'}). Each must have
                                                            % a 'direction' category in its conditions.

    % --- Responsiveness filter ---
    params.threshold        double  = 0.05                  % p-value cutoff for responsiveness mask

    % --- Tuning curve construction ---
    params.aggMethod        string  = "mean"                % "mean" or "max" passed to DirectionTuning

    % --- ResponseWindow / Statistics recomputation ---
    params.overwriteRW      logical = false                 % force recompute ResponseWindow
    params.overwriteStats   logical = false                 % force recompute StatisticsPerNeuron

    % --- Hierarchical bootstrap ---
    params.nBoot            double  = 10000                 % number of bootstrap resamples
    params.rngSeed          double  = 42                    % fixed RNG seed for reproducibility

    % --- Plotting ---
    params.plot             logical = true                   % generate swarm plots
    params.PaperFig         logical = false                  % save publication-quality figures via printFig
    params.Alpha            double  = 0.4                    % dot transparency in swarm plot
    params.yMaxOSI          double  = 1                      % y-axis upper limit for OSI plot
    params.yMaxDSI          double  = 1                      % y-axis upper limit for DSI plot

    % --- Saving ---
    params.overwrite        logical = false                  % overwrite previously saved pooled results
end

% =========================================================================
%  0.  INITIALISE
% =========================================================================

nStim = numel(params.stimuli);                              % number of stimuli to pool
figs  = struct();                                           % output figure handles

% Fix RNG seed for reproducible bootstrap confidence intervals across runs
rng(params.rngSeed, 'twister');                             % Mersenne Twister, fixed seed

% -------------------------------------------------------------------------
% Build save path using the first experiment (same convention as AllExpAnalysis)
% -------------------------------------------------------------------------
NP_first    = loadNPclassFromTable(expList(1));             % load first experiment for path extraction
firstObj    = createStimulusObject(NP_first, params.stimuli{1}, 0);

% Extract the shared analysis directory from the first experiment
pathStr = extractBefore(firstObj.getAnalysisFileName, 'lizards');
saveDir = fullfile([pathStr 'lizards'], 'Combined_lizard_analysis');

% Create directory if it does not exist
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% Construct save file name incorporating experiment range and stimuli
stimLabel  = strjoin(params.stimuli, '-');                   % e.g. 'SDGm-MB'
nameOfFile = sprintf('Ex_%d-%d_DirectionTuning_%s.mat', ...
    expList(1), expList(end), stimLabel);
savePath   = fullfile(saveDir, nameOfFile);

% -------------------------------------------------------------------------
% Check for existing results (load if available and not overwriting)
% -------------------------------------------------------------------------
if exist(savePath, 'file') && ~params.overwrite
    fprintf('Loading existing results: %s\n', savePath);
    S = load(savePath, 'tblOut');                            % load the saved long table
    tblOut = S.tblOut;
else
    % =====================================================================
    %  1.  BUILD LONG TABLE — LOOP OVER EXPERIMENTS AND STIMULI
    % =====================================================================

    % Pre-allocate the long table with typed empty columns.
    % Columns: animal, insertion, stimulus, NeurID, metric, value
    tblOut = table( ...
        categorical.empty(0,1), ...                         % animal   — e.g. PV123
        categorical.empty(0,1), ...                         % insertion — stable numeric ID
        categorical.empty(0,1), ...                         % stimulus — e.g. SDGm, MB
        categorical.empty(0,1), ...                         % NeurID   — unique neuron identifier
        categorical.empty(0,1), ...                         % metric   — 'OSI' or 'DSI'
        double.empty(0,1), ...                              % value    — the index value
        'VariableNames', {'animal','insertion','stimulus','NeurID','metric','value'});

    % Maps for stable numeric IDs across experiments (same as AllExpAnalysis)
    animalMap       = containers.Map();                      % animalKey → unique integer
    insertionMap    = containers.Map();                      % insKey → unique integer
    nextAnimalIdx   = 0;                                    % counter for new animals
    nextInsertionIdx = 0;                                   % counter for new insertions

    for ei = 1:numel(expList)
        ex = expList(ei);                                   % current experiment ID
        fprintf('\n======= Experiment %d (%d/%d) =======\n', ex, ei, numel(expList));

        % -----------------------------------------------------------------
        %  1a. Load NP class for this experiment
        % -----------------------------------------------------------------
        try
            NP = loadNPclassFromTable(ex);                  % load Neuropixels data object
        catch ME
            warning('Could not load experiment %d: %s', ex, ME.message);
            continue                                        % skip to next experiment
        end

        % Extract animal ID from recording name (supports PV### and SA### patterns)
        animalID = string(regexp(NP.recordingName, 'PV\d+', 'match', 'once'));
        if animalID == ""
            animalID = string(regexp(NP.recordingName, 'SA\d+', 'match', 'once'));
        end
        if animalID == ""
            warning('Cannot extract animal ID from "%s". Skipping.', NP.recordingName);
            continue
        end

        % Extract insertion number from the directory path
        insStr = regexp(NP.recordingDir, 'Insertion\d+', 'match', 'once');
        insNum = str2double(regexp(insStr, '\d+', 'match'));
        if isempty(insNum) || isnan(insNum)
            insNum = ei;                                    % fallback: use loop index
        end

        % Register animal and insertion in maps for stable unique IDs
        animalKey = char(animalID);
        if ~animalMap.isKey(animalKey)
            nextAnimalIdx = nextAnimalIdx + 1;              % assign new index
            animalMap(animalKey) = nextAnimalIdx;
        end
        insKey = sprintf('%s__Ins%d', animalKey, insNum);
        if ~insertionMap.isKey(insKey)
            nextInsertionIdx = nextInsertionIdx + 1;        % assign new index
            insertionMap(insKey) = nextInsertionIdx;
        end
        insertionIdx = insertionMap(insKey);                % stable numeric ID for this insertion

        % -----------------------------------------------------------------
        %  1b. Get phy IDs for good somatic units (shared across stimuli)
        % -----------------------------------------------------------------
        firstStimObj = createStimulusObject(NP, params.stimuli{1}, 0);
        if isempty(firstStimObj)
            warning('Cannot create stimulus object for %s. Skipping experiment %d.', ...
                params.stimuli{1}, ex);
            continue
        end
        p_sort  = NP.convertPhySorting2tIc(firstStimObj.spikeSortingFolder, 0, 1, 1);
        phy_IDg = p_sort.phy_ID(string(p_sort.label') == 'good');

        % -----------------------------------------------------------------
        %  1c. Loop over stimulus types
        % -----------------------------------------------------------------
        for s = 1:nStim
            stimName = params.stimuli{s};                   % e.g. 'SDGm', 'MB'
            fprintf('  — Stimulus: %s\n', stimName);

            % Create the analysis object for this stimulus
            vsObj = createStimulusObject(NP, stimName, 0);
            if isempty(vsObj) || isempty(vsObj.VST)
                fprintf('    Stimulus %s not found. Skipping.\n', stimName);
                continue                                    % stimulus not present in this session
            end

            % Run DirectionTuning for this stimulus and experiment
            try
                res = DirectionTuning(vsObj, ...
                    'aggMethod',      params.aggMethod, ...
                    'threshold',      params.threshold, ...
                    'overwriteRW',    params.overwriteRW, ...
                    'overwriteStats', params.overwriteStats, ...
                    'save',           true, ...
                    'overwrite',      params.overwrite);
            catch ME
                warning('DirectionTuning failed for %s in exp %d: %s', ...
                    stimName, ex, ME.message);
                continue
            end

            % Apply responsiveness filter: only include neurons with p < threshold
            respIdx = find(res.respMask);                   % indices of responsive neurons
            nResp   = numel(respIdx);                       % number of responsive neurons
            fprintf('    Responsive: %d / %d\n', nResp, numel(res.respMask));

            if nResp == 0
                continue                                    % no responsive neurons — skip
            end

            % Build unique NeurID strings for NeurID-based pairing across stimuli.
            % Format: "ins<X>_phy<Y>" where X = stable insertion ID, Y = phy cluster ID.
            % This ensures the SAME neuron is matched across stimulus types.
            neurIDs = arrayfun(@(pid) sprintf('ins%d_phy%d', insertionIdx, pid), ...
                res.phyID(respIdx), 'UniformOutput', false);

            % Append OSI rows to the long table
            nRows = nResp;                                  % one row per responsive neuron per metric
            osiRows = table( ...
                repmat(categorical(animalID),  nRows, 1), ...   % animal
                repmat(categorical(insertionIdx), nRows, 1), ...% insertion (numeric ID)
                repmat(categorical(string(stimName)), nRows, 1), ...% stimulus
                categorical(neurIDs(:)), ...                    % NeurID
                repmat(categorical("OSI"), nRows, 1), ...       % metric
                res.OSI(respIdx), ...                           % value
                'VariableNames', tblOut.Properties.VariableNames);

            % Append DSI rows to the long table
            dsiRows = table( ...
                repmat(categorical(animalID),  nRows, 1), ...
                repmat(categorical(insertionIdx), nRows, 1), ...
                repmat(categorical(string(stimName)), nRows, 1), ...
                categorical(neurIDs(:)), ...
                repmat(categorical("DSI"), nRows, 1), ...
                res.DSI(respIdx), ...
                'VariableNames', tblOut.Properties.VariableNames);

            % Concatenate to master table
            tblOut = [tblOut; osiRows; dsiRows];            %#ok<AGROW>

        end   % stimulus loop
    end   % experiment loop

    % Remove unused categories from categorical columns
    tblOut.animal    = removecats(tblOut.animal);
    tblOut.insertion = removecats(tblOut.insertion);
    tblOut.stimulus  = removecats(tblOut.stimulus);
    tblOut.NeurID    = removecats(tblOut.NeurID);
    tblOut.metric    = removecats(tblOut.metric);

    % Save the pooled long table
    save(savePath, 'tblOut', 'params');
    fprintf('\nSaved pooled results: %s\n', savePath);
end

% =========================================================================
%  2.  GENERATE COMPARISON PAIRS (IF MORE THAN ONE STIMULUS)
% =========================================================================

% Build all pairwise combinations of stimuli for statistical testing
if nStim >= 2
    pairs = nchoosek(params.stimuli, 2);                    % [nPairs × 2] cell array
else
    pairs = {};                                             % nothing to compare
end

% =========================================================================
%  3.  PLOT OSI AND DSI WITH HIERARCHICAL BOOTSTRAP
% =========================================================================

if params.plot

    % -----------------------------------------------------------------
    %  3a. OSI plot
    % -----------------------------------------------------------------
    tblOSI = tblOut(tblOut.metric == 'OSI', :);             % filter to OSI rows only
    if ~isempty(tblOSI) && height(tblOSI) > 0

        % Compute hierBoot p-values for each stimulus pair
        psOSI = computeHierBootPvals(tblOSI, pairs, params.nBoot);

        % Generate swarm plot with bootstrap mean/CI and significance brackets
        [figOSI, ~] = plotSwarmBootstrapWithComparisons(tblOSI, pairs, psOSI, {'value'}, ...
            yLegend     = 'OSI', ...
            yMaxVis     = params.yMaxOSI, ...
            diff        = (nStim >= 2), ...                 % paired diff plot only if comparing
            Alpha       = params.Alpha, ...
            plotMeanSem = true);

        % Add title and format
        title('Orientation Selectivity Index', 'FontSize', 10);
        applyPaperAxes(gca);                                % standardise axis formatting

        % Save figure if PaperFig mode is enabled
        if params.PaperFig
            firstObj.printFig(figOSI, sprintf('OSI-comparison-%s', stimLabel), ...
                PaperFig = params.PaperFig);
        end

        figs.OSI = figOSI;                                  % store figure handle
    end

    % -----------------------------------------------------------------
    %  3b. DSI plot
    % -----------------------------------------------------------------
    tblDSI = tblOut(tblOut.metric == 'DSI', :);             % filter to DSI rows only
    if ~isempty(tblDSI) && height(tblDSI) > 0

        % Compute hierBoot p-values for each stimulus pair
        psDSI = computeHierBootPvals(tblDSI, pairs, params.nBoot);

        % Generate swarm plot with bootstrap mean/CI and significance brackets
        [figDSI, ~] = plotSwarmBootstrapWithComparisons(tblDSI, pairs, psDSI, {'value'}, ...
            yLegend     = 'DSI', ...
            yMaxVis     = params.yMaxDSI, ...
            diff        = (nStim >= 2), ...
            Alpha       = params.Alpha, ...
            plotMeanSem = true);

        title('Direction Selectivity Index', 'FontSize', 10);
        applyPaperAxes(gca);

        if params.PaperFig
            firstObj.printFig(figDSI, sprintf('DSI-comparison-%s', stimLabel), ...
                PaperFig = params.PaperFig);
        end

        figs.DSI = figDSI;
    end
end

fprintf('\nDone. Table: %d rows, %d experiments, %d stimuli.\n', ...
    height(tblOut), numel(expList), nStim);

end   % ===== END OF MAIN FUNCTION =====


% =========================================================================
%                      LOCAL HELPER FUNCTIONS
% =========================================================================


function ps = computeHierBootPvals(tbl, pairs, nBoot)
% computeHierBootPvals  Compute p-values for stimulus pair comparisons
%   using hierarchical bootstrap on paired differences.
%
%   Pairing is NeurID-based: only neurons present in BOTH stimuli of a pair
%   contribute. The hierarchy is: neurons nested within insertions nested
%   within animals, matching the mixed-model random-effects structure.
%
%   INPUTS
%     tbl    — long table with columns: animal, insertion, stimulus, NeurID, value
%     pairs  — [nPairs × 2] cell of stimulus name pairs
%     nBoot  — number of bootstrap resamples
%
%   OUTPUTS
%     ps     — [nPairs × 1] two-sided p-values

    nPairs = size(pairs, 1);                                % number of comparisons
    ps     = nan(nPairs, 1);                                % pre-allocate

    for i = 1:nPairs
        stim1 = pairs{i, 1};                                % first stimulus in pair
        stim2 = pairs{i, 2};                                % second stimulus in pair

        % Accumulate paired differences, insertion IDs, and animal IDs
        diffs   = [];                                       % paired differences (stim1 − stim2)
        insers  = [];                                       % insertion grouping variable
        animals = [];                                       % animal grouping variable

        % Loop over unique insertions to enforce NeurID-based pairing
        for ins = unique(tbl.insertion)'
            % Find rows matching this insertion and each stimulus
            idx1 = tbl.insertion == ins & tbl.stimulus == stim1;
            idx2 = tbl.insertion == ins & tbl.stimulus == stim2;

            % Extract NeurIDs present in both stimuli for this insertion
            nIDs1 = tbl.NeurID(idx1);                       % NeurIDs for stim1
            nIDs2 = tbl.NeurID(idx2);                       % NeurIDs for stim2
            shared = intersect(nIDs1, nIDs2);               % neurons in both

            if isempty(shared)
                continue                                    % no paired data for this insertion
            end

            % Extract paired values, matched by NeurID (NOT by row order)
            [~, loc1] = ismember(shared, nIDs1);            % indices into stim1 rows
            [~, loc2] = ismember(shared, nIDs2);            % indices into stim2 rows
            vals1_all = tbl.value(idx1);                    % all values for stim1 at this insertion
            vals2_all = tbl.value(idx2);                    % all values for stim2 at this insertion
            V1 = vals1_all(loc1);                           % paired values stim1
            V2 = vals2_all(loc2);                           % paired values stim2

            % Extract the animal ID for this insertion
            animal = unique(tbl.animal(idx1));               % should be a single animal

            % Append paired differences and grouping variables
            nShared = numel(shared);                        % number of paired neurons
            diffs   = [diffs;   V1 - V2];                  %#ok<AGROW>
            insers  = [insers;  double(repmat(ins, nShared, 1))];     %#ok<AGROW>
            animals = [animals; double(repmat(animal, nShared, 1))];  %#ok<AGROW>
        end

        if isempty(diffs)
            % No paired neurons found — cannot test
            ps(i) = NaN;
            fprintf('  Pair %s vs %s: no paired neurons. p = NaN.\n', stim1, stim2);
        else
            % Run hierarchical bootstrap on paired differences
            % hierBoot resamples neurons within insertions within animals
            bootDiff = hierBoot(diffs, nBoot, insers, animals);

            % Two-sided p-value: proportion of bootstrap means ≤ 0
            % (tests H0: no difference between stimuli)
            ps(i) = min(mean(bootDiff <= 0), mean(bootDiff >= 0)) * 2;

            fprintf('  Pair %s vs %s: n=%d neurons, p = %.4f\n', ...
                stim1, stim2, numel(diffs), ps(i));
        end
    end
end


function vsObj = createStimulusObject(NP, stimName, session)
% createStimulusObject  Create an analysis object for a given stimulus.
%
%   Maps stimulus abbreviations to their analysis class constructors.
%   Returns [] if the constructor throws an error.
%
%   INPUTS
%     NP        – Neuropixels data object (from loadNPclassFromTable)
%     stimName  – stimulus abbreviation (e.g. 'MB', 'SDGm', 'RG')
%     session   – session number (0 = default, 1+ = specific session)

    vsObj = [];                                             % default: empty on failure
    try
        key = getObjKey(stimName);                          % map abbreviation to shared key
        if session == 0
            % Default session — no Session argument passed to constructor
            switch key
                case 'MB',  vsObj = linearlyMovingBallAnalysis(NP);
                case 'RG',  vsObj = rectGridAnalysis(NP);
                case 'MBR', vsObj = linearlyMovingBarAnalysis(NP);
                case 'SDG', vsObj = StaticDriftingGratingAnalysis(NP);
                case 'NI',  vsObj = imageAnalysis(NP);
                case 'NV',  vsObj = movieAnalysis(NP);
                case 'FFF', vsObj = fullFieldFlashAnalysis(NP);
                otherwise
                    error('Unknown stimName: %s (key: %s)', stimName, key);
            end
        else
            % Explicit session number passed to constructor
            switch key
                case 'MB',  vsObj = linearlyMovingBallAnalysis(NP, 'Session', session);
                case 'RG',  vsObj = rectGridAnalysis(NP, 'Session', session);
                case 'MBR', vsObj = linearlyMovingBarAnalysis(NP, 'Session', session);
                case 'SDG', vsObj = StaticDriftingGratingAnalysis(NP, 'Session', session);
                case 'NI',  vsObj = imageAnalysis(NP, 'Session', session);
                case 'NV',  vsObj = movieAnalysis(NP, 'Session', session);
                case 'FFF', vsObj = fullFieldFlashAnalysis(NP, 'Session', session);
                otherwise
                    error('Unknown stimName: %s (key: %s)', stimName, key);
            end
        end
    catch ME
        fprintf('  Could not create %s Session=%d: %s\n', stimName, session, ME.message);
        vsObj = [];                                         % return empty on failure
    end
end


function key = getObjKey(stimName)
% getObjKey  Map stimulus abbreviation to shared analysis-object key.
%   SDGm and SDGs share a single StaticDriftingGratingAnalysis object.
    switch stimName
        case {'SDGm', 'SDGs'},  key = 'SDG';               % both use the same constructor
        otherwise,              key = stimName;             % all others map to themselves
    end
end