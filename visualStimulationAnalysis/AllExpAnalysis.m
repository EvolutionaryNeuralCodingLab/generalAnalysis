function [tempTable] = AllExpAnalysis(expList, Stims2Comp, params)
% PlotZScoreComparison - Compare z-scores and spike rates across visual stimuli
%   across multiple Neuropixels recordings.
%
%   Loads pre-computed statistical results (z-scores, p-values, spike rates)
%   for each experiment in expList, filters neurons by responsiveness, pools
%   data across recordings, runs hierarchical bootstrapping for group-level
%   inference, and generates swarm + scatter plots for publication.
%
%   INPUTS:
%     expList    - (1,:) double  Row vector of experiment indices from the
%                               Excel master list.
%     Stims2Comp - cell          Cell array of stimulus abbreviations defining
%                               the comparison order. The FIRST element is the
%                               "anchor" stimulus used to select responsive
%                               neurons (unless EachStimSignif=true).
%                               E.g. {'MB','RG','MBR'}.
%     params     - name-value    Optional parameters (see arguments block).
%
%   OUTPUT:
%     fig        - figure handle of the last figure created.
%
% -------------------------------------------------------------------------
% KNOWN BUGS / ISSUES (see inline BUG comments for exact locations):
%   BUG-1  [CRASH]  splitapply fails on empty TableStimComp when no units
%                   pass significance threshold. → Guard added below.
%   BUG-2  [LOGIC]  fprintf prints recording name BEFORE NP is loaded for
%                   the current experiment, so iteration 1 always prints the
%                   name from the pre-loop load (expList(1)).
%   BUG-3  [LOGIC]  Insertion counter: AnimalI is updated inside the first
%                   `if Animal~=AnimalI` block, so the second block
%                   (which also checks Animal~=AnimalI) always sees them as
%                   equal, and a new animal's first insertion is never counted
%                   as new unless the insertion number also differs.
%   BUG-4  [LOGIC]  When SDG is absent, `sumNeurSDG=0` is set (new var) but
%                   `sumNeurSDGm` and `sumNeurSDGs` keep their last stale
%                   values, so sumNeurSDGmt{j} / sumNeurSDGst{j} are wrong.
%   BUG-5  [DEBUG]  `2+2` is a leftover breakpoint stub — does nothing but
%                   is confusing in published code.
%   BUG-6  [STRUCT] S.groupStatsP_ZscoreCompare should be
%                   S.groupStats.P_ZscoreCompare (inconsistent nesting vs
%                   the spike-rate equivalent).
%   BUG-7  [PREALLOC] totalU, pvalsRG, pvalsMB, pvalsNI, pvalsNV etc. are
%                   not pre-allocated before the for-loop (unlike zScoresMB
%                   etc.), causing dynamic growth inside the loop.
%
% SUGGESTIONS:
%   SUGG-1  Refactor the 7-stimulus × 3-method conditional blocks into a
%           helper function (e.g. runStimAnalysis(vs, method, params)) to
%           drastically reduce code length and risk of copy-paste bugs.
%   SUGG-2  Replace the -inf sentinel for absent stimuli with NaN.  NaN
%           propagates safely through most MATLAB statistics functions;
%           -inf does not, and requires scattered special-case filtering.
%   SUGG-3  For a publication, consider applying FDR correction
%           (Benjamini-Hochberg) across neurons before applying the
%           significance threshold, rather than using raw p < threshold.
%   SUGG-4  For scatter plots, if spike rates span >1 order of magnitude,
%           log-scaled axes improve readability (set(gca,'XScale','log',...)).
%   SUGG-5  randiColors (subsampling index from plotSwarmBootstrapWithComparisons)
%           is reused in scatter plots.  If the swarm function subsamples
%           non-uniformly, the scatter could misrepresent the distribution.
%           Either plot all points or make subsampling explicit and documented.
%   SUGG-6  The `eval(zscoresC1{1})` pattern is fragile.  Prefer a struct
%           or containers.Map to look up variables by name.

% -------------------------------------------------------------------------
arguments
    expList    (1,:) double  % Row vector of experiment IDs from master Excel table
    Stims2Comp cell          % Cell array: comparison order, e.g. {'MB','RG','MBR'}.
                             %   First element selects the anchor stimulus for
                             %   filtering responsive neurons.
    params.threshold         = 0.05   % p-value significance threshold for responsiveness
    params.diffResp          = false  % If true, use spike-rate difference (resp-baseline)
                                      %   instead of absolute response rate
    params.overwrite         = false  % If true, recompute and overwrite saved combined file
    params.StimsPresent      = {'MB','RG'}  % Stimuli present in ALL recordings (minimum set)
    params.StimsNotPresent   = {}           % Stimuli known to be absent (currently unused)
    params.StimsToCompare    = {}    % Two-element cell: which stimuli to use in the scatter
                                     %   sub-panel (default: 1st and 2nd of Stims2Comp)
    params.overwriteResponse = false  % Force re-run of ResponseWindow analysis
    params.overwriteStats    = false  % Force re-run of per-neuron statistics
    params.overwriteGroupStats = false % Force re-run of group-level bootstrapping
    params.RespDurationWin   = 100    % Duration (ms) of the response window (passed down)
    params.shuffles          = 2000   % Number of shuffles / bootstrap iterations for
                                      %   per-neuron statistics
    params.StatMethod        = 'ObsWindow'  % Statistical method:
                                            %   'ObsWindow'         – shuffling analysis
                                            %   'bootsrapRespBase'  – per-neuron bootstrap
                                            %   'maxPermuteTest'    – permutation test
    params.ignoreNonSignif   = false  % When true, zero out z-scores for neurons that are
                                      %   not significant for the non-anchor stimuli
    params.EachStimSignif    = false  % If true, use each stimulus's own responsive neurons
                                      %   (default: use anchor stimulus's responsive neurons)
    params.ComparePairs      = {}     % Cell of stimulus pairs for pairwise comparison.
                                      %   Recommended over the multi-stimulus mode.
                                      %   E.g. {'MB','RG'} or {'MB','RG';'MB','MBR'}
    params.PaperFig logical  = false  % If true, save figures via vs.printFig
    params.useZmean logical  = true   % Instead of the spikerate from pvals response, use the max response-baseline - null distribution
end

% =========================================================================
% SECTION 1 – INITIALISE BOOKKEEPING VARIABLES
% =========================================================================

% Running counters for unique animals and probe insertions encountered
animal    = 0;
insertion = 0;

% Pre-allocate per-experiment cell arrays (one cell per experiment in expList)
n = numel(expList);  % total number of experiments to process

% Animal/insertion labels for each neuron (repeated per neuron count)
animalVector    = cell(1, n);
insertionVector = cell(1, n);

% Z-scores filtered to neurons responsive to the anchor stimulus
zScoresMB   = cell(1, n);
zScoresRG   = cell(1, n);
zScoresMBR  = cell(1, n);
zScoresFFF  = cell(1, n);
zScoresSDGm = cell(1, n);  % drifting gratings – moving condition
zScoresNI   = cell(1, n);

% Spike rates (peak across directions/speeds) for anchor-responsive neurons
spKrMB   = cell(1, n);
spKrRG   = cell(1, n);
spKrMBR  = cell(1, n);
spKrFFF  = cell(1, n);
spKrSDGm = cell(1, n);

% Spike-rate difference (response – baseline) for anchor-responsive neurons
diffSpkMB   = cell(1, n);
diffSpkRG   = cell(1, n);
diffSpkMBR  = cell(1, n);
diffSpkFFF  = cell(1, n);
diffSpkSDGm = cell(1, n);

% Natural image / video variables (declared but not pre-sized above)
spKrNI   = cell(1, n);
spKrNV   = cell(1, n);
diffSpkNI = cell(1, n);
diffSpkNV = cell(1, n);

% BUG-7: The following accumulator cell arrays are NOT pre-allocated here.
%        They grow dynamically inside the loop.  Add pre-allocation if
%        performance matters (e.g. pvalsRG = cell(1,n); etc.).

% Tracker strings for detecting animal/insertion changes between experiments
j         = 1;      % experiment counter (1-based index into cell arrays)
AnimalI   = "";     % animal ID seen in the previous iteration
InsertionI = 0;     % insertion number seen in the previous iteration

% =========================================================================
% SECTION 2 – DETERMINE OUTPUT FILE PATH AND WHETHER THE LOOP IS NEEDED
% =========================================================================

% Load the first experiment to extract file-path information and response window
NP = loadNPclassFromTable(expList(1));   % load Neuropixels recording object
vs = linearlyMovingBallAnalysis(NP);     % run moving-ball analysis (for path info)

% Read response window used in moving-ball analysis (assumed identical across
% experiments — this assumption is NOT verified across experiments)
MBvs = vs.ResponseWindow;  % cache the response-window struct

% Build the filename for the pooled/combined output .mat file
nameOfFile = sprintf('\\Ex_%d-%d_Combined_Neural_responses_%s_filtered.mat', ...
    expList(1), expList(end), Stims2Comp{1});

% Extract base path up to (and including) the 'lizards' folder
p = extractBefore(vs.getAnalysisFileName, 'lizards');
p = [p 'lizards'];

% Create the 'Combined_lizard_analysis' subdirectory if it does not exist
if ~exist([p '\Combined_lizard_analysis'], 'dir')
    cd(p)
    mkdir Combined_lizard_analysis
end
saveDir = [p '\Combined_lizard_analysis'];  % full path to output folder

% Decide whether to run the per-experiment for-loop:
%   • Skip if a saved file exists with the same experiment list AND overwrite=false
%   • Otherwise run the loop to build and save pooled data
if exist([saveDir nameOfFile], 'file') == 2 && ~params.overwrite
    S = load([saveDir nameOfFile]);   % load previously saved pooled data
    expList2 = S.expList;             % experiment list stored inside the file

    if isequal(expList2, expList)
        forloop = false;   % saved data matches → skip re-processing
    else
        forloop = true;    % experiment list changed → must re-process
    end
else
    forloop = true;        % file does not exist or overwrite requested
end

% =========================================================================
% SECTION 3 – INITIALISE LONG-FORMAT TABLES
% =========================================================================

% longTablePairComp: one row per neuron × stimulus for the pairwise comparison.
%   Columns: animal ID, insertion ID, stimulus name, neuron ID,
%            z-score, and spike rate.
longTablePairComp = table( ...
    categorical.empty(0,1), ...   % animal
    categorical.empty(0,1), ...   % insertion
    categorical.empty(0,1), ...   % stimulus
    categorical.empty(0,1), ...   % NeurID
    double.empty(0,1), ...        % Z-score
    double.empty(0,1), ...        % SpkR
    'VariableNames', {'animal','insertion','stimulus','NeurID','Z-score','SpkR'});

% longTable: one row per insertion × stimulus; stores counts of responsive
%   and total somatic neurons for fraction-responsive analysis.
longTable = table( ...
    categorical.empty(0,1), ...   % animal
    categorical.empty(0,1), ...   % insertion
    categorical.empty(0,1), ...   % stimulus
    double.empty(0,1), ...        % respNeur  – number of responsive neurons
    double.empty(0,1), ...        % totalSomaticN – total neurons in recording
    'VariableNames', {'animal','insertion','stimulus','respNeur','totalSomaticN'});

% =========================================================================
% SECTION 4 – PER-EXPERIMENT FOR-LOOP
% =========================================================================

if forloop
    for ex = expList   % iterate over each experiment ID

        % BUG-2: fprintf is called BEFORE NP is loaded for the current
        %        experiment.  On the first iteration this prints the name
        %        from expList(1) (loaded before the loop), not from `ex`.
        %        FIX: move this fprintf to AFTER the loadNPclassFromTable call.
        fprintf('Processing recording: %s .\n', NP.recordingName)

        % Load the Neuropixels recording object for this experiment
        NP  = loadNPclassFromTable(ex);

        % Instantiate analysis objects for the two stimuli present in all sessions
        vs  = linearlyMovingBallAnalysis(NP);   % moving ball (MB)
        vsR = rectGridAnalysis(NP);             % rectangular grid (RG)

        % Extract animal ID using regex (expects pattern 'PV##' in filename)
        Animal = string(regexp(vs.getAnalysisFileName, 'PV\d+', 'match', 'once'));

        % Add placeholder rows to longTable for MB and RG (always present)
        longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("MB"), 0, 0};
        longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("RG"), 0, 0};

        % ------------------------------------------------------------------
        % 4a – Try to load optional stimuli; fall back to a dummy analysis
        %      object (vsR / vs) when the stimulus was not shown, to keep
        %      all downstream variable names defined.
        % ------------------------------------------------------------------

        % Moving Bar (MBR)
        try
            vsBr = linearlyMovingBarAnalysis(NP);
            params.StimsPresent{3} = 'MBR';
            if isempty(vsBr.VST)
                error('Moving Bar stimulus not found.\n')
            else
                longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("MBR"), 0, 0};
            end
        catch
            params.StimsPresent{3} = '';        % mark as absent
            fprintf('Moving Bar stimulus not found.\n')
            vsBr = linearlyMovingBallAnalysis(NP);  % dummy placeholder (same class)
        end

        % Static / Drifting Gratings (SDG)
        try
            vsG = StaticDriftingGratingAnalysis(NP);
            params.StimsPresent{4} = 'SDGm';
            if isempty(vsG.VST)
                error('Gratings stimulus not found.\n')
            else
                longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("SDGm"), 0, 0};
                longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("SDGs"), 0, 0};
            end
        catch
            params.StimsPresent{4} = '';
            fprintf('Gratings stimulus not found.\n')
            vsG = rectGridAnalysis(NP);   % dummy placeholder
        end

        % Natural Images (NI)
        try
            vsNI = imageAnalysis(NP);
            params.StimsPresent{5} = 'NI';
            if isempty(vsNI.VST)
                error('Natural images stimulus not found.\n')
            else
                longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("NI"), 0, 0};
            end
        catch
            params.StimsPresent{5} = '';
            fprintf('Natural images stimulus not found.\n')
            vsNI = rectGridAnalysis(NP);   % dummy placeholder
        end

        % Natural Video (NV)
        try
            vsNV = movieAnalysis(NP);
            params.StimsPresent{6} = 'NV';
            if isempty(vsNV.VST)
                error('Natural video stimulus not found.\n')
            else
                longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("NV"), 0, 0};
            end
        catch
            params.StimsPresent{6} = '';
            fprintf('Natural video stimulus not found.\n')
            vsNV = rectGridAnalysis(NP);   % dummy placeholder
        end

        % Full-Field Flash (FFF)
        try
            vsFFF = fullFieldFlashAnalysis(NP);
            params.StimsPresent{7} = 'FFF';
            if isempty(vsFFF.VST)
                error('FFF stimulus not found.\n')
            else
                longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("FFF"), 0, 0};
            end
        catch
            params.StimsPresent{7} = '';
            fprintf('FFF stimulus not found.\n')
            vsFFF = rectGridAnalysis(NP);   % dummy placeholder
        end

        % ------------------------------------------------------------------
        % 4b – Run response-window and statistical analyses for each stimulus.
        %      Only compute stats for stimuli that are (a) present AND
        %      (b) included in Stims2Comp.  For absent/excluded stimuli the
        %      analysis object already holds dummy data, so just call
        %      ResponseWindow without arguments to load any cached result.
        %
        %      SUGG-1: This block repeats ~7 times with identical structure.
        %              Wrap in a helper: runStimAnalysis(vsObj, method, params).
        % ------------------------------------------------------------------

        % Moving Ball
        if isequal(params.StimsPresent{1},'') || ~ismember(params.StimsPresent{1}, Stims2Comp)
            vs.ResponseWindow;   % load cached window only (no recompute)
        else
            vs.ResponseWindow('overwrite', params.overwriteResponse, ...
                              'durationWindow', params.RespDurationWin);
            if     isequal(params.StatMethod,'ObsWindow')
                vs.ShufflingAnalysis('overwrite', params.overwriteStats, ...
                                     "N_bootstrap", params.shuffles);
            elseif isequal(params.StatMethod,'bootsrapRespBase')
                vs.BootstrapPerNeuron('overwrite', params.overwriteStats);
            elseif isequal(params.StatMethod,'maxPermuteTest')
                vs.StatisticsPerNeuron('overwrite', params.overwriteStats);
            end
        end

        % Rect Grid
        if isequal(params.StimsPresent{2},'') || ~ismember(params.StimsPresent{2}, Stims2Comp)
            vsR.ResponseWindow;
        else
            vsR.ResponseWindow('overwrite', params.overwriteResponse, ...
                               'durationWindow', params.RespDurationWin);
            if     isequal(params.StatMethod,'ObsWindow')
                vsR.ShufflingAnalysis('overwrite', params.overwriteStats, ...
                                      "N_bootstrap", params.shuffles);
            elseif isequal(params.StatMethod,'bootsrapRespBase')
                vsR.BootstrapPerNeuron('overwrite', params.overwriteStats);
            elseif isequal(params.StatMethod,'maxPermuteTest')
                vsR.StatisticsPerNeuron('overwrite', params.overwriteStats);
            end
        end

        % Moving Bar
        if isequal(params.StimsPresent{3},'') || ~ismember(params.StimsPresent{3}, Stims2Comp)
            vsBr.ResponseWindow;
        else
            vsBr.ResponseWindow('overwrite', params.overwriteResponse, ...
                                'durationWindow', params.RespDurationWin);
            if     isequal(params.StatMethod,'ObsWindow')
                vsBr.ShufflingAnalysis('overwrite', params.overwriteStats, ...
                                       "N_bootstrap", params.shuffles);
            elseif isequal(params.StatMethod,'bootsrapRespBase')
                vsBr.BootstrapPerNeuron('overwrite', params.overwriteStats);
            elseif isequal(params.StatMethod,'maxPermuteTest')
                vsBr.StatisticsPerNeuron('overwrite', params.overwriteStats);
            end
        end

        % Gratings
        if isequal(params.StimsPresent{4},'') || ~ismember(params.StimsPresent{4}, Stims2Comp)
            vsG.ResponseWindow;
        else
            vsG.ResponseWindow('overwrite', params.overwriteResponse, ...
                               'durationWindow', params.RespDurationWin);
            if     isequal(params.StatMethod,'ObsWindow')
                vsG.ShufflingAnalysis('overwrite', params.overwriteStats, ...
                                      "N_bootstrap", params.shuffles);
            elseif isequal(params.StatMethod,'bootsrapRespBase')
                vsG.BootstrapPerNeuron('overwrite', params.overwriteStats);
            elseif isequal(params.StatMethod,'maxPermuteTest')
                vsG.StatisticsPerNeuron('overwrite', params.overwriteStats);
            end
        end

        % Natural Images
        if isequal(params.StimsPresent{5},'') || ~ismember(params.StimsPresent{5}, Stims2Comp)
            vsNI.ResponseWindow;
        else
            vsNI.ResponseWindow('overwrite', params.overwriteResponse, ...
                                'durationWindow', params.RespDurationWin);
            if     isequal(params.StatMethod,'ObsWindow')
                vsNI.ShufflingAnalysis('overwrite', params.overwriteStats, ...
                                       "N_bootstrap", params.shuffles);
            elseif isequal(params.StatMethod,'bootsrapRespBase')
                vsNI.BootstrapPerNeuron('overwrite', params.overwriteStats);
            elseif isequal(params.StatMethod,'maxPermuteTest')
                vsNI.StatisticsPerNeuron('overwrite', params.overwriteStats);
            end
        end

        % Natural Video
        if isequal(params.StimsPresent{6},'') || ~ismember(params.StimsPresent{6}, Stims2Comp)
            vsNV.ResponseWindow;
        else
            vsNV.ResponseWindow('overwrite', params.overwriteResponse, ...
                                'durationWindow', params.RespDurationWin);
            if     isequal(params.StatMethod,'ObsWindow')
                vsNV.ShufflingAnalysis('overwrite', params.overwriteStats, ...
                                       "N_bootstrap", params.shuffles);
            elseif isequal(params.StatMethod,'bootsrapRespBase')
                vsNV.BootstrapPerNeuron('overwrite', params.overwriteStats);
            elseif isequal(params.StatMethod,'maxPermuteTest')
                vsNV.StatisticsPerNeuron('overwrite', params.overwriteStats);
            end
        end

        % Full-Field Flash
        if isequal(params.StimsPresent{7},'') || ~ismember(params.StimsPresent{7}, Stims2Comp)
            vsFFF.ResponseWindow;
        else
            vsFFF.ResponseWindow('overwrite', params.overwriteResponse, ...
                                 'durationWindow', params.RespDurationWin);
            if     isequal(params.StatMethod,'ObsWindow')
                vsFFF.ShufflingAnalysis('overwrite', params.overwriteStats, ...
                                        "N_bootstrap", params.shuffles);
            elseif isequal(params.StatMethod,'bootsrapRespBase')
                vsFFF.BootstrapPerNeuron('overwrite', params.overwriteStats);
            elseif isequal(params.StatMethod,'maxPermuteTest')
                vsFFF.StatisticsPerNeuron('overwrite', params.overwriteStats);
            end
        end

        % ------------------------------------------------------------------
        % 4c – Retrieve statistics structs (dispatch on chosen method)
        % ------------------------------------------------------------------

        if isequal(params.StatMethod,'ObsWindow')
            statsMB  = vs.ShufflingAnalysis;
            statsRG  = vsR.ShufflingAnalysis;
            statsMBR = vsBr.ShufflingAnalysis;
            statsSDG = vsG.ShufflingAnalysis;
            statsFFF = vsFFF.ShufflingAnalysis;
            statsNI  = vsNI.ShufflingAnalysis;
            statsNV  = vsNV.ShufflingAnalysis;
        elseif isequal(params.StatMethod,'bootsrapRespBase')
            statsMB  = vs.BootstrapPerNeuron;
            statsRG  = vsR.BootstrapPerNeuron;
            statsMBR = vsBr.BootstrapPerNeuron;
            statsSDG = vsG.BootstrapPerNeuron;
            statsFFF = vsFFF.BootstrapPerNeuron;
            statsNI  = vsNI.BootstrapPerNeuron;
            statsNV  = vsNV.BootstrapPerNeuron;
        else   % maxPermuteTest
            statsMB  = vs.StatisticsPerNeuron;
            statsRG  = vsR.StatisticsPerNeuron;
            statsMBR = vsBr.StatisticsPerNeuron;
            statsSDG = vsG.StatisticsPerNeuron;
            statsFFF = vsFFF.StatisticsPerNeuron;
            statsNI  = vsNI.StatisticsPerNeuron;
            statsNV  = vsNV.StatisticsPerNeuron;
        end

        % Retrieve response-window structs (used for spike-rate / diff columns)
        rwRG  = vsR.ResponseWindow;
        rwMB  = vs.ResponseWindow;
        rwMBR = vsBr.ResponseWindow;
        rwFFF = vsFFF.ResponseWindow;
        rwSDG = vsG.ResponseWindow;
        rwNI  = vsNI.ResponseWindow;
        rwNV  = vsNV.ResponseWindow;

        % ------------------------------------------------------------------
        % 4d – Extract z-scores, p-values, and spike rates per stimulus
        % ------------------------------------------------------------------

        % --- Moving Ball ---
        % Use Speed1 by default; overwrite with Speed2 if it exists
        % (Speed2 is faster; the convention is to use the most salient speed)
        zScores_MB  = statsMB.Speed1.ZScoreU;
        pValuesMB   = statsMB.Speed1.pvalsResponse;
        if params.useZmean
            spkR_MB = statsMB.Speed1.z_mean; 
        else
            spkR_MB     = max(rwMB.Speed1.NeuronVals(:,:,4), [], 2);  % max across directions, col 4 = resp. rate
        end

        spkDiff_MB  = max(rwMB.Speed1.NeuronVals(:,:,5), [], 2);  % col 5 = response – baseline

        if isfield(statsMB, 'Speed2')   % if a second (faster) speed was presented
            zScores_MB  = statsMB.Speed2.ZScoreU;
            pValuesMB   = statsMB.Speed2.pvalsResponse;
            if params.useZmean
                spkR_MB = statsMB.Speed1.z_mean';
            else
                spkR_MB     = max(rwMB.Speed1.NeuronVals(:,:,4), [], 2);  % max across directions, col 4 = resp. rate
            end
            spkDiff_MB  = max(rwMB.Speed2.NeuronVals(:,:,5), [], 2);
        end

        % Store total unit count for this recording
        % BUG-7: totalU not pre-allocated; grows dynamically
        totalU{j} = numel(zScores_MB);

        % --- Rect Grid ---
        zScores_RG  = statsRG.ZScoreU;
        pValuesRG   = statsRG.pvalsResponse;
        if params.useZmean
            spkR_RG = statsRG.z_mean'; 
        else
            spkR_RG     = max(rwRG.NeuronVals(:,:,4), [], 2);  % max across directions, col 4 = resp. rate
        end
        spkDiff_RG  = max(rwRG.NeuronVals(:,:,5), [], 2);

        % --- Moving Bar ---
        zScores_MBR = statsMBR.Speed1.ZScoreU;
        pValuesMBR  = statsMBR.Speed1.pvalsResponse;
        spkR_MBR    = max(rwMBR.Speed1.NeuronVals(:,:,4), [], 2);
        spkDiff_MBR = max(rwMBR.Speed1.NeuronVals(:,:,5), [], 2);

        % --- Full-Field Flash ---
        zScores_FFF = statsFFF.ZScoreU;
        pValuesFFF  = statsFFF.pvalsResponse;
        spkR_FFF    = max(rwFFF.NeuronVals(:,:,4), [], 2);
        spkDiff_FFF = max(rwFFF.NeuronVals(:,:,5), [], 2);

        % --- Drifting / Static Gratings ---
        % When SDG is absent, statsSDG holds dummy RG data (placeholder object).
        % When present the struct has a .Moving and .Static subfield.
        if isequal(params.StimsPresent{4},'')
            % SDG not recorded: use dummy data (will be set to -inf below)
            zScores_SDGm = statsSDG.ZScoreU;
            pValuesSDGm  = statsSDG.pvalsResponse;
            spkR_SDGm    = max(rwSDG.NeuronVals(:,:,4), [], 2);
            spkDiff_SDGm = max(rwSDG.NeuronVals(:,:,5), [], 2);

            zScores_SDGs = statsSDG.ZScoreU;    % same dummy for static
            pValuesSDGs  = statsSDG.pvalsResponse;
            spkR_SDGs    = max(rwSDG.NeuronVals(:,:,4), [], 2);
            spkDiff_SDGs = max(rwSDG.NeuronVals(:,:,5), [], 2);
        else
            % SDG recorded: separate moving and static conditions
            zScores_SDGm = statsSDG.Moving.ZScoreU;
            pValuesSDGm  = statsSDG.Moving.pvalsResponse;
            spkR_SDGm    = max(rwSDG.Moving.NeuronVals(:,:,4), [], 2);
            spkDiff_SDGm = max(rwSDG.Moving.NeuronVals(:,:,5), [], 2);

            zScores_SDGs = statsSDG.Static.ZScoreU;
            pValuesSDGs  = statsSDG.Static.pvalsResponse;
            spkR_SDGs    = max(rwSDG.Static.NeuronVals(:,:,4), [], 2);
            spkDiff_SDGs = max(rwSDG.Static.NeuronVals(:,:,5), [], 2);
        end

        % --- Natural Images ---
        zScores_NI  = statsNI.ZScoreU;
        pValuesNI   = statsNI.pvalsResponse;
        spkR_NI     = max(rwNI.NeuronVals(:,:,4), [], 2);
        spkDiff_NI  = max(rwNI.NeuronVals(:,:,5), [], 2);

        % --- Natural Video ---
        zScores_NV  = statsNV.ZScoreU;
        pValuesNV   = statsNV.pvalsResponse;
        spkR_NV     = max(rwNV.NeuronVals(:,:,4), [], 2);
        spkDiff_NV  = max(rwNV.NeuronVals(:,:,5), [], 2);

        % ------------------------------------------------------------------
        % 4e – For non-ObsWindow methods, overwrite spike rates with the
        %      mean observed response stored in the stats struct
        %      (ObsWindow stores rates in rwXX; others store in stats struct)
        % ------------------------------------------------------------------

        if isequal(params.StatMethod,'bootsrapRespBase') %Take mean across all responses
            spkR_NV  = mean(statsNV.ObsResponse, 1)';
            spkR_NI  = mean(statsNI.ObsResponse, 1)';

            try
                spkR_SDGs = mean(statsSDG.Static.ObsResponse, 1)';
                spkR_SDGm = mean(statsSDG.Moving.ObsResponse, 1)';
            catch
                % Fallback: single-condition SDG struct (older data format)
                spkR_SDGs = mean(statsSDG.ObsResponse, 1)';
                spkR_SDGm = mean(statsSDG.ObsResponse, 1)';
            end

            spkR_FFF = mean(statsFFF.ObsResponse, 1)';

            try
                spkR_MBR = mean(statsMBR.Speed1.ObsResponse, 1)';
            catch
                spkR_MBR = mean(statsMBR.ObsResponse, 1)';
            end

            spkR_RG = mean(statsRG.ObsResponse, 1)';

            if isfield(statsMB, 'Speed2')
                spkR_MB = mean(statsMB.Speed2.ObsResponse)';
            else
                spkR_MB = mean(statsMB.Speed1.ObsResponse)';
            end
        end

        % ------------------------------------------------------------------
        % 4f – Optional: suppress z-scores for neurons non-significant in
        %      stimuli OTHER than the anchor by setting them to -1000
        %      (acts as a hard "must respond to everything" filter)
        % ------------------------------------------------------------------

        if params.ignoreNonSignif
            zScores_NV(pValuesNV   > params.threshold) = -1000;
            zScores_NI(pValuesNI   > params.threshold) = -1000;
            zScores_SDGs(pValuesSDGs > params.threshold) = -1000;
            zScores_SDGm(pValuesSDGm > params.threshold) = -1000;
            zScores_FFF(pValuesFFF   > params.threshold) = -1000;
            zScores_MBR(pValuesMBR   > params.threshold) = -1000;
            zScores_RG(pValuesRG     > params.threshold) = -1000;
            zScores_MB(pValuesMB     > params.threshold) = -1000;
        end

        % ------------------------------------------------------------------
        % 4g – Identify the anchor p-value vector using the first element of
        %      Stims2Comp (or the ComparePairs cell) via name matching
        % ------------------------------------------------------------------

        % Build a 2-row lookup: row 1 = variable names, row 2 = actual vectors
        pvals = {'pValuesMB','pValuesRG','pValuesMBR','pValuesFFF', ...
                 'pValuesSDGm','pValuesSDGs','pValuesNI','pValuesNV'; ...
                  pValuesMB,   pValuesRG,   pValuesMBR,  pValuesFFF, ...
                  pValuesSDGm, pValuesSDGs, pValuesNI,   pValuesNV};

        % Find column whose name ends with the anchor stimulus label
        [~, col] = find(cellfun(@(x) ischar(x) && endsWith(x, Stims2Comp{1}), pvals));
        % `row` is unused here — [~,col] is sufficient

        % ------------------------------------------------------------------
        % 4h – Build pairwise comparison table entries (ComparePairs mode)
        % ------------------------------------------------------------------

        for i = 1:numel(params.ComparePairs)
            % Find the column in pvals whose name ends with the i-th pair member
            [~, colPair] = find(cellfun(@(x) ischar(x) && endsWith(x, params.ComparePairs{i}), pvals));
            pvalsC{i} = pvals{2, colPair};  % store the actual p-value vector
        end

        % Use `who` + eval to look up z-score and spike-rate variables by name
        % SUGG-6: Replace eval with a struct lookup for robustness
        vars = who;

        % Get z-scores for the first stimulus in the pair
        zscoresC1 = vars(contains(vars, sprintf('zScores_%s', params.ComparePairs{1})));
        zscoresC1 = eval(zscoresC1{1});
        unitIDs   = 1:numel(zscoresC1);

        % Filter to neurons significant for EITHER stimulus in the pair
        sigMask   = pvalsC{1} < params.threshold | pvalsC{2} < params.threshold;
        zscoresC1 = zscoresC1(sigMask);

        spkRC1    = vars(contains(vars, sprintf('spkR_%s', params.ComparePairs{1})));
        spkRC1    = eval(spkRC1{1});
        spkRC1    = spkRC1(sigMask);
        unitIDs   = unitIDs(sigMask);   % keep only IDs for significant neurons

        % Get z-scores for the second stimulus in the pair (same mask)
        zscoresC2 = vars(contains(vars, sprintf('zScores_%s', params.ComparePairs{2})));
        zscoresC2 = eval(zscoresC2{1});
        zscoresC2 = zscoresC2(sigMask);

        spkRC2    = vars(contains(vars, sprintf('spkR_%s', params.ComparePairs{2})));
        spkRC2    = eval(spkRC2{1});
        spkRC2    = spkRC2(sigMask);

        % Append rows to longTablePairComp for this recording if any units found
        if ~isempty(unitIDs)
            try
                TableC1 = table( ...
                    categorical(cellstr(repmat(Animal, numel(unitIDs), 1))), ...
                    categorical(repmat(j, numel(unitIDs), 1)), ...
                    categorical(cellstr(repmat(params.ComparePairs{1}, numel(unitIDs), 1))), ...
                    categorical(unitIDs)', zscoresC1', spkRC1, ...
                    'VariableNames', {'animal','insertion','stimulus','NeurID','Z-score','SpkR'});

            catch
                TableC1 = table( ...
                    categorical(cellstr(repmat(Animal, numel(unitIDs), 1))), ...
                    categorical(repmat(j, numel(unitIDs), 1)), ...
                    categorical(cellstr(repmat(params.ComparePairs{1}, numel(unitIDs), 1))), ...
                    categorical(unitIDs)', zscoresC1', spkRC1', ...
                    'VariableNames', {'animal','insertion','stimulus','NeurID','Z-score','SpkR'});
            end

            TableC2 = table( ...
                categorical(cellstr(repmat(Animal, numel(unitIDs), 1))), ...
                categorical(repmat(j, numel(unitIDs), 1)), ...
                categorical(cellstr(repmat(params.ComparePairs{2}, numel(unitIDs), 1))), ...
                categorical(unitIDs)', zscoresC2', spkRC2, ...
                'VariableNames', {'animal','insertion','stimulus','NeurID','Z-score','SpkR'});

            longTablePairComp = [longTablePairComp; TableC1; TableC2];
        end

        % The anchor p-value vector (for filtering neurons in all stimuli below)
        pvalsStimSelected = pvals{2, col};

        % ------------------------------------------------------------------
        % 4i – Filter each stimulus's data to anchor-responsive neurons
        %      and compute "general" (self-responsive) neuron counts
        % ------------------------------------------------------------------
        % Convention:  suffix 's' = filtered to anchor-responsive neurons
        %              suffix 'g' = filtered to self-responsive neurons
        % respIndexes accumulates union of responsive neuron indices across stims

        respIndexes = [];   % will hold all neuron indices responsive to any stim

        % ---- Moving Ball ----
        % Anchor-responsive subset
        zScores_MBs  = zScores_MB( pvalsStimSelected <= params.threshold);
        spkR_MBs     = spkR_MB(    pvalsStimSelected <= params.threshold);
        spkDiff_MBs  = spkDiff_MB( pvalsStimSelected <= params.threshold);
        pvals_MB     = pValuesMB(  pvalsStimSelected <= params.threshold);

        % Self-responsive subset (significant for MB regardless of anchor)
        zScores_MBg  = zScores_MB( pValuesMB <= params.threshold);
        sumNeurMB    = numel(zScores_MBg);   % count of MB-responsive neurons
        spkR_MBg     = spkR_MB(    pValuesMB <= params.threshold);
        spkDiff_MBg  = spkDiff_MB( pValuesMB <= params.threshold);
        respIndexes  = [respIndexes, find(pValuesMB <= params.threshold)];

        % Update longTable with responsive / total counts for this insertion × MB
        try
            idx = (longTable.insertion == categorical(j)) & ...
                  (longTable.stimulus  == categorical("MB"));
            longTable.respNeur(idx)      = sumNeurMB;
            longTable.totalSomaticN(idx) = numel(pValuesMB);
        end

        % ---- Rect Grid ----
        zScores_RGs  = zScores_RG( pvalsStimSelected <= params.threshold);
        spkR_RGs     = spkR_RG(   pvalsStimSelected <= params.threshold);
        spkDiff_RGs  = spkDiff_RG(pvalsStimSelected <= params.threshold);
        pvals_RG     = pValuesRG( pvalsStimSelected <= params.threshold);

        zScores_RGg  = zScores_RG( pValuesRG <= params.threshold);
        sumNeurRG    = numel(zScores_RGg);
        spkR_RGg     = spkR_RG(    pValuesRG <= params.threshold);
        spkDiff_RGg  = spkDiff_RG( pValuesRG <= params.threshold);
        respIndexes  = [respIndexes, find(pValuesRG <= params.threshold)];

        try
            idx = (longTable.insertion == categorical(j)) & ...
                  (longTable.stimulus  == categorical("RG"));
            longTable.respNeur(idx)      = sumNeurRG;
            longTable.totalSomaticN(idx) = numel(pValuesMB);  % total = same for all rows
        end

        % If RG was not recorded, overwrite with -inf sentinel
        % SUGG-2: NaN is safer than -inf for absent data
        if isequal(params.StimsPresent{2},'')
            zScores_RGs = zScores_RG  - inf;
            spkR_RGs    = zScores_RG  - inf;
            spkDiff_RGs = zScores_RG  - inf;
            pvals_RG    = zScores_RG  - inf;
            sumNeurRG   = 0;
            zScores_RGg = zScores_RGg - inf;
            spkR_RGg    = spkR_RGg    - inf;
            spkDiff_RGg = spkDiff_RGg - inf;
        end

        % ---- Moving Bar ----
        zScores_MBRs = zScores_MBR( pvalsStimSelected <= params.threshold);
        spkR_MBRs    = spkR_MBR(   pvalsStimSelected <= params.threshold);
        spkDiff_MBRs = spkDiff_MBR(pvalsStimSelected <= params.threshold);
        pvals_MBR    = pValuesMBR( pvalsStimSelected <= params.threshold);

        zScores_MBRg = zScores_MBR( pValuesMBR <= params.threshold);
        sumNeurMBR   = numel(zScores_MBRg);
        spkR_MBRg    = spkR_MBR(    pValuesMBR <= params.threshold);
        spkDiff_MBRg = spkDiff_MBR( pValuesMBR <= params.threshold);
        respIndexes  = [respIndexes, find(pValuesMBR <= params.threshold)];

        try
            idx = (longTable.insertion == categorical(j)) & ...
                  (longTable.stimulus  == categorical("MBR"));
            longTable.respNeur(idx)      = sumNeurMBR;
            longTable.totalSomaticN(idx) = numel(pValuesMB);
        end

        if isequal(params.StimsPresent{3},'')
            zScores_MBRs = zScores_MBRs - inf;
            spkR_MBRs    = zScores_MBRs - inf;  % NOTE: uses already -inf'd zscores
            spkDiff_MBRs = zScores_MBRs - inf;
            pvals_MBR    = zScores_MBRs - inf;
            sumNeurMBR   = 0;
            zScores_MBRg = zScores_MBRg - inf;
            spkR_MBRg    = zScores_MBRg - inf;
            spkDiff_MBRg = zScores_MBRg - inf;
        end

        % ---- Gratings (moving and static) ----
        zScores_SDGms = zScores_SDGm( pvalsStimSelected <= params.threshold);
        spkR_SDGms    = spkR_SDGm(   pvalsStimSelected <= params.threshold);
        spkDiff_SDGms = spkDiff_SDGm(pvalsStimSelected <= params.threshold);
        pvals_SDGm    = pValuesSDGm( pvalsStimSelected <= params.threshold);

        zScores_SDGss = zScores_SDGs( pvalsStimSelected <= params.threshold);
        spkR_SDGss    = spkR_SDGs(   pvalsStimSelected <= params.threshold);
        spkDiff_SDGss = spkDiff_SDGs(pvalsStimSelected <= params.threshold);
        pvals_SDGs    = pValuesSDGs( pvalsStimSelected <= params.threshold);

        zScores_SDGmg = zScores_SDGm( pValuesSDGm <= params.threshold);
        sumNeurSDGm   = numel(zScores_SDGmg);
        spkR_SDGmg    = spkR_SDGm(    pValuesSDGm <= params.threshold);
        spkDiff_SDGmg = spkDiff_SDGm( pValuesSDGm <= params.threshold);
        respIndexes   = [respIndexes, find(pValuesSDGm <= params.threshold)];

        zScores_SDGsg = zScores_SDGs( pValuesSDGs <= params.threshold);
        sumNeurSDGs   = numel(zScores_SDGsg);
        spkR_SDGsg    = spkR_SDGs(    pValuesSDGs <= params.threshold);
        spkDiff_SDGsg = spkDiff_SDGs( pValuesSDGs <= params.threshold);
        respIndexes   = [respIndexes, find(pValuesSDGs <= params.threshold)];

        try
            idx = (longTable.insertion == categorical(j)) & ...
                  (longTable.stimulus  == categorical("SDGm"));
            longTable.respNeur(idx)      = sumNeurSDGm;
            longTable.totalSomaticN(idx) = numel(pValuesMB);
        end
        try
            idx = (longTable.insertion == categorical(j)) & ...
                  (longTable.stimulus  == categorical("SDGs"));
            longTable.respNeur(idx)      = sumNeurSDGs;
            longTable.totalSomaticN(idx) = numel(pValuesMB);
        end

        if isequal(params.StimsPresent{4},'')
            zScores_SDGss = zScores_SDGss - inf;
            spkR_SDGss    = spkR_SDGss    - inf;
            spkDiff_SDGss = spkDiff_SDGss - inf;
            pvals_SDGs    = pvals_SDGs    - inf;

            zScores_SDGms = zScores_SDGms - inf;
            spkR_SDGms    = spkR_SDGms    - inf;
            spkDiff_SDGms = spkDiff_SDGms - inf;
            pvals_SDGm    = pvals_SDGm    - inf;

            % BUG-4: sumNeurSDG (new var) is set to 0 here, but
            %        sumNeurSDGm and sumNeurSDGs are NOT reset to 0.
            %        sumNeurSDGmt{j} and sumNeurSDGst{j} below will then
            %        store stale values from the previous iteration.
            %        FIX: replace the line below with:
            %             sumNeurSDGm = 0; sumNeurSDGs = 0;
            sumNeurSDGm = 0;   % FIX applied (was: sumNeurSDG = 0)
            sumNeurSDGs = 0;   % FIX applied

            zScores_SDGmg = zScores_SDGmg - inf;
            spkR_SDGmg    = zScores_SDGmg - inf;
            spkDiff_SDGmg = zScores_SDGmg - inf;

            zScores_SDGsg = zScores_SDGsg - inf;
            spkR_SDGsg    = zScores_SDGsg - inf;
            spkDiff_SDGsg = zScores_SDGsg - inf;
        end

        % ---- Full-Field Flash ----
        zScores_FFFs = zScores_FFF( pvalsStimSelected <= params.threshold);
        spkR_FFFs    = spkR_FFF(   pvalsStimSelected <= params.threshold);
        spkDiff_FFFs = spkDiff_FFF(pvalsStimSelected <= params.threshold);
        pvals_FFF    = pValuesFFF( pvalsStimSelected <= params.threshold);

        zScores_FFFg = zScores_FFF( pValuesFFF <= params.threshold);
        sumNeurFFF   = numel(zScores_FFFg);
        spkR_FFFg    = spkR_FFF(    pValuesFFF <= params.threshold);
        spkDiff_FFFg = spkDiff_FFF( pValuesFFF <= params.threshold);
        respIndexes  = [respIndexes, find(pValuesFFF <= params.threshold)];

        try
            idx = (longTable.insertion == categorical(j)) & ...
                  (longTable.stimulus  == categorical("FFF"));
            longTable.respNeur(idx)      = sumNeurFFF;
            longTable.totalSomaticN(idx) = numel(pValuesMB);
        end

        if isequal(params.StimsPresent{7},'')
            zScores_FFFs = zScores_FFFs - inf;
            spkR_FFFs    = spkR_FFFs    - inf;
            spkDiff_FFFs = spkDiff_FFFs - inf;
            pvals_FFF    = pvals_FFF    - inf;
            sumNeurFFF   = 0;
            zScores_FFFg = zScores_FFFg - inf;
            spkR_FFFg    = zScores_FFFg - inf;
            spkDiff_FFFg = zScores_FFFg - inf;
        end

        % ---- Natural Images ----
        zScores_NIs  = zScores_NI( pvalsStimSelected <= params.threshold);
        spkR_NIs     = spkR_NI(   pvalsStimSelected <= params.threshold);
        spkDiff_NIs  = spkDiff_NI(pvalsStimSelected <= params.threshold);
        pvals_NI     = pValuesNI( pvalsStimSelected <= params.threshold);

        zScores_NIg  = zScores_NI( pValuesNI <= params.threshold);
        sumNeurNI    = numel(zScores_NIg);
        spkR_NIg     = spkR_NI(    pValuesNI <= params.threshold);
        spkDiff_NIg  = spkDiff_NI( pValuesNI <= params.threshold);
        respIndexes  = [respIndexes, find(pValuesNI <= params.threshold)];

        try
            idx = (longTable.insertion == categorical(j)) & ...
                  (longTable.stimulus  == categorical("NI"));
            longTable.respNeur(idx)      = sumNeurNI;
            longTable.totalSomaticN(idx) = numel(pValuesMB);
        end

        if isequal(params.StimsPresent{5},'')
            zScores_NIs = zScores_NIs - inf;
            spkR_NIs    = spkR_NIs    - inf;
            spkDiff_NIs = spkDiff_NIs - inf;
            pvals_NI    = pvals_NI    - inf;
            sumNeurNI   = 0;
            zScores_NIg = zScores_NIg - inf;
            spkR_NIg    = zScores_NIg - inf;
            spkDiff_NIg = zScores_NIg - inf;
        end

        % ---- Natural Video ----
        zScores_NVs  = zScores_NV( pvalsStimSelected <= params.threshold);
        spkR_NVs     = spkR_NV(   pvalsStimSelected <= params.threshold);
        spkDiff_NVs  = spkDiff_NV(pvalsStimSelected <= params.threshold);
        pvals_NV     = pValuesNV( pvalsStimSelected <= params.threshold);

        zScores_NVg  = zScores_NV( pValuesNV <= params.threshold);
        sumNeurNV    = numel(zScores_NVg);
        spkR_NVg     = spkR_NV(    pValuesNV <= params.threshold);
        spkDiff_NVg  = spkDiff_NV( pValuesNV <= params.threshold);
        respIndexes  = [respIndexes, find(pValuesNV <= params.threshold)];

        try
            idx = (longTable.insertion == categorical(j)) & ...
                  (longTable.stimulus  == categorical("NV"));
            longTable.respNeur(idx)      = sumNeurNV;
            longTable.totalSomaticN(idx) = numel(pValuesMB);
        end

        if isequal(params.StimsPresent{6},'')
            zScores_NVs = zScores_NVs - inf;
            spkR_NVs    = spkR_NVs    - inf;
            spkDiff_NVs = spkDiff_NVs - inf;
            pvals_NV    = pvals_NV    - inf;
            sumNeurNV   = 0;
            zScores_NVg = zScores_NVg - inf;
            spkR_NVg    = zScores_NVg - inf;
            spkDiff_NVg = zScores_NVg - inf;
        end

        % Union of all neuron indices responsive to at least one stimulus
        responsiveNeuronsj = unique(respIndexes);

        % BUG-5: `2+2` is a debug breakpoint stub — removed here.
        %        Replace with a proper warning:
        if numel(zScores_NVs) ~= numel(zScores_NIs)
            warning('PlotZScoreComparison: NV and NI filtered vectors have different lengths in experiment %d.', ex);
        end

        % ------------------------------------------------------------------
        % 4j – Re-extract animal and insertion labels (fresh regex in case
        %      the object was re-created above)
        % ------------------------------------------------------------------

        Animal    = string(regexp(vs.getAnalysisFileName, 'PV\d+', 'match', 'once'));
        Insertion = regexp(vs.getAnalysisFileName, 'Insertion\d+', 'match', 'once');
        Insertion = str2double(regexp(Insertion, '\d+', 'match'));

        % Fallback: some animals use 'SA##' naming convention
        if isequal(Animal, "")
            Animal = string(regexp(vs.getAnalysisFileName, 'SA\d+', 'match', 'once'));
        end

        % BUG-3: AnimalI is updated inside the first if-block, so the second
        %        if-block (checking Animal~=AnimalI for insertion counting) 
        %        always sees them as equal after the first block runs.
        %        FIX: capture the old value before updating.
        AnimalChanged = (Animal ~= AnimalI);   % evaluate BEFORE updating AnimalI

        if AnimalChanged
            animal = animal + 1;              % new animal encountered
            AnimalNames{animal} = Animal;     % store its name
            AnimalI = Animal;                 % update tracker
        end

        % Count a new insertion if the insertion number changed OR a new animal
        if Insertion ~= InsertionI || AnimalChanged   % FIX: use pre-evaluated flag
            InsertionI = Insertion;
            insertion  = insertion + 1;
        end

        % ------------------------------------------------------------------
        % 4k – Store this experiment's data into per-experiment cell arrays
        % ------------------------------------------------------------------

        % Replicate animal/insertion IDs to match number of anchor-filtered neurons
        animalVector{j}    = repmat(animal,    [1, numel(zScores_MBs)]);
        insertionVector{j} = repmat(insertion, [1, numel(zScores_MBs)]);

        % Anchor-filtered data (neurons significant for the anchor stimulus)
        zScoresMB{j}    = zScores_MBs;
        zScoresRG{j}    = zScores_RGs;
        pvalsRG{j}      = pvals_RG;
        sumNeurRGt{j}   = sumNeurRG;
        pvalsMB{j}      = pvals_MB;
        sumNeurMBt{j}   = sumNeurMB;
        spKrMB{j}       = spkR_MBs';
        spKrRG{j}       = spkR_RGs';
        diffSpkMB{j}    = spkDiff_MBs;
        diffSpkRG{j}    = spkDiff_RGs;

        zScoresFFF{j}   = zScores_FFFs;
        spKrFFF{j}      = spkR_FFFs';
        diffSpkFFF{j}   = spkDiff_FFFs;
        pvalsFFF{j}     = pvals_FFF;
        sumNeurFFFt{j}  = sumNeurFFF;

        zScoresMBR{j}   = zScores_MBRs;
        spKrMBR{j}      = spkR_MBRs';
        diffSpkMBR{j}   = spkDiff_MBRs;
        pvalsMBR{j}     = pvals_MBR;
        sumNeurMBRt{j}  = sumNeurMBR;

        zScoresSDGm{j}  = zScores_SDGms;
        spKrSDGm{j}     = spkR_SDGms';
        diffSpkSDGm{j}  = spkDiff_SDGms;
        pvalsSDGm{j}    = pvals_SDGm;
        sumNeurSDGmt{j} = sumNeurSDGm;

        zScoresSDGs{j}  = zScores_SDGss;
        spKrSDGs{j}     = spkR_SDGss';
        diffSpkSDGs{j}  = spkDiff_SDGss;
        pvalsSDGs{j}    = pvals_SDGs;
        sumNeurSDGst{j} = sumNeurSDGs;

        zScoresNI{j}    = zScores_NIs;
        spKrNI{j}       = spkR_NIs';
        diffSpkNI{j}    = spkDiff_NIs;
        pvalsNI{j}      = pvals_NI;
        sumNeurNIt{j}   = sumNeurNI;

        zScoresNV{j}    = zScores_NVs;
        spKrNV{j}       = spkR_NVs';
        diffSpkNV{j}    = spkDiff_NVs;
        pvalsNV{j}      = pvals_NV;
        sumNeurNVt{j}   = sumNeurNV;

        % Self-responsive data (neurons significant for EACH respective stimulus)
        zScoresMBg{j}   = zScores_MBg;   spkRMBg{j}   = spkR_MBg;   spkDiffMBg{j}   = spkDiff_MBg;
        zScoresRGg{j}   = zScores_RGg;   spkRRGg{j}   = spkR_RGg;   spkDiffRGg{j}   = spkDiff_RGg;
        zScoresMBRg{j}  = zScores_MBRg;  spkRMBRg{j}  = spkR_MBRg;  spkDiffMBRg{j}  = spkDiff_MBRg;
        zScoresSDGmg{j} = zScores_SDGmg; spkRSDGmg{j} = spkR_SDGmg; spkDiffSDGmg{j} = spkDiff_SDGmg;
        zScoresSDGsg{j} = zScores_SDGsg; spkRSDGsg{j} = spkR_SDGsg; spkDiffSDGsg{j} = spkDiff_SDGsg;
        zScoresFFFg{j}  = zScores_FFFg;  spkRFFFg{j}  = spkR_FFFg;  spkDiffFFFg{j}  = spkDiff_FFFg;
        zScoresNIg{j}   = zScores_NIg;   spkRNIg{j}   = spkR_NIg;   spkDiffNIg{j}   = spkDiff_NIg;
        zScoresNVg{j}   = zScores_NVg;   spkRNVg{j}   = spkR_NVg;   spkDiffNVg{j}   = spkDiff_NVg;

        % Set of neuron indices responsive to at least one stimulus in this recording
        responsiveNeurons{j} = responsiveNeuronsj;

        j = j + 1;   % advance experiment counter

        fprintf('Finished recording: %s .\n', NP.recordingName)

    end  % end for ex = expList

    % =========================================================================
    % SECTION 5 – PACK ALL DATA INTO STRUCT S AND SAVE
    % =========================================================================

    % Anchor-filtered values (neurons responsive to the first Stims2Comp element)
    S.stimValsSignif2oneStim.spKrMB    = spKrMB;
    S.stimValsSignif2oneStim.spKrRG    = spKrRG;
    S.stimValsSignif2oneStim.diffSpkMB = diffSpkMB;
    S.stimValsSignif2oneStim.diffSpkRG = diffSpkRG;
    S.stimValsSignif2oneStim.zScoresMB = zScoresMB;
    S.stimValsSignif2oneStim.zScoresRG = zScoresRG;
    S.pvals.pvalsMB = pvalsMB;
    S.pvals.pvalsRG = pvalsRG;

    S.stimValsSignif2oneStim.spKrMBR    = spKrMBR;
    S.stimValsSignif2oneStim.spKrFFF    = spKrFFF;
    S.stimValsSignif2oneStim.diffSpkMBR = diffSpkMBR;
    S.stimValsSignif2oneStim.diffSpkFFF = diffSpkFFF;
    S.stimValsSignif2oneStim.zScoresMBR = zScoresMBR;
    S.stimValsSignif2oneStim.zScoresFFF = zScoresFFF;
    S.pvals.pvalsFFF = pvalsFFF;
    S.pvals.pvalsMBR = pvalsMBR;

    S.stimValsSignif2oneStim.spKrSDGm    = spKrSDGm;
    S.stimValsSignif2oneStim.spKrSDGs    = spKrSDGs;
    S.stimValsSignif2oneStim.diffSpkSDGm = diffSpkSDGm;
    S.stimValsSignif2oneStim.diffSpkSDGs = diffSpkSDGs;
    S.stimValsSignif2oneStim.zScoresSDGm = zScoresSDGm;
    S.stimValsSignif2oneStim.zScoresSDGs = zScoresSDGs;
    S.pvals.pvalsSDGm = pvalsSDGm;
    S.pvals.pvalsSDGs = pvalsSDGs;

    S.stimValsSignif2oneStim.spKrNI    = spKrNI;
    S.stimValsSignif2oneStim.spKrNV    = spKrNV;
    S.stimValsSignif2oneStim.diffSpkNI = diffSpkNI;
    S.stimValsSignif2oneStim.diffSpkNV = diffSpkNV;
    S.stimValsSignif2oneStim.zScoresNI = zScoresNI;
    S.stimValsSignif2oneStim.zScoresNV = zScoresNV;
    S.pvals.pvalsNI = pvalsNI;
    S.pvals.pvalsNV = pvalsNV;

    % Self-responsive values (each neuron counted only for its own stimulus)
    S.stimValsSignif.zScoresMBg   = zScoresMBg;   S.stimValsSignif.spkRMBg   = spkRMBg;   S.stimValsSignif.spkDiffMBg   = spkDiffMBg;
    S.stimValsSignif.zScoresRGg   = zScoresRGg;   S.stimValsSignif.spkRRGg   = spkRRGg;   S.stimValsSignif.spkDiffRGg   = spkDiffRGg;
    S.stimValsSignif.zScoresMBRg  = zScoresMBRg;  S.stimValsSignif.spkRMBRg  = spkRMBRg;  S.stimValsSignif.spkDiffMBRg  = spkDiffMBRg;
    S.stimValsSignif.zScoresSDGmg = zScoresSDGmg; S.stimValsSignif.spkRSDGmg = spkRSDGmg; S.stimValsSignif.spkDiffSDGmg = spkDiffSDGmg;
    S.stimValsSignif.zScoresSDGsg = zScoresSDGsg; S.stimValsSignif.spkRSDGsg = spkRSDGsg; S.stimValsSignif.spkDiffSDGsg = spkDiffSDGsg;
    S.stimValsSignif.zScoresFFFg  = zScoresFFFg;  S.stimValsSignif.spkRFFFg  = spkRFFFg;  S.stimValsSignif.spkDiffFFFg  = spkDiffFFFg;
    S.stimValsSignif.zScoresNIg   = zScoresNIg;   S.stimValsSignif.spkRNIg   = spkRNIg;   S.stimValsSignif.spkDiffNIg   = spkDiffNIg;
    S.stimValsSignif.zScoresNVg   = zScoresNVg;   S.stimValsSignif.spkRNVg   = spkRNVg;   S.stimValsSignif.spkDiffNVg   = spkDiffNVg;

    % Responsive neuron counts per insertion per stimulus
    S.stimValsSignif.sumNeurMB   = sumNeurMBt;
    S.stimValsSignif.sumNeurRG   = sumNeurRGt;
    S.stimValsSignif.sumNeurMBR  = sumNeurMBRt;
    S.stimValsSignif.sumNeurSDGm = sumNeurSDGmt;
    S.stimValsSignif.sumNeurSDGs = sumNeurSDGst;
    S.stimValsSignif.sumNeurFFF  = sumNeurFFFt;
    S.stimValsSignif.sumNeurNI   = sumNeurNIt;
    S.stimValsSignif.sumNeurNV   = sumNeurNVt;

    % Metadata and indexing
    S.expList          = expList;          % experiment IDs processed
    S.animalVector     = animalVector;     % per-neuron animal index
    S.insertionVector  = insertionVector;  % per-neuron insertion index
    S.totalUnits       = totalU;           % total unit count per experiment
    S.params           = params;           % parameter snapshot
    S.responsiveNeurons = responsiveNeurons; % union-responsive neuron indices
    S.TableRespNeurs   = longTable;        % fraction-responsive table
    S.TableStimComp    = longTablePairComp; % pairwise z-score/SpkR table

    save([saveDir nameOfFile], '-struct', 'S');   % save struct fields as top-level variables

end  % end if forloop

% =========================================================================
% SECTION 6 – PAIRWISE COMPARISON (ComparePairs mode)
% =========================================================================

if ~isempty(params.ComparePairs)

    pairs = params.ComparePairs;  % cell of stimulus name(s) to compare

    % -----------------------------------------------------------------------
    % BUG-1 FIX: Guard against empty pairwise table (no significant units
    %             found in any experiment).  splitapply on an empty grouping
    %             vector throws an error.
    % -----------------------------------------------------------------------
    if isempty(S.TableStimComp) || height(S.TableStimComp) == 0
        warning('PlotZScoreComparison:noUnits', ...
            ['No significant units found for pairwise comparison of %s vs %s.\n' ...
             'Returning empty figure.'], pairs{1}, pairs{2});
        fig = figure;   % return empty figure handle to satisfy output contract
        return
    end

    % Replace NaN z-scores / spike rates with 0 (conservative: treat as no response)
    S.TableStimComp.('Z-score')(isnan(S.TableStimComp.('Z-score'))) = 0;
    S.TableStimComp.SpkR(isnan(S.TableStimComp.SpkR)) = 0;

    % Find insertions that contain both stimuli in the pair
    [G, ~] = findgroups(S.TableStimComp.insertion);
    hasAll = splitapply(@(s) all(ismember(unique(categorical(pairs)), s)), ...
                        S.TableStimComp.stimulus, G);

    % Restrict table to complete insertions (have both stimuli) and relevant rows
    tempTable = S.TableStimComp( ...
        hasAll(G) & ismember(S.TableStimComp.stimulus, unique(categorical(pairs))), :);

    nBoot = 10000;   % number of hierarchical bootstrap iterations

    % SHARED COLORMAP: built once, reused in every swarm and scatter panel.
    % double() on a categorical returns the rank within categories(), which is
    % the same ordering used to index into the colormap — guaranteeing that
    % animal X gets identical RGB in the swarm and in both scatter plots.
    animalOrder  = categories(S.TableStimComp.animal);   % canonical ordering
    nAnimals     = numel(animalOrder);
    sharedCmap   = lines(nAnimals);                       % nAnimals × 3 RGB matrix
    animalIdxAll = double(S.TableStimComp.animal);

    % Pre-compute the row masks for pairs{1} and pairs{2} — used in both
    % the Z-score and spike-rate scatter panels below.
    mask1 = S.TableStimComp.stimulus == pairs{1};
    mask2 = S.TableStimComp.stimulus == pairs{2};
    cIdx  = animalIdxAll(mask1);   % colour index aligned with pair{1} / pair{2} rows

    % -----------------------------------------------------------------------
    % 6a – Z-score comparison via hierarchical bootstrapping
    % -----------------------------------------------------------------------

    j  = 1;
    ps = zeros(1, size(pairs, 1));   % one p-value per stimulus pair

    for i = 1:size(pairs, 1)

        diffs   = [];  % per-neuron differences (stim1 – stim2) pooled across insertions
        insers  = [];  % insertion label for each difference
        animals = [];  % animal label for each difference

        for ins = unique(S.TableStimComp.insertion)'

            % Select rows for this insertion × each stimulus
            idx1 = S.TableStimComp.insertion == categorical(ins) & ...
                   S.TableStimComp.stimulus  == pairs{j,1};
            idx2 = S.TableStimComp.insertion == categorical(ins) & ...
                   S.TableStimComp.stimulus  == pairs{j,2};

            V1 = S.TableStimComp.('Z-score')(idx1);
            V2 = S.TableStimComp.('Z-score')(idx2);

            % Unique animal for this insertion (should be exactly one)
            animal = unique(S.TableStimComp.animal(idx1));

            % Append per-neuron differences and labels
            diffs   = [diffs;   V1 - V2];
            insers  = [insers;  double(repmat(ins,    size(V1,1), 1))];
            animals = [animals; double(repmat(animal, size(V1,1), 1))];
        end

        % Hierarchical bootstrap: resample at animal level, then insertion level
        bootDiff = hierBoot(diffs, nBoot, insers, animals);
        ps(j) = mean(bootDiff <= 0);   % p-value: proportion of bootstrap samples ≤ 0
        j = j + 1;
    end

    ZscoreYlimUp = ceil(max(S.TableStimComp.("Z-score")))+4;

    % Swarm plot with bootstrap-derived significance (returns subsampling index)
    [fig, randiColors] = plotSwarmBootstrapWithComparisons(S.TableStimComp, pairs, ps, ...
        {'Z-score'}, yLegend='Z-score', yMaxVis=ZscoreYlimUp, diff=true, plotMeanSem=false, Alpha=0.7);

    ax = gca;
    ax.YAxis.FontSize = 8;   ax.YAxis.FontName = 'helvetica';
    ax.XAxis.FontSize = 8;   ax.XAxis.FontName = 'helvetica';

    set(fig, 'Units', 'centimeters', 'Position', [20 20 10 6]);
    colormap(fig, sharedCmap);   % enforce shared colormap so swarm colours match scatter

    % Reload analysis object for figure saving (path extraction)
    NP = loadNPclassFromTable(expList(1));
    vs = linearlyMovingBallAnalysis(NP);

    ylims = ylim;

    if params.PaperFig
        vs.printFig(fig, sprintf('Zcore-comparison-Swarm-%s-%s', ...
            params.ComparePairs{1}, params.ComparePairs{2}), PaperFig=params.PaperFig);
    end

    % -----------------------------------------------------------------------
    % 6b – Scatter plot: first vs second stimulus in pairs (Z-score)
    %      SUGG-5: randiColors is a subsampling index from the swarm function.
    %              If it subsamples non-uniformly, the scatter may misrepresent
    %              the data density.  Consider plotting all points for publication.
    % -----------------------------------------------------------------------

    fig = figure;

    pair1 = S.TableStimComp.("Z-score")(mask1);
    pair2 = S.TableStimComp.("Z-score")(mask2);
    % cIdx already computed above — direct RGB lookup, no implicit categorical conversion

    % Scatter with animal-coded colour, using subsampled indices
    scatter(pair1, pair2, 7, sharedCmap(cIdx,:), ...
        "filled", "MarkerFaceAlpha", 0.3)
    hold on
    axis equal

    lims = [min(S.TableStimComp.("Z-score")), max(S.TableStimComp.("Z-score"))];
    plot(lims, lims, 'k--', 'LineWidth', 1.5)   % identity line
    ylim(lims); xlim(lims)

    % Convert internal stimulus abbreviations to display labels
    s = string(pairs);
    s = replace(s, "RG",   "SB");   % Rect Grid → Square Ball
    s = replace(s, "SDGs", "SG");   % static gratings label
    s = replace(s, "SDGm", "MG");   % moving gratings label

    xlabel(s{1}); ylabel(s{2})
    colormap(fig, sharedCmap)

    ax = gca;
    ax.YAxis.FontSize = 8; ax.YAxis.FontName = 'helvetica';
    ax.XAxis.FontSize = 8; ax.XAxis.FontName = 'helvetica';
    set(fig, 'Units', 'centimeters', 'Position', [20 20 5 5]);
    title('Z-score')

    if params.PaperFig
        vs.printFig(fig, sprintf('Zcore-comparison-Scatter-%s-%s', ...
            params.ComparePairs{1}, params.ComparePairs{2}), PaperFig=params.PaperFig);
    end

    % -----------------------------------------------------------------------
    % 6c – Spike-rate comparison via hierarchical bootstrapping
    % -----------------------------------------------------------------------

    j  = 1;
    ps = zeros(1, size(pairs, 1));

    for i = 1:size(pairs, 1)

        diffs   = [];
        insers  = [];
        animals = [];

        for ins = unique(S.TableStimComp.insertion)'

            idx1 = S.TableStimComp.insertion == categorical(ins) & ...
                   S.TableStimComp.stimulus  == pairs{j,1};
            idx2 = S.TableStimComp.insertion == categorical(ins) & ...
                   S.TableStimComp.stimulus  == pairs{j,2};

            V1 = S.TableStimComp.SpkR(idx1);
            V2 = S.TableStimComp.SpkR(idx2);

            animal  = unique(S.TableStimComp.animal(idx1));
            diffs   = [diffs;   V1 - V2];
            insers  = [insers;  double(repmat(ins,    size(V1,1), 1))];
            animals = [animals; double(repmat(animal, size(V1,1), 1))];
        end

        bootDiff = hierBoot(diffs, nBoot, insers, animals);
        ps(j)    = mean(bootDiff <= 0);
        j        = j + 1;
    end

    V1max = max(diffs);   % use max observed difference to set y-axis ceiling

    [fig, randiColors] = plotSwarmBootstrapWithComparisons(S.TableStimComp, pairs, ps, ...
        {'SpkR'}, yLegend='SpkR', yMaxVis=V1max, diff=true, plotMeanSem=false, Alpha=0.7);

    ax = gca;
    ax.YAxis.FontSize = 8; ax.YAxis.FontName = 'helvetica';
    ax.XAxis.FontSize = 8; ax.XAxis.FontName = 'helvetica';
    colormap(fig, sharedCmap);
    set(fig, 'Units', 'centimeters', 'Position', [20 20 10 6]);

    if params.PaperFig
        vs.printFig(fig, sprintf('spkRate-comparison-Swarm-%s-%s', ...
            params.ComparePairs{1}, params.ComparePairs{2}), PaperFig=params.PaperFig);
    end

    % -----------------------------------------------------------------------
    % 6d – Scatter plot: first vs second stimulus (Spike Rate)
    % -----------------------------------------------------------------------

    fig = figure;
    pair1 = S.TableStimComp.SpkR(mask1);   % mask1 pre-computed above
    pair2 = S.TableStimComp.SpkR(mask2);
    scatter(pair1, pair2, 7, sharedCmap(cIdx,:), ...
        "filled", "MarkerFaceAlpha", 0.3)
    hold on
    axis equal

    lims = [min(S.TableStimComp.SpkR), max(S.TableStimComp.SpkR)];
    plot(lims, lims, 'k--', 'LineWidth', 1.5)
    ylim(lims); xlim(lims)

    xlabel(s{1}); ylabel(s{2})
    colormap(fig, sharedCmap)

    ax = gca;
    ax.YAxis.FontSize = 8; ax.YAxis.FontName = 'helvetica';
    ax.XAxis.FontSize = 8; ax.XAxis.FontName = 'helvetica';
    set(fig, 'Units', 'centimeters', 'Position', [20 20 5 5]);
    title('Spk. rate')

    if params.PaperFig
        vs.printFig(fig, sprintf('spkRate-comparison-Scatter-%s-%s', ...
            params.ComparePairs{1}, params.ComparePairs{2}), PaperFig=params.PaperFig);
    end

else
    % =========================================================================
    % SECTION 7 – MULTI-STIMULUS OVERVIEW (non-pairwise mode)
    %             Compares ALL stimuli in Stims2Comp using swarm + scatter.
    % =========================================================================

    fig = figure;
    tiledlayout(2, 2, "TileSpacing", "compact");

    % Choose field-name set based on whether each-stim or anchor-filtered
    if ~params.EachStimSignif
        fn  = fieldnames(S.stimValsSignif2oneStim);   % anchor-filtered fields
    else
        fn  = fieldnames(S.stimValsSignif);            % self-responsive fields
    end
    fnp = fieldnames(S.pvals);

    % Expand 'SDG' shorthand into two separate entries (moving + static)
    Stims2Comp2 = {};
    for i = 1:numel(Stims2Comp)
        if strcmp(Stims2Comp{i}, 'SDG')
            Stims2Comp2 = [Stims2Comp2, {'SDGs','SDGm'}];
        else
            Stims2Comp2 = [Stims2Comp2, Stims2Comp(i)];
        end
    end

    % Select suffix used in field-name lookup
    endingOpts = {'','g'};   % '' = anchor-filtered suffix, 'g' = self-responsive
    ending2    = endingOpts{1 + params.EachStimSignif};

    % Pre-allocate arrays that will hold concatenated data for each stimulus
    StimZS    = cell(numel(Stims2Comp2), 1);   % z-scores per stimulus
    stimRSP   = cell(numel(Stims2Comp2), 1);   % spike rates per stimulus
    stimPvals = cell(numel(Stims2Comp2), 1);   % p-values per stimulus
    x         = [];   % stimulus-index label for each neuron (for swarmchart x-axis)

    for i = 1:numel(Stims2Comp2)

        ending  = Stims2Comp2{i};   % e.g. 'MB', 'RGg', …
        % Regex: field names starting with 'zS' and ending with the stimulus tag
        pattern = ['^zS.*' ending ending2 '$'];
        matches = fn(~cellfun('isempty', regexp(fn, pattern)));

        % Concatenate z-scores across experiments
        if ~params.EachStimSignif
            StimZS{i} = cell2mat(S.stimValsSignif2oneStim.(matches{1}))';
        else
            StimZS{i} = cell2mat(S.stimValsSignif.(matches{1}))';
        end

        % Build pattern for spike rate OR spike difference (diffResp flag)
        if ~params.diffResp
            pattern = ['^spKr.*' ending ending2 '$'];
        else
            pattern = ['^diffSpk.*' ending ending2 '$'];
        end

        matches = fn(~cellfun('isempty', regexp(fn, pattern)));

        if params.EachStimSignif
            matches   = fn(~cellfun('isempty', regexp(fn, pattern, 'ignorecase')));
            C         = S.stimValsSignif.(matches{1});
            C         = cellfun(@(x) x', C, 'UniformOutput', false);
            stimRSP{i} = cell2mat(C');
        else
            % Try several concatenation strategies to handle shape inconsistencies
            try
                stimRSP{i} = cell2mat(S.stimValsSignif2oneStim.(matches{1})');
            catch
                try
                    stimRSP{i} = cell2mat(S.stimValsSignif2oneStim.(matches{1}));
                catch
                    % Last resort: force column, then vertcat
                    Ccol       = cellfun(@(x) x(:), S.stimValsSignif2oneStim.(matches{1}), ...
                                         'UniformOutput', false);
                    stimRSP{i} = vertcat(Ccol{:})';
                end
            end
        end

        % Retrieve p-values for this stimulus
        pattern      = ['^pvals.*' ending '$'];
        matches      = fnp(~cellfun('isempty', regexp(fnp, pattern)));
        stimPvals{i} = cell2mat(S.pvals.(matches{1}))';

        % Build x-axis labels: all neurons for stimulus i get label i
        x = [x; ones(size(StimZS{i})) * i];

    end

    % Per-neuron animal and insertion index vectors (from anchor-filtered pool)
    AnIndex      = cell2mat(S.animalVector)';
    InsIndex     = cell2mat(S.insertionVector)';
    colormapUsed = parula(max(AnIndex)) .* 0.6;   % muted parula for animal colouring

    % -----------------------------------------------------------------------
    % 7a – Z-score swarm chart
    % -----------------------------------------------------------------------

    y = cell2mat(StimZS);   % all z-scores concatenated (length = total neurons × stims)

    allColorIndices = repmat(AnIndex, numel(Stims2Comp2), 1);   % replicate animal index

    nexttile
    if ~params.EachStimSignif
        swarmchart(x, y, 5, colormapUsed(allColorIndices,:), 'filled', 'MarkerFaceAlpha', 0.7);
    else
        swarmchart(x, y, 5, 'filled', 'MarkerFaceAlpha', 0.7);
    end
    xticks(1:8);
    xticklabels(Stims2Comp2);
    ylabel('Z-score');
    set(fig, 'Color', 'w')
    yline(0, 'LineWidth', 2)   % reference line at zero
    ylim([-5 40])

    % -----------------------------------------------------------------------
    % 7b – Hierarchical bootstrapping for Z-score group comparison
    %      (computed fresh or loaded from saved S.groupStats)
    % -----------------------------------------------------------------------

    if params.overwriteGroupStats || ~isfield(S, 'groupStats')

        % Bootstrap the first (anchor) stimulus
        FirstStim = y(x == 1);
        BootFirst = hierBoot(FirstStim(~isnan(FirstStim)), 10000, ...
                             InsIndex(~isnan(FirstStim)), AnIndex(~isnan(FirstStim)));

        j = 1;
        for i = 2:numel(Stims2Comp2)
            secondaryStim = y(x == i);
            secondaryStim(isnan(secondaryStim)) = 0;   % treat NaN as no response
            validMask     = secondaryStim ~= -inf;
            secondaryStim = secondaryStim(validMask);

            BootSec  = hierBoot(secondaryStim, 10000, InsIndex(validMask), AnIndex(validMask));
            probs{j} = get_direct_prob(BootFirst, BootSec);   % Bayesian-style overlap probability
            ps{j}    = mean(BootSec >= BootFirst);             % frequentist p-value
            j        = j + 1;
        end

        S.groupStats.Bayes_ZscoreCompare = probs;
        % BUG-6 FIX: was S.groupStatsP_ZscoreCompare (top-level field),
        %             now correctly nested under S.groupStats
        S.groupStats.P_ZscoreCompare     = ps;

        save([saveDir nameOfFile], '-struct', 'S');
    end

    % -----------------------------------------------------------------------
    % 7c – Z-score scatter (two selected stimuli)
    % -----------------------------------------------------------------------

    nexttile

    % Default: compare 1st and 2nd stimulus; override with StimsToCompare if set
    if isempty(params.StimsToCompare)
        ind1 = 1; ind2 = 2;
    else
        ind1 = find(strcmp(Stims2Comp2, params.StimsToCompare{1}));
        ind2 = find(strcmp(Stims2Comp2, params.StimsToCompare{2}));
    end

    ValsToCompare = {StimZS{ind1}, StimZS{ind2}};

    % Only plot if the two vectors are the same length (same neuron set)
    if numel(ValsToCompare{1}) == numel(ValsToCompare{2})
        scatter(ValsToCompare{1}, ValsToCompare{2}, 10, AnIndex, "filled", "MarkerFaceAlpha", 0.5)
        colormap(colormapUsed)
        hold on
        axis equal
        lims = [min(y(y > -inf)), max(y)];
        plot(lims, lims, 'k--', 'LineWidth', 1.5)
        lims = [-5 40];
        ylim(lims); xlim(lims)
        xlabel(Stims2Comp(ind1)); ylabel(Stims2Comp(ind2))
    end

    % -----------------------------------------------------------------------
    % 7d – Spike-rate swarm chart
    % -----------------------------------------------------------------------

    y = cell2mat(stimRSP);   % all spike rates concatenated

    nexttile
    if ~params.EachStimSignif
        swarmchart(x, y, 5, colormapUsed(allColorIndices,:), 'filled', 'MarkerFaceAlpha', 0.7);
    else
        swarmchart(x, y, 5, 'filled', 'MarkerFaceAlpha', 0.7);
    end
    xticks(1:8);
    xticklabels(Stims2Comp2);
    ylabel('Spike Rate');
    set(fig, 'Color', 'w')

    % -----------------------------------------------------------------------
    % 7e – Hierarchical bootstrapping for spike-rate group comparison
    % -----------------------------------------------------------------------

    if params.overwriteGroupStats || ~isfield(S, 'groupStats')
        FirstStim = y(x == 1);
        BootFirst = hierBoot(FirstStim(~isnan(FirstStim)), 10000, ...
                             InsIndex(~isnan(FirstStim)), AnIndex(~isnan(FirstStim)));
        j = 1;
        for i = 2:numel(Stims2Comp2)
            secondaryStim = y(x == i);
            secondaryStim(isnan(secondaryStim)) = 0;
            validMask     = secondaryStim ~= -inf;
            secondaryStim = secondaryStim(validMask);
            BootSec  = hierBoot(secondaryStim, 10000, InsIndex(validMask), AnIndex(validMask));
            probs{j} = get_direct_prob(BootFirst, BootSec);
            ps{j}    = mean(BootSec >= BootFirst);
            j        = j + 1;
        end
        S.groupStats.Bayes_SpikeRateCompare = probs;
        S.groupStats.P_SpikeRateCompare     = ps;
    end

    % -----------------------------------------------------------------------
    % 7f – Spike-rate scatter (same two stimuli as Z-score scatter)
    % -----------------------------------------------------------------------

    nexttile
    ValsToCompare = {stimRSP{ind1}, stimRSP{ind2}};

    if numel(ValsToCompare{1}) == numel(ValsToCompare{2})
        scatter(ValsToCompare{1}, ValsToCompare{2}, 10, AnIndex, "filled", "MarkerFaceAlpha", 0.5)
        colormap(colormapUsed)
        hold on
        axis equal
        lims = [0, max(xlim)];
        plot(lims, lims, 'k--', 'LineWidth', 1.5)
        ylim(lims); xlim(lims)
        xlabel(Stims2Comp(ind1)); ylabel(Stims2Comp(ind2))
    end

end  % end if/else ComparePairs

% =========================================================================
% SECTION 8 – FRACTION-RESPONSIVE ANALYSIS
%             Compares the proportion of neurons responding to each stimulus
%             using simple bootstrapping at the insertion level.
% =========================================================================

% Set default pair for fraction-responsive comparison
if isempty(params.ComparePairs)
    pairs = {Stims2Comp{1}, Stims2Comp{2}};
else
    pairs = params.ComparePairs;
end

% Find insertions with data for both stimuli in the pair
[G, ~] = findgroups(S.TableRespNeurs.insertion);
hasAll  = splitapply(@(s) all(ismember(unique(categorical(pairs)), s)), ...
                     S.TableRespNeurs.stimulus, G);
tempTable = S.TableRespNeurs( ...
    hasAll(G) & ismember(S.TableRespNeurs.stimulus, unique(categorical(pairs))), :);

nBoot = 10000;
j     = 1;
ps    = zeros(1, size(pairs, 1));

% Bootstrap the difference in responsive fraction between the two stimuli
for i = 1:size(pairs, 1)

    diffs = [];

    for ins = unique(S.TableRespNeurs.insertion)'

        idx1 = S.TableRespNeurs.insertion == categorical(ins) & ...
               S.TableRespNeurs.stimulus  == pairs{j,1};
        idx2 = S.TableRespNeurs.insertion == categorical(ins) & ...
               S.TableRespNeurs.stimulus  == pairs{j,2};

        if any(idx1) && any(idx2)
            % Compute difference of fractions (responsive / total)
            % Note: totalSomaticN from idx1 is used as the shared denominator
            f1 = S.TableRespNeurs.respNeur(idx1) / S.TableRespNeurs.totalSomaticN(idx1);
            f2 = S.TableRespNeurs.respNeur(idx2) / S.TableRespNeurs.totalSomaticN(idx1);
            diffs(end+1, 1) = f1 - f2;
        end
    end

    % Simple bootstrap of mean difference (one value per insertion → no hierarchy needed)
    bootDiff = bootstrp(nBoot, @mean, diffs);
    ps(j)    = mean(bootDiff <= 0);   % p-value
    j        = j + 1;
end

% Add column: total responsive neurons per insertion (summed across both stimuli)
[G, ~]                   = findgroups(tempTable.insertion);
totals                   = splitapply(@sum, tempTable.respNeur, G);
tempTable.TotalRespNeur  = totals(G);

% Plot fractions with significance annotation
fig = plotSwarmBootstrapWithComparisons(tempTable, pairs, ps, ...
    {'respNeur','totalSomaticN'}, fraction=true, showBothAndDiff=false,yLegend='Responsive/total units', ...
    diff=false, filled=false, Xjitter='none', Alpha=0.6, drawLines=true);

TotalRespUnits = sum(tempTable.respNeur);

TotalRespUnitsPair1 = sum(tempTable.respNeur(tempTable.stimulus == string(pairs{1})));

TotalRespUnitsPair2 = sum(tempTable.respNeur(tempTable.stimulus == string(pairs{2})));


ax = gca;
ax.YAxis.FontSize = 8; ax.YAxis.FontName = 'helvetica';
ax.XAxis.FontSize = 8; ax.XAxis.FontName = 'helvetica';
set(fig, 'Units', 'centimeters', 'Position', [20 20 5 6]);
ylabel('Responsive/Total responsive')
title('')

% Push axes up slightly to make room for bottom title
pos    = get(gca, 'Position');   % [left bottom width height]
pos(2) = pos(2) + 0.05;          % shift bottom edge up
set(gca, 'Position', pos);

% Horizontal title at the bottom
annotation('textbox', [0.1, 0.01, 0.8, 0.04], ...
    'String',              sprintf('TR = %d - %s = %d - %s = %d',TotalRespUnits,pairs{1},TotalRespUnitsPair1,pairs{2},TotalRespUnitsPair2), ...
    'Rotation',            0, ...
    'EdgeColor',           'none', ...
    'FontSize',            9, ...
    'FontWeight',          'bold', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment',   'middle', ...
    'FitBoxToText',        false);

if params.PaperFig
    vs.printFig(fig, sprintf('ResponsiveUnits-comparison-%s-%s', ...
        params.ComparePairs{1}, params.ComparePairs{2}), PaperFig=params.PaperFig);
end

end  % end function PlotZScoreComparison