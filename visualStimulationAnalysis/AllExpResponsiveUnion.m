function [T, summary] = AllExpResponsiveUnion(expList, params)
% AllExpResponsiveUnion  Count somatic units responsive to AT LEAST ONE visual
%   stimulus, pooling recordings that did not all run the same stimulus set.
%
%   WHY THIS EXISTS (and is not a mode of AllExpAnalysis):
%     AllExpAnalysis skips an entire experiment when ANY stimulus in
%     ComparePairs is absent, and its fraction-responsive section keeps only
%     insertions holding every compared item.  With six stimuli that collapses
%     the dataset to the few recordings that ran all six.  It also applies a
%     single BaseRespWindow to every stimulus, whereas the published calls use
%     250 ms (MB/RG), 500 ms (MBR) and 1000 ms (SDGm/SDGs).
%
%     This function instead visits each experiment once, takes each stimulus
%     family's p-values from THAT FAMILY'S OWN cached StatisticsPerNeuron
%     result (so every family keeps the window and focus it was published
%     with), ORs the per-neuron significance masks across whichever families
%     that recording actually has, and counts.
%
%   HOW STATS ARE READ (do not change without reading this):
%     Cached results are fetched with a BARE call, obj.StatisticsPerNeuron.
%     StatisticsPerNeuron returns the saved struct verbatim only when
%     nargin==1; passing any name-value pair triggers computationParamsMatch,
%     which — against the focused cache AllExpAnalysis wrote — fails, then
%     RECOMPUTES with different parameters and OVERWRITES the cache on disk.
%     That was the findFalseNegAndPos Stage-2 bug.  The only call here that
%     passes parameters is the one that creates a cache which does not exist
%     yet (nothing to clobber in that case).
%
%   INPUTS
%     expList  (1,:) double   Row vector of experiment IDs from master Excel.
%     params   Name-value     See arguments block below.
%
%   OUTPUTS
%     T        table   One row per experiment (insertion numbers preserved):
%                      animal, insertion, realInsertion, recordingName,
%                      nSomatic, nUnion, nUnionFDR, and per stimulus
%                      present_<stim>, nResp_<stim>, nUnique_<stim>.
%     summary  struct  Pooled totals, somatic-per-responsive ratio,
%                      per-stimulus contributions, hierarchical-bootstrap CI
%                      on the union fraction, and the cache-parameter audit.
%
%   EXAMPLES
%     % Cross-check against the Figure-1 MB/RG number (must equal its TR=):
%     T = AllExpResponsiveUnion([1, 8, 40:43,49:54,64:66,68:85 87:97], ...
%             Stimuli = {'MB','RG'}, computeMissing = false);
%
%     % Full six-stimulus union, creating the FFF caches that do not exist yet:
%     [T, S] = AllExpResponsiveUnion([1, 8, 40:43,49:54,64:66,68:85 87:97], ...
%             PaperFig = true);
%
%   KNOWN LIMITATIONS (state these in Methods)
%     - A union of k tests at alpha=0.05 inflates the per-neuron false-positive
%       rate (~26% at k=6).  nUnionFDR is the Benjamini-Hochberg-corrected
%       companion count; nUnion is the uncorrected number that matches the
%       per-stimulus counts in the existing figures.
%     - Each family keeps its own analysis window, so "responsive" is not an
%       equal-sensitivity criterion across stimuli.
%     - Default session only (Session=1, as in AllExpAnalysis Mode 1).  A
%       stimulus recorded only in a second session is counted as absent.
%
%   See also: AllExpAnalysis, findFalseNegAndPos, StatisticsPerNeuron,
%             hierBootMatchFreq, plotSwarmBootstrapWithComparisons

% =========================================================================
%                          ARGUMENTS BLOCK
% =========================================================================
arguments
    expList (1,:) double                                    % row vector of experiment IDs from master Excel

    % --- Stimuli to union ---
    params.Stimuli          cell    = {'MB','RG','MBR','SDGm','SDGs','FFF'}  % stimulus families ORed together

    % --- Responsiveness criterion ---
    params.threshold        double  = 0.05                  % per-stimulus p-value cutoff
    params.useTtest         logical = false                 % true -> pValTTest, false -> pvalsResponse (permutation; primary criterion)

    % --- Missing-cache policy (the only path that ever writes stats) ---
    params.computeMissing   logical = true                  % compute + cache StatisticsPerNeuron when a family has no cache yet
    params.MissingBaseRespWindow struct = struct( ...       % BaseRespWindow used ONLY when creating a missing cache (ms)
        'MB',   250, ...                                    % matches the published MB/RG call
        'RG',   250, ...                                    % matches the published MB/RG call
        'MBR',  500, ...                                    % matches the published MB/MBR call
        'SDGm', 1000, ...                                   % matches the published grating calls
        'SDGs', 1000, ...                                   % matches the published grating calls
        'FFF',  300)                                        % StatisticsPerNeuron default; FFF has no published call yet
    params.RespDurationWin  double  = 100                   % ResponseWindow durationWindow, missing-cache path only

    % --- Inference / output ---
    params.nBoot            double  = 10000                 % hierBootMatchFreq replicates for the union-fraction CI
    params.plotFig          logical = true                  % draw the per-insertion fraction swarm
    params.PaperFig         logical = false                 % export the figure via printFig
    params.overwrite        logical = false                 % force rerun of the per-experiment loop (ignore pooled cache)
end

% =========================================================================
% SECTION 1 — VALIDATE THE REQUESTED STIMULUS LIST
% =========================================================================

stims  = params.Stimuli(:)';                                % row cell of stimulus abbreviations
nStim  = numel(stims);                                      % number of families in the union
assert(nStim >= 1, 'Stimuli must contain at least one stimulus abbreviation.');
assert(numel(unique(stims)) == nStim, 'Stimuli contains duplicates: %s', strjoin(stims, ', '));

% Every requested family must be constructible and must have a missing-cache window
known = {'MB','RG','MBR','SDGm','SDGs','FFF','NI','NV'};    % families buildObj knows how to construct
for si = 1:nStim                                            % check each requested family up front
    assert(ismember(stims{si}, known), ...                  % unknown abbreviation -> fail now, not mid-loop
        'Unknown stimulus "%s". Known: %s', stims{si}, strjoin(known, ', '));
    assert(isfield(params.MissingBaseRespWindow, stims{si}), ...   % needed only for the missing-cache path
        'MissingBaseRespWindow has no field "%s" — add one before requesting that stimulus.', stims{si});
end

% Sanitised names for use as table variable suffixes (identical here, but explicit)
stimVar = matlab.lang.makeValidName(stims);                 % e.g. 'SDGm' -> 'SDGm' (already valid)

% =========================================================================
% SECTION 2 — DIRECTORY SETUP AND POOLED-CACHE MANAGEMENT
% =========================================================================

NP0      = loadNPclassFromTable(expList(end));              % load one recording just to resolve the data root
rootPath = extractBefore(NP0.recordingDir, 'lizards');      % path up to (excluding) the 'lizards' folder
rootPath = [rootPath 'lizards'];                            % include the 'lizards' folder itself
saveDir  = fullfile(rootPath, 'Combined_lizard_analysis');  % pooled-analysis directory (same as AllExpAnalysis)
if ~exist(saveDir, 'dir')                                   % create it on first use
    mkdir(saveDir);
end

% Filename encodes the experiment range and the stimulus set so different
% unions never collide on disk
nameOfFile = sprintf('Ex_%d-%d_ResponsiveUnion_%s.mat', ...
    expList(1), expList(end), strjoin(stims, '-'));
savePath   = fullfile(saveDir, nameOfFile);                 % full path of the pooled cache

% Skip the per-experiment loop only if a cache exists for exactly this
% experiment list, this stimulus set and the same responsiveness criterion
runLoop = true;                                             % assume we must compute
if exist(savePath, 'file') == 2 && ~params.overwrite
    S = load(savePath);                                     % peek at the cached pooled result
    sameCrit = isfield(S, 'threshold') && isequal(S.threshold, params.threshold) && ...
               isfield(S, 'useTtest')  && isequal(S.useTtest,  params.useTtest);
    if isfield(S, 'expList') && isequal(S.expList, expList) && sameCrit
        runLoop = false;                                    % cache is valid for this call
        fprintf('Loaded pooled union from cache: %s\n', savePath);
    else
        fprintf('Pooled cache exists but expList/criterion differ — recomputing.\n');
    end
end

% =========================================================================
% SECTION 3 — PER-EXPERIMENT LOOP
% =========================================================================

if runLoop

    % Map-based indices guarantee that the same animal / (animal,insertion)
    % always receives the same number regardless of expList ordering, and
    % that insertion indices are globally unique across animals — required
    % by hierBootMatchFreq.
    animalMap    = containers.Map('KeyType','char','ValueType','double');   % animal name -> index
    insertionMap = containers.Map('KeyType','char','ValueType','double');   % 'PV103__Ins7' -> index
    nextAnimalIdx    = 0;                                   % running animal counter
    nextInsertionIdx = 0;                                   % running insertion counter

    rows      = {};                                         % per-experiment result structs (converted to T at the end)
    auditRows = {};                                         % per (experiment,stimulus) cache-parameter audit entries
    insUsage  = containers.Map('KeyType','double','ValueType','any');       % insertion index -> recording names using it

    % Configuration template per family: the params of the FIRST cache seen
    % for that stimulus.  A family whose cache is missing in some recordings
    % is then computed with the SAME configuration as its existing caches,
    % instead of falling back to defaults and silently mixing analysis
    % windows within one stimulus.
    templateParams = containers.Map('KeyType','char','ValueType','any');    % stimulus -> saved params struct

    % ---- 3-pre: discover the configuration template of each family ----
    % Filling templates DURING the main loop would make the result depend on
    % expList order: a recording without a cache placed before every cached
    % one would be computed with fallback defaults, permanently caching a
    % configuration that differs from the rest of its family.  This pre-scan
    % removes that ordering dependency.  It only peeks at the saved params
    % (a plain load of the .mat file) — never a bare StatisticsPerNeuron
    % call — so no cache is read into memory, recomputed or rewritten here.
    if params.computeMissing
        fprintf('Scanning for per-family configuration templates (skip with computeMissing=false)...\n');
        for ex = expList
            if templateParams.Count == nStim, break; end     % every requested family resolved
            NPt = loadNPclassFromTable(ex);                  % recording to probe
            for si = 1:nStim
                sn = stims{si};                              % family still needing a template
                if templateParams.isKey(sn), continue; end
                oT = buildObj(NPt, sn);                      % [] when the stimulus is absent here
                if isempty(oT), continue; end
                fT = fullfile(oT.visualStimAnalysisFolder, 'StatisticsPerNeuron.mat');
                if isfile(fT)
                    cT = load(fT, 'params');                 % read the stored parameters only
                    if isfield(cT, 'params')
                        templateParams(sn) = cT.params;      % this family's existing configuration
                        fprintf('  template for %s taken from exp %d.\n', sn, ex);
                    end
                end
            end
        end
        for si = 1:nStim                                     % families with no cache anywhere
            if ~templateParams.isKey(stims{si})
                fprintf('  no existing cache for %s — missing ones will use MissingBaseRespWindow=%g ms.\n', ...
                    stims{si}, params.MissingBaseRespWindow.(stims{si}));
            end
        end
    end

    for ex = expList

        % ---- 3a: Load the recording ----
        NP = loadNPclassFromTable(ex);                      % Neuropixels recording object
        fprintf('\n=== Exp %d: %s ===\n', ex, NP.recordingName);

        % ---- 3b: Metadata (animal, insertion, real insertion) ----
        % Identical parsing to AllExpAnalysis so insertion numbers match
        % between this function and the swarm/raster figures.
        animalID = string(regexp(NP.recordingName, 'PV\d+', 'match', 'once'));    % Pogona ID
        if animalID == ""                                                          % not a PV recording
            animalID = string(regexp(NP.recordingName, 'SA\d+', 'match', 'once')); % Laudakia ID
        end
        assert(animalID ~= "", 'Could not parse an animal ID (PV\\d+ / SA\\d+) from %s', NP.recordingName);

        insStr = regexp(NP.recordingDir, 'Insertion\d+', 'match', 'once');   % e.g. 'Insertion7'
        insNum = str2double(regexp(insStr, '\d+', 'match'));                 % numeric insertion from the path
        assert(~isnan(insNum), 'Could not parse Insertion\\d+ from %s', NP.recordingDir);

        nameTok         = strsplit(NP.recordingName, '_');   % {'PV103','Experiment','12','6','23','1'}
        realInsertionID = str2double(nameTok{end});          % trailing token = human-facing insertion number
        assert(~isnan(realInsertionID), ...
            'recordingName does not end in a numeric insertion token: %s', NP.recordingName);

        animalKey = char(animalID);                          % map key for the animal
        if ~animalMap.isKey(animalKey)                       % first time we see this animal
            nextAnimalIdx      = nextAnimalIdx + 1;          % allocate a new index
            animalMap(animalKey) = nextAnimalIdx;
        end
        animalIdx = animalMap(animalKey);                    % stable animal index

        insKey = sprintf('%s__Ins%d', animalKey, insNum);    % globally unique insertion key
        if ~insertionMap.isKey(insKey)                       % first time we see this insertion
            nextInsertionIdx        = nextInsertionIdx + 1;  % allocate a new index
            insertionMap(insKey)    = nextInsertionIdx;
        end
        insertionIdx = insertionMap(insKey);                 % stable, globally unique insertion index

        % Track which recordings share an insertion index (warned about below)
        if insUsage.isKey(insertionIdx)                      % this insertion already contributed a recording
            insUsage(insertionIdx) = [insUsage(insertionIdx), {NP.recordingName}];
        else
            insUsage(insertionIdx) = {NP.recordingName};     % first recording for this insertion
        end

        % ---- 3c: Build one analysis object per unique class key ----
        % SDGm and SDGs share a single StaticDriftingGratingAnalysis object
        % (and therefore a single StatisticsPerNeuron.mat holding both
        % .Static and .Moving), so build/load once and extract twice.
        objByKey = containers.Map();                         % class key -> analysis object
        for si = 1:nStim
            k = getObjKey(stims{si});                        % 'SDGm'/'SDGs' -> 'SDG'
            if ~objByKey.isKey(k)                            % build each class at most once
                objByKey(k) = buildObj(NP, stims{si});       % [] when the stimulus was not recorded
            end
        end

        % ---- 3d: Somatic unit universe for this recording ----
        % 'good' units ARE the somatic set (StatisticsPerNeuron keeps
        % p.ic(:, label=='good')), and every stats vector is in that order.
        refObj = [];                                         % first object that constructed successfully
        objKeys = keys(objByKey);                            % all class keys attempted
        for k = 1:numel(objKeys)
            if ~isempty(objByKey(objKeys{k}))                % found a usable object
                refObj = objByKey(objKeys{k});
                break
            end
        end
        if isempty(refObj)                                   % no requested stimulus exists in this recording
            warning('AllExpResponsiveUnion:noStimuli', ...
                'Exp %d (%s): none of the requested stimuli are present — experiment skipped.', ...
                ex, NP.recordingName);
            continue
        end

        ps       = NP.convertPhySorting2tIc(refObj.spikeSortingFolder);        % phy sorting for this recording
        phy_IDg  = ps.phy_ID(string(ps.label') == 'good');                     % phy IDs of somatic ('good') units
        nSomatic = numel(phy_IDg);                                             % denominator for this recording
        assert(nSomatic > 0, 'Exp %d: no somatic ("good") units found.', ex);
        fprintf('  %d somatic units.\n', nSomatic);

        % Every family must come from the SAME spike sorting, otherwise the
        % per-neuron masks refer to different unit sets and cannot be ORed.
        for k = 1:numel(objKeys)
            o = objByKey(objKeys{k});                        % candidate object
            if isempty(o), continue; end                     % stimulus absent — nothing to check
            assert(strcmp(o.spikeSortingFolder, refObj.spikeSortingFolder), ...
                ['Exp %d: %s uses a different spike sorting folder than %s. ' ...
                 'Per-neuron masks cannot be unioned across different sortings.'], ...
                ex, objKeys{k}, class(refObj));
        end

        % ---- 3e: Per-family p-values ----
        P       = nan(nSomatic, nStim);                      % p-value per (neuron, family); NaN = family absent
        present = false(1, nStim);                           % did this recording yield usable stats for the family?

        for si = 1:nStim
            sn  = stims{si};                                 % stimulus abbreviation
            obj = objByKey(getObjKey(sn));                   % shared analysis object for its class

            if isempty(obj) || isempty(obj.VST)              % stimulus was never recorded here
                fprintf('  [%s] not present.\n', sn);
                continue
            end

            % Cache path must be built by hand: getAnalysisFileName derives
            % the filename from the CALLING method via dbstack, so calling it
            % from here would resolve to 'AllExpResponsiveUnion.mat'.
            statFile = fullfile(obj.visualStimAnalysisFolder, 'StatisticsPerNeuron.mat');

            Stats = [];                                      % stats struct for this family
            src   = '';                                      % provenance string for the audit
            if isfile(statFile)
                % BARE call: returns the saved struct verbatim (nargin==1
                % passthrough).  Passing parameters here would recompute and
                % overwrite the cache the published figures were built from.
                try
                    Stats = obj.StatisticsPerNeuron;         % load cache exactly as written
                    src   = 'cache';
                    if isfield(Stats, 'params') && ~templateParams.isKey(sn)
                        templateParams(sn) = Stats.params;   % first cache of this family becomes the template
                    end
                catch ME
                    warning('AllExpResponsiveUnion:statsLoadFailed', ...
                        'Exp %d [%s]: could not load cached stats (%s) — family skipped.', ex, sn, ME.message);
                    continue
                end
            elseif params.computeMissing
                % No cache exists, so nothing can be clobbered.  Prefer the
                % configuration this family already uses in other recordings
                % (templateParams) so one stimulus never mixes windows; fall
                % back to MissingBaseRespWindow only when no cache of this
                % family has been seen yet.
                if templateParams.isKey(sn)
                    statArgs = paramsToArgs(templateParams(sn));             % replicate the family's existing configuration
                    srcNote  = 'template from an existing cache of this family';
                else
                    statArgs = {'BaseRespWindow', params.MissingBaseRespWindow.(sn)};   % first-ever cache for this family
                    srcNote  = sprintf('MissingBaseRespWindow=%g ms', params.MissingBaseRespWindow.(sn));
                end
                fprintf('  [%s] no cache — computing (%s).\n', sn, srcNote);
                try
                    obj.ResponseWindow('overwrite', false, ...               % response traces + condition matrix
                        'durationWindow', params.RespDurationWin);
                    Stats = obj.StatisticsPerNeuron(statArgs{:});            % explicit params: safe, no prior cache
                    src   = 'computed';
                    if ~templateParams.isKey(sn) && isfield(Stats, 'params')
                        templateParams(sn) = Stats.params;                   % this becomes the family template
                    end
                    warning('AllExpResponsiveUnion:statsComputed', ...
                        'Exp %d [%s]: stats did not exist and were computed fresh (%s).', ex, sn, srcNote);
                catch ME
                    % A template taken from an older cache may name a
                    % parameter the current StatisticsPerNeuron no longer
                    % accepts.  Retry once with the plain fallback so one
                    % stale field does not silently drop the recording; the
                    % audit will then show the configuration difference.
                    retried = false;                                          % did the fallback attempt run?
                    if templateParams.isKey(sn)
                        warning('AllExpResponsiveUnion:templateArgsFailed', ...
                            'Exp %d [%s]: computing with the family template failed (%s) — retrying with BaseRespWindow=%g ms.', ...
                            ex, sn, ME.message, params.MissingBaseRespWindow.(sn));
                        try
                            Stats   = obj.StatisticsPerNeuron('BaseRespWindow', params.MissingBaseRespWindow.(sn));
                            src     = 'computed';
                            srcNote = sprintf('MissingBaseRespWindow=%g ms (template rejected)', params.MissingBaseRespWindow.(sn));
                            retried = true;
                            warning('AllExpResponsiveUnion:statsComputed', ...
                                'Exp %d [%s]: stats did not exist and were computed fresh (%s).', ex, sn, srcNote);
                        catch ME2
                            ME = ME2;                                         % report the fallback's failure
                        end
                    end
                    if ~retried
                        warning('AllExpResponsiveUnion:statsComputeFailed', ...
                            'Exp %d [%s]: computing stats failed (%s) — family skipped.', ex, sn, ME.message);
                        continue
                    end
                end
            else
                fprintf('  [%s] no stats cache and computeMissing=false — family skipped.\n', sn);
                continue
            end

            % Pull the per-neuron p-values for this family.  Returns empty
            % when the required sub-struct is missing (e.g. a grating session
            % with no moving block), which counts as "family absent".
            pv = extractP(Stats, sn, params.useTtest);
            if isempty(pv)
                warning('AllExpResponsiveUnion:noPvals', ...
                    'Exp %d [%s]: stats exist but hold no p-values for this family — skipped.', ex, sn);
                continue
            end

            % Alignment invariant: stats vectors are indexed by good unit.
            assert(numel(pv) == nSomatic, ...
                ['Exp %d [%s]: %d p-values but %d somatic units. ' ...
                 'The stats cache is stale relative to the current spike sorting — recompute it.'], ...
                ex, sn, numel(pv), nSomatic);

            nNaN = sum(isnan(pv));                           % NaN p-values never count as responsive
            if nNaN > 0
                warning('AllExpResponsiveUnion:nanPvals', ...
                    'Exp %d [%s]: %d/%d p-values are NaN and are treated as non-responsive.', ...
                    ex, sn, nNaN, nSomatic);
            end

            P(:, si)    = pv(:);                             % store this family's p-values
            present(si) = true;                              % family contributed
            fprintf('  [%s] %d/%d responsive (p<%.3g, %s).\n', ...
                sn, sum(pv < params.threshold), nSomatic, params.threshold, src);

            % Audit entry: which parameters produced the numbers just used
            if isfield(Stats, 'params')                      % every save writes S.params
                auditRows{end+1} = struct( ...               %#ok<AGROW>
                    'exp', ex, 'stimulus', string(sn), 'source', string(src), ...
                    'signature', string(paramSignature(Stats.params)));
            else
                auditRows{end+1} = struct( ...               %#ok<AGROW> pre-params cache
                    'exp', ex, 'stimulus', string(sn), 'source', string(src), ...
                    'signature', "<no params field in cache>");
            end
        end

        if ~any(present)                                     % nothing usable in this recording
            warning('AllExpResponsiveUnion:noUsableStats', ...
                'Exp %d (%s): no usable stats for any requested stimulus — experiment skipped.', ...
                ex, NP.recordingName);
            continue
        end

        % ---- 3f: Union masks ----
        sigMat    = P < params.threshold;                    % [nSomatic x nStim]; NaN comparisons are false
        maskUnion = any(sigMat, 2);                          % responsive to at least one stimulus (uncorrected)
        maskFDR   = bhAnyRejection(P, params.threshold);     % same, after BH across the families present per neuron

        % Marginal contribution: units flagged by exactly one family
        onlyOne   = sum(sigMat, 2) == 1;                     % neurons with a single significant family
        nUnique   = zeros(1, nStim);                         % per-family unique contribution
        for si = 1:nStim
            nUnique(si) = sum(sigMat(:, si) & onlyOne);      % this family was that neuron's only hit
        end

        % ---- 3g: Assemble the row ----
        r = struct();                                        % one struct per experiment
        r.exp             = ex;                              % experiment ID
        r.animal          = string(animalID);                % animal name
        r.animalIdx       = animalIdx;                       % numeric animal index (bootstrap level 2)
        r.insertion       = insertionIdx;                    % running insertion index (globally unique)
        r.realInsertion   = realInsertionID;                 % human-facing insertion number
        r.recordingName   = string(NP.recordingName);        % full recording name
        r.nSomatic        = nSomatic;                        % somatic units in this recording
        r.nUnion          = sum(maskUnion);                  % responsive to >=1 stimulus (uncorrected)
        r.nUnionFDR       = sum(maskFDR);                    % responsive to >=1 stimulus (BH across stimuli)
        r.nStimPresent    = sum(present);                    % how many families this recording contributed
        for si = 1:nStim                                     % per-family columns
            r.(['present_' stimVar{si}]) = present(si);      % family available?
            r.(['nResp_'   stimVar{si}]) = sum(sigMat(:, si));   % responsive to this family
            r.(['nUnique_' stimVar{si}]) = nUnique(si);          % responsive to ONLY this family
        end
        rows{end+1} = r;                                     %#ok<AGROW>

        fprintf('  -> union %d / %d somatic (%.1f%%), FDR-corrected %d.\n', ...
            r.nUnion, nSomatic, 100*r.nUnion/nSomatic, r.nUnionFDR);
    end

    assert(~isempty(rows), 'No experiment produced usable data — nothing to summarise.');
    T = struct2table([rows{:}]);                             % wide per-experiment table

    % Loud warning if two recordings map to the same insertion index: their
    % somatic units are counted twice in the pooled totals, and their phy IDs
    % come from different sortings so they cannot be merged.
    dupKeys = keys(insUsage);                                % all insertion indices seen
    for k = 1:numel(dupKeys)
        names = insUsage(dupKeys{k});                        % recordings sharing this index
        if numel(names) > 1
            warning('AllExpResponsiveUnion:duplicateInsertion', ...
                ['Insertion index %d is used by %d recordings (%s). Their somatic units are ' ...
                 'summed independently in the pooled totals — verify this is intended.'], ...
                dupKeys{k}, numel(names), strjoin(names, ', '));
        end
    end

    % Cache the pooled result together with the criterion it was built under
    S = struct();                                            % save struct
    S.expList   = expList;                                   % experiment list this result covers
    S.Stimuli   = stims;                                     % stimulus set (cell of abbreviations)
    S.T         = T;                                         % per-experiment table
    S.audit     = struct2table([auditRows{:}]);              % cache-parameter audit
    S.threshold = params.threshold;                          % responsiveness cutoff
    S.useTtest  = params.useTtest;                           % which p-value field was used
    S.params    = params;                                    % full parameter set for reproducibility
    save(savePath, '-struct', 'S');                          % write pooled cache
    fprintf('\nSaved pooled union to %s\n', savePath);
end

% Restore from cache when the loop was skipped
T     = S.T;                                                 % per-experiment table
audit = S.audit;                                             % cache-parameter audit

% =========================================================================
% SECTION 4 — PARAMETER AUDIT
% =========================================================================
% Each StatisticsPerNeuron cache stores the parameters that produced it.
% Different windows ACROSS families are intended; more than one signature
% WITHIN a family means the pooled count mixes analysis windows.

fprintf('\n--- Cache parameter audit ---\n');
auditStims = unique(audit.stimulus, 'stable');               % families that contributed
for si = 1:numel(auditStims)
    sn   = auditStims(si);                                   % this family
    sigs = audit.signature(audit.stimulus == sn);            % its parameter signatures
    [uSig, ~, ic] = unique(sigs, 'stable');                  % distinct signatures
    cnt  = accumarray(ic, 1);                                % experiments per signature
    fprintf('  %s (%d experiments):\n', sn, numel(sigs));
    for k = 1:numel(uSig)
        fprintf('     %3d x  %s\n', cnt(k), uSig(k));        % how many experiments used each signature
    end
    if numel(uSig) > 1                                       % inconsistent within a family — real problem
        warning('AllExpResponsiveUnion:mixedStatsParams', ...
            ['Stimulus %s was analysed with %d DIFFERENT parameter sets across experiments. ' ...
             'The pooled count mixes analysis windows — rerun AllExpAnalysis for this family ' ...
             'so every experiment shares one configuration.'], sn, numel(uSig));
    end
end

% =========================================================================
% SECTION 5 — POOLED SUMMARY
% =========================================================================

summary = struct();                                          % output summary
summary.nRecordings   = height(T);                           % recordings contributing
summary.nInsertions   = numel(unique(T.insertion));          % distinct insertions
summary.nAnimals      = numel(unique(T.animal));             % distinct animals
summary.totalSomatic  = sum(T.nSomatic);                     % pooled somatic units
summary.totalUnion    = sum(T.nUnion);                       % pooled units responsive to >=1 stimulus
summary.totalUnionFDR = sum(T.nUnionFDR);                    % same, BH-corrected across stimuli
summary.unionFraction = summary.totalUnion / summary.totalSomatic;              % pooled fraction
summary.somaticPerResponsive = summary.totalSomatic / summary.totalUnion;       % somatic units per responsive unit

% Per-family totals, unique contributions, and how many recordings had them
summary.perStimulus = table('Size', [nStim 4], ...
    'VariableTypes', {'string','double','double','double'}, ...
    'VariableNames', {'stimulus','nRecordings','nResponsive','nUniqueContribution'});
for si = 1:nStim
    summary.perStimulus.stimulus(si)            = string(stims{si});                 % family name
    summary.perStimulus.nRecordings(si)         = sum(T.(['present_' stimVar{si}])); % recordings with this family
    summary.perStimulus.nResponsive(si)         = sum(T.(['nResp_'   stimVar{si}])); % units it flagged
    summary.perStimulus.nUniqueContribution(si) = sum(T.(['nUnique_' stimVar{si}])); % units only it flagged
end

% Hierarchical bootstrap CI on the union fraction: resample animals, then
% insertions (never plain hierBoot, whose top level resamples observations).
fracPerRec = T.nUnion ./ T.nSomatic;                         % per-recording union fraction
animalIdxV = double(T.animalIdx);                            % animal level
insIdxV    = double(T.insertion);                            % insertion level (globally unique)
bootFrac   = hierBootMatchFreq(fracPerRec, params.nBoot, animalIdxV, insIdxV);   % bootstrap distribution of the mean
summary.unionFractionMean = mean(fracPerRec);                % unweighted mean across recordings
summary.ciUnionFraction   = prctile(bootFrac, [2.5 97.5]);   % 95% percentile CI
summary.audit             = audit;                           % carry the audit with the summary

% ---- Printed report ----
fprintf('\n--- Responsive-unit union (%s) ---\n', strjoin(stims, ', '));
fprintf('  Recordings / insertions / animals : %d / %d / %d\n', ...
    summary.nRecordings, summary.nInsertions, summary.nAnimals);
fprintf('  Somatic units                     : %d\n', summary.totalSomatic);
fprintf('  Responsive to >=1 stimulus        : %d (%.1f%% pooled)\n', ...
    summary.totalUnion, 100*summary.unionFraction);
fprintf('  Same, BH-corrected across stimuli : %d (%.1f%% pooled)\n', ...
    summary.totalUnionFDR, 100*summary.totalUnionFDR/summary.totalSomatic);
fprintf('  Somatic units per responsive unit : %.2f\n', summary.somaticPerResponsive);
fprintf('  Mean per-recording fraction       : %.3f  [95%% CI %.3f - %.3f]\n', ...
    summary.unionFractionMean, summary.ciUnionFraction(1), summary.ciUnionFraction(2));
fprintf('  Per stimulus (recordings | responsive | unique):\n');
for si = 1:nStim
    fprintf('     %-5s %4d | %5d | %5d\n', summary.perStimulus.stimulus(si), ...
        summary.perStimulus.nRecordings(si), summary.perStimulus.nResponsive(si), ...
        summary.perStimulus.nUniqueContribution(si));
end

% ---- Sanity invariants (cheap, catch extraction bugs immediately) ----
respCols = T{:, cellfun(@(s) ['nResp_' s], stimVar, 'UniformOutput', false)};   % per-family counts
assert(all(T.nUnion >= max(respCols, [], 2)), 'Union count is below a single-stimulus count — extraction bug.');
assert(all(T.nUnion <= T.nSomatic),          'Union count exceeds the somatic unit count — alignment bug.');
assert(all(T.nUnionFDR <= T.nUnion),         'FDR-corrected count exceeds the uncorrected count — BH bug.');

% =========================================================================
% SECTION 6 — PER-INSERTION FRACTION SWARM
% =========================================================================

if params.plotFig

    % Long format for plotSwarmBootstrapWithComparisons: one row per
    % (recording, column), where columns are the individual stimuli plus the
    % union.  Rows exist only where the stimulus was actually recorded —
    % that heterogeneity is the point of this figure.
    % Column order is set by ROW order, not by reordercats: the plotting
    % function calls reorderStimulusByLevel, which — with no numeric level
    % tokens in these labels — reorders categories to their first appearance
    % among the rows and would discard any category order set here.  So each
    % row carries a sort key (0 = Union, then the requested Stimuli order)
    % and the rows are sorted by it, which puts Union in the first column.
    longRows = {};                                           % accumulator
    for r = 1:height(T)
        longRows{end+1} = {T.animal(r), T.insertion(r), 'Union', ...      %#ok<AGROW>
            T.nUnion(r), T.nSomatic(r), 0};                               % union column, sort key 0
        for si = 1:nStim
            if ~T.(['present_' stimVar{si}])(r), continue; end            % stimulus absent in this recording
            longRows{end+1} = {T.animal(r), T.insertion(r), stims{si}, ...  %#ok<AGROW>
                T.(['nResp_' stimVar{si}])(r), T.nSomatic(r), si};         % sort key = position in Stimuli
        end
    end
    L = vertcat(longRows{:});                                % cell matrix [nRows x 6]
    tblLong = table( ...
        categorical(cellstr(string([L{:,1}]'))), ...         % animal   (categorical, string-based)
        categorical(cellstr(string([L{:,2}]'))), ...         % insertion(categorical, string-based)
        categorical(string(L(:,3))), ...                     % stimulus (categorical)
        cell2mat(L(:,4)), cell2mat(L(:,5)), cell2mat(L(:,6)), ...   % responsive, somatic, column sort key
        'VariableNames', {'animal','insertion','stimulus','respNeur','totalSomaticN','colOrder'});
    tblLong = sortrows(tblLong, 'colOrder');                 % first appearance = Union, then Stimuli order
    tblLong.colOrder = [];                                   % drop the helper column before plotting

    yMax = min(1, 1.15 * max(tblLong.respNeur ./ tblLong.totalSomaticN)); % headroom without wasting axis

    % Empty pairs => no significance brackets are drawn (the plotting
    % function guards on ~isempty(pairs)); showBothAndDiff must be off since
    % there is no pair to difference.
    % Xjitter='density' spreads each column by local point density, so the
    % shape of the column reads as the distribution (Xjitter='none' would
    % stack every recording on one vertical line and hide it).  drawLines is
    % omitted: it only applies to <=2 columns, and pairing lines across
    % columns with different memberships would be misleading anyway.
    % 4th output is the stats block behind the panel; previously discarded,
    % captured so it can be persisted next to the figure for the supplementary
    % table. With no pairs requested, only stats.groups is populated.
    [fig, ~, ~, stUnion] = plotSwarmBootstrapWithComparisons( ...
        tblLong, {}, [], {'respNeur','totalSomaticN'}, ...
        fraction = true, showBothAndDiff = false, diff = false, ...
        yLegend = 'Responsive / somatic units', yMaxVis = yMax, ...
        filled = true, Xjitter = 'density', swarmWidth = 0.7, ...
        dotSize = 14, Alpha = 0.55, ...
        nBoot = params.nBoot);

    summary.stats = stUnion;                                 % also returned for interactive use

    figure(fig);                                             % bring to front
    set(fig, 'Units', 'centimeters', 'Position', [20 20 8 6]);   % publication size
    set(gca, 'Clipping', 'off');                             % keep CI bars visible at the axes edge

    % Bottom annotation with the pooled counts behind the figure
    annotation('textbox', [0.1, 0.01, 0.8, 0.04], ...
        'String', sprintf('TS=%d - TR=%d - ratio=%.2f', ...
            summary.totalSomatic, summary.totalUnion, summary.somaticPerResponsive), ...
        'EdgeColor', 'none', 'FontSize', 5, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FitBoxToText', false);

    if params.PaperFig
        % printFig is a method of the analysis object, so build one for the
        % last experiment purely to resolve the paper-figure directory.
        vsPrint = firstAvailableObj(loadNPclassFromTable(expList(end)), stims);
        if isempty(vsPrint)
            warning('AllExpResponsiveUnion:noPrintObj', ...
                'Could not build an analysis object for printFig — figure not exported.');
        else
            [~, figTag, ~] = fileparts(nameOfFile);          % reuse the cache name as the figure tag
            statsName = sprintf('ResponsiveUnion-%s', figTag);
            vsPrint.printFig(fig, statsName, PaperFig = params.PaperFig);

            % Persist the numbers next to the figure so the supplementary
            % table can pick them up once both are filed into the figure folder.
            saveStatsStruct(vsPrint, statsName, ...
                {'Responsive fraction per stimulus (union set)', stUnion, 'insertions'}, ...
                params, params.PaperFig);
        end
    end
end

end % end function AllExpResponsiveUnion


% #########################################################################
%                        LOCAL HELPER FUNCTIONS
% #########################################################################

function key = getObjKey(stimName)
% getObjKey  Map a stimulus abbreviation to its shared analysis-object key.
%   SDGm and SDGs are two conditions of ONE StaticDriftingGratingAnalysis
%   object and therefore share a single StatisticsPerNeuron cache.
    switch stimName
        case {'SDGm','SDGs'},  key = 'SDG';                  % gratings share one object
        otherwise,             key = stimName;               % everything else is its own class
    end
end


function obj = buildObj(NP, stim)
% buildObj  Construct the analysis object for a stimulus abbreviation.
%   Returns [] when the stimulus was not recorded (constructor throws or the
%   visual-stimulation struct is empty).  Mirrors findFalseNegAndPos/buildObj.
    obj = [];                                                % default: stimulus absent
    try
        switch stim
            case 'MB',            obj = linearlyMovingBallAnalysis(NP);      % moving ball
            case 'RG',            obj = rectGridAnalysis(NP);                % rectangular grid (static ball)
            case 'MBR',           obj = linearlyMovingBarAnalysis(NP);       % moving bar
            case {'SDGm','SDGs'}, obj = StaticDriftingGratingAnalysis(NP);   % static + drifting gratings
            case 'NI',            obj = imageAnalysis(NP);                   % natural images
            case 'NV',            obj = movieAnalysis(NP);                   % natural movies
            case 'FFF',           obj = fullFieldFlashAnalysis(NP);          % full-field flash
            otherwise, error('Unknown stimulus "%s".', stim);
        end
        if isempty(obj.VST)                                  % object built but no stimulus data
            obj = [];
        end
    catch
        obj = [];                                            % constructor failed -> treat as absent
    end
end


function p = extractP(Stats, stim, useTtest)
% extractP  Per somatic-unit p-value for one stimulus family.
%   MB/MBR: minimum p across speeds (matches AllExpAnalysis/extractStimData
%   and findFalseNegAndPos/extractPZ).  SDGs/SDGm read the .Static/.Moving
%   sub-structs of the shared grating cache.  Returns [] when the required
%   field is absent, which the caller treats as "family not available".

    pField = 'pvalsResponse';                                % permutation-test p (primary criterion)
    if useTtest, pField = 'pValTTest'; end                    % t-test p (must match any figure it is compared to)
    p = [];                                                   % default: nothing extractable

    switch stim
        case {'MB','MBR'}
            spF = fieldnames(Stats);                          % all fields of the stats struct
            spF = spF(startsWith(string(spF), "Speed"));      % per-speed sub-structs
            if isempty(spF)                                   % single-condition fallback
                if isfield(Stats, pField), p = Stats.(pField)(:); end
                return
            end
            if ~isfield(Stats.(spF{1}), pField), return; end   % requested p-field not stored
            Pm = nan(numel(Stats.(spF{1}).(pField)), numel(spF));   % [nUnits x nSpeeds]
            for k = 1:numel(spF)
                Pm(:,k) = Stats.(spF{k}).(pField)(:);         % p-values for this speed
            end
            p = min(Pm, [], 2);                               % best (lowest-p) speed per neuron

        case 'SDGs'
            if isfield(Stats, 'Static') && isfield(Stats.Static, pField)
                p = Stats.Static.(pField)(:);                 % static grating block
            end

        case 'SDGm'
            if isfield(Stats, 'Moving') && isfield(Stats.Moving, pField)
                p = Stats.Moving.(pField)(:);                 % drifting grating block
            end

        otherwise
            if isfield(Stats, pField)
                p = Stats.(pField)(:);                        % flat struct (RG, FFF, NI, NV)
            end
    end
end


function mask = bhAnyRejection(P, alpha)
% bhAnyRejection  Per-neuron Benjamini-Hochberg across the stimuli tested.
%   P is [nNeurons x nStim] with NaN where a stimulus was not tested for that
%   recording.  A neuron is called responsive when BH rejects at least one of
%   its own tests at level alpha, i.e. when any sorted p_(i) <= alpha*i/k.
%   Correcting per neuron (rather than across neurons) is what addresses the
%   inflation caused by ORing k stimulus tests together.

    nNeurons = size(P, 1);                                    % neurons in this recording
    mask     = false(nNeurons, 1);                            % default: not responsive
    for n = 1:nNeurons
        pv = P(n, :);                                         % this neuron's p-values
        pv = sort(pv(~isnan(pv)));                            % drop untested stimuli, sort ascending
        k  = numel(pv);                                       % number of tests for this neuron
        if k == 0, continue; end                              % nothing tested
        mask(n) = any(pv(:)' <= alpha * (1:k) / k);           % BH step-up rejection rule
    end
end


function args = paramsToArgs(p)
% paramsToArgs  Turn a saved StatisticsPerNeuron params struct back into a
%   name-value argument list, so a missing cache can be computed with the
%   SAME configuration as the other recordings of that stimulus family.
%   Bookkeeping fields are dropped: 'overwrite' would change behaviour and
%   'inputParams'/'analysisTime' are not computation parameters (the same
%   exclusions computationParamsMatch applies).

    skip   = {'overwrite', 'inputParams', 'analysisTime'};   % non-computational bookkeeping fields
    fn     = fieldnames(p);                                  % every stored parameter
    fn     = fn(~ismember(fn, skip));                        % keep the computational ones
    args   = cell(1, 2*numel(fn));                           % preallocate the name-value list
    for k = 1:numel(fn)
        args{2*k-1} = fn{k};                                 % parameter name
        args{2*k}   = p.(fn{k});                             % parameter value
    end
end


function sig = paramSignature(p)
% paramSignature  Compact, comparable string of the computation parameters
%   stored alongside a StatisticsPerNeuron cache.  Used to detect whether all
%   experiments of one stimulus family were analysed identically.

    f = @(name, def) getFieldOr(p, name, def);                % shorthand with a default
    sig = sprintf(['BaseRespWindow=%s maxCategory=%s SpatialGridMode=%s ' ...
                   'MovingWindowPVal=%s durationWindow=%s useSegments=%s ' ...
                   'UseLOO=%s ApplyFDR=%s nBoot=%s CategoryLevels=%s Lock=%s'], ...
        num2str(f('BaseRespWindow', NaN)),  num2str(f('maxCategory', NaN)), ...
        num2str(f('SpatialGridMode', NaN)), num2str(f('MovingWindowPVal', NaN)), ...
        num2str(f('durationWindow', [])),   num2str(f('useSegments', NaN)), ...
        num2str(f('UseLOO', NaN)),          num2str(f('ApplyFDR', NaN)), ...
        num2str(f('nBoot', NaN)),           flatten(f('CategoryLevels', "all")), ...
        flatten(f('Lock', {})));
end


function v = getFieldOr(s, name, default)
% getFieldOr  Field value if present, otherwise a default (older caches may
%   predate some parameters).
    if isstruct(s) && isfield(s, name)
        v = s.(name);                                         % stored value
    else
        v = default;                                          % parameter did not exist yet
    end
end


function str = flatten(v)
% flatten  Render a nested CategoryLevels / Lock spec as a short string so two
%   specifications can be compared by eye and by string equality.
    if isempty(v)
        str = '{}';                                           % no focus applied
    elseif iscell(v)
        parts = cellfun(@flatten, v, 'UniformOutput', false);  % recurse into nested cells
        str   = ['{' strjoin(parts, ',') '}'];
    elseif isstring(v) || ischar(v)
        str = char(strjoin(string(v), '|'));                   % string array -> 'a|b'
    else
        str = num2str(v(:)');                                  % numeric vector
    end
end


function obj = firstAvailableObj(NP, stims)
% firstAvailableObj  First constructible analysis object among stims, used
%   only to reach the printFig method (its output folder is recording-derived).
    obj = [];                                                 % default: none available
    for si = 1:numel(stims)
        candidate = buildObj(NP, stims{si});                  % [] when absent
        if ~isempty(candidate)
            obj = candidate;                                  % first one wins
            return
        end
    end
end
