function T = extractExampleNeurons(cacheFile, params)
% extractExampleNeurons  Pick MIDDLE-of-distribution example neurons per stimulus.
%
%   From the per-neuron table AllExpAnalysis caches (TableStimComp, which holds
%   only responsive neurons), select the few neurons whose response metric sits
%   closest to the MEDIAN of each stimulus's distribution — i.e. representative,
%   not extreme, examples for a figure. For each pick it resolves the Phy cluster
%   ID and the globally-unique experiment ID (by loading that neuron's NP), so you
%   can plot the rasters yourself and then hand-build a markUnits cell.
%
%   IMPORTANT — responsiveness gate: TableStimComp stores the OR-mask neuron set
%   (responsive to at least ONE compared stimulus) replicated across every stimulus
%   column. So a stimulus's block contains neurons responsive only to the OTHER
%   stimulus, whose z here is ~0. Before taking the median, each stimulus is gated to
%   neurons with Z-score >= params.minZ (default 1.64, ~one-tailed p=0.05 for the
%   right-tailed enhancement test) so the "middle" is a genuinely responsive example.
%
%   The returned table exposes, per neuron:
%     animal, realInsertion (original insertion, resets across animals),
%     Exp (experiment ID, globally unique — the loadNPclassFromTable key),
%     PhyID (Phy cluster ID), plus NeurID and the metric value for reference.
%
%   EXAMPLE
%     T = extractExampleNeurons( ...
%           "W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\Ex_1-97_Combined_MB-RG.mat");
%     % then plot T.PhyID rasters yourself and build markUnits from the winners.

arguments
    cacheFile                                       % path to AllExpAnalysis *.mat cache, OR a preloaded struct, OR a TableStimComp table
    params.metric      string  = "Z-score"          % ranking column: "Z-score" (default) or "SpkR"
    params.nPer        double  = 5                   % candidate neurons to return per stimulus (around the median)
    params.percentile  double  = 60                 % "middle" of the distribution = 50th percentile (exposed for flexibility)
    params.minZ        double  = 1.5                % responsiveness gate on Z-score for THIS stimulus (~one-tailed p=0.05, right-tailed enhancement). See note below.
    params.stimuli     string  = ["MB","RG"]          % stimulus categories to sample
    params.expList     double  = []                  % experiment integers; only needed if a bare table/struct without expList is passed
    params.plot        logical = false               % also draw each candidate's raster via that experiment's plotRaster
    params.PaperFig    logical = false               % paper-style figures; also saves them into a "neuron candidates_stimA_stimB" subfolder
    params.outRoot     string  = ""                  % root for the save subfolder (default: folder of cacheFile, else pwd)
end

% =========================================================================
% 1. Load TableStimComp and the experiment list
% =========================================================================
if istable(cacheFile)                                % caller passed the table directly
    T0      = cacheFile;                             % use it as-is
    expList = params.expList;                        % expList must then come from the argument
elseif isstruct(cacheFile)                           % caller passed a preloaded cache struct
    assert(isfield(cacheFile,'TableStimComp'), 'Struct lacks TableStimComp.'); % must contain the per-neuron table
    T0      = cacheFile.TableStimComp;               % pull the per-neuron table
    expList = pickExpList(cacheFile, params.expList);% prefer struct.expList, else the argument
else                                                 % caller passed a path to the .mat cache
    S = load(cacheFile);                             % load the cached AllExpAnalysis output
    assert(isfield(S,'TableStimComp'), 'Cache lacks TableStimComp — not an AllExpAnalysis cache?'); % validate contents
    T0      = S.TableStimComp;                       % per-neuron table (responsive neurons only)
    expList = pickExpList(S, params.expList);        % experiment integers saved alongside the table
end

% =========================================================================
% 2. Validate columns and inputs (fail loud, per project standard)
% =========================================================================
metric   = params.metric;                            % ranking column name (may contain a hyphen)
needCols = ["animal","realInsertion","insertion","stimulus","NeurID", metric]; % columns this function relies on
missing  = needCols(~ismember(needCols, string(T0.Properties.VariableNames))); % any required column absent?
assert(isempty(missing), 'TableStimComp is missing required column(s): %s', strjoin(missing, ', ')); % clear error listing them
assert(~isempty(expList), ...                        % Phy/Exp resolution cannot proceed without the experiment list
    'expList is empty — pass a cache with S.expList or supply params.expList so Phy IDs can be resolved.');

% =========================================================================
% 3. Select the middle-of-distribution neurons for each stimulus
% =========================================================================
picks = table();                                     % accumulator for chosen rows across all stimuli
for stim = params.stimuli                            % one stimulus at a time (e.g. "MB", then "RG")
    rows = T0(T0.stimulus == stim, :);               % all responsive neurons for this stimulus
    if isempty(rows)                                 % stimulus not present (e.g. category/level-mode composite labels)
        warning('No rows for stimulus "%s" — skipping.', stim); %#ok<*WNTAG> % tell the user and move on
        continue                                     % nothing to sample for this stimulus
    end

    % Responsiveness gate: TableStimComp replicates the OR-mask neuron set across
    % EVERY stimulus column, so a stimulus's block includes neurons responsive only
    % to the OTHER stimulus (their z here can be ~0). Restrict to neurons genuinely
    % responsive to THIS stimulus before taking the median, else the "middle" lands
    % on a non-responsive unit (e.g. RG median z ~0.2 without this gate).
    zcol  = rows.('Z-score');                        % responsiveness measure for THIS stimulus (always Z-score, even if ranking by SpkR)
    m     = rows.(char(metric));                     % ranking metric per neuron (char() handles the hyphen in "Z-score")
    keep  = ~isnan(m) & ~isnan(zcol) & zcol >= params.minZ; % usable AND responsive to this stimulus
    nDrop = sum((zcol < params.minZ) & ~isnan(zcol));% how many were gated out as non-responsive to this stimulus
    rows  = rows(keep, :);                           % keep only responsive, usable rows
    m     = m(keep);                                 % keep matching ranking values
    if nDrop > 0                                     % report the gating so the pool size is transparent
        fprintf('  %s: gated out %d neuron(s) with Z-score < %.2f (not responsive to %s).\n', stim, nDrop, params.minZ, stim);
    end
    if isempty(m)                                    % nothing passed the gate
        warning('No neuron passes the Z-score >= %.2f responsiveness gate for "%s" — skipping.', params.minZ, stim); % report and skip
        continue                                     % cannot rank an empty set
    end

    target = prctile(m, params.percentile);          % the "middle" value of this stimulus's distribution
    [~, ord] = sort(abs(m - target), 'ascend');      % order neurons by closeness to that middle value
    n      = min(params.nPer, numel(ord));           % cannot return more than exist
    if n < params.nPer                               % fewer neurons than requested
        warning('Only %d neuron(s) available for "%s" (requested %d).', n, stim, params.nPer); % note the shortfall
    end
    sel    = ord(1:n);                               % indices of the n neurons nearest the median

    % --- build the per-pick block; str2double(string()) NEVER double() on categoricals ---
    blk = table();                                                    % this stimulus's chosen block
    blk.stimulus      = repmat(string(stim), n, 1);                   % raw stimulus label (e.g. "MB")
    blk.animal        = string(rows.animal(sel));                     % animal ID (categorical -> string)
    blk.realInsertion = str2double(string(rows.realInsertion(sel)));  % original insertion number (resets across animals)
    blk.insertionSeq  = str2double(string(rows.insertion(sel)));      % table's sequential insertion index (reference only; != Exp)
    blk.NeurID        = str2double(string(rows.NeurID(sel)));         % good-unit index (needed later for markUnits)
    blk.(matlab.lang.makeValidName(metric)) = m(sel);                % metric value (column named e.g. "Zscore" or "SpkR")
    blk.distToMedian  = abs(m(sel) - target);                        % |metric - median|; 0 = most representative

    picks = [picks; blk];                                            %#ok<AGROW> % append this stimulus's block
end

assert(~isempty(picks), 'No candidate neurons selected for any requested stimulus.'); % nothing to resolve -> stop clearly

% =========================================================================
% 4. Resolve experiment ID (Exp) and Phy cluster ID (PhyID) per pick
% =========================================================================
picks.Exp   = nan(height(picks), 1);                 % globally-unique experiment ID (loadNPclassFromTable key)
picks.PhyID = nan(height(picks), 1);                 % Phy cluster ID

% Target (animal, realInsertion) pairs we must locate among the experiments.
pairKey   = @(a,ri) sprintf('%s|%g', a, ri);         % canonical key for one (animal, realInsertion) pair
wantKeys  = arrayfun(@(k) pairKey(picks.animal(k), picks.realInsertion(k)), ...
                     (1:height(picks))', 'UniformOutput', false);    % key per pick row
uniqueKeys = unique(wantKeys);                        % distinct experiments we need to load
resolved   = containers.Map('KeyType','char','ValueType','any');     % key -> struct(exp, phy) once found

for e = expList(:)'                                  % scan the experiment list
    if resolved.Count == numel(uniqueKeys), break; end % stop early once every needed experiment is found
    try                                              % loading may fail for a stale/missing recording
        NP = loadNPclassFromTable(e);                % Neuropixels object for this experiment
    catch ME                                         % could not load this one
        warning('Exp %g: could not load NP (%s) — skipping.', e, ME.message); % report and continue scanning
        continue
    end

    % --- parse animal / realInsertion exactly as AllExpAnalysis / findFalseNegAndPos do ---
    a = string(regexp(NP.recordingName, 'PV\d+', 'match', 'once')); % Pogona animal token
    if a == "", a = string(regexp(NP.recordingName, 'SA\d+', 'match', 'once')); end % else Laudakia token
    tok    = strsplit(NP.recordingName, '_');         % split the recording name on underscores
    realIns = str2double(tok{end});                   % insertion = trailing recording-name token
    key    = pairKey(a, realIns);                     % this experiment's (animal, realInsertion) key

    if ~ismember(key, uniqueKeys) || isKey(resolved, key) % not needed, or already resolved
        continue                                      % skip experiments that hold none of our picks
    end

    % --- good-unit Phy IDs (order matches NeurID); shared across stimuli of this recording ---
    stimHere = picks.stimulus(find(strcmp(wantKeys, key), 1)); % a stimulus present for this experiment's picks
    try                                                        % object build / sorting read may fail
        obj     = buildStimObj(NP, stimHere);                 % stimulus-specific analysis object
        ps      = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder); % Phy sorting table
        phy_IDg = ps.phy_ID(string(ps.label') == 'good');     % Phy IDs of good units, in NeurID order
        resolved(key) = struct('exp', e, 'phy', phy_IDg(:), 'NP', NP); % cache resolution + NP handle (reused for plotting)
    catch ME                                                   % resolution failed for this experiment
        warning('Exp %g (%s): could not resolve Phy IDs (%s).', e, a, ME.message); % report; picks stay NaN
    end
end

% --- fill Exp / PhyID for every pick from the resolved map ---
for k = 1:height(picks)                              % each chosen neuron
    key = wantKeys{k};                               % its (animal, realInsertion) key
    if ~isKey(resolved, key)                         % its experiment was never resolved
        warning('Pick %s ins %g NeurID %g: experiment not found in expList — Exp/PhyID left NaN.', ...
            picks.animal(k), picks.realInsertion(k), picks.NeurID(k)); % flag the gap
        continue                                     % leave Exp/PhyID as NaN
    end
    r = resolved(key);                               % struct(exp, phy) for this experiment
    assert(picks.NeurID(k) >= 1 && picks.NeurID(k) <= numel(r.phy), ... % ALIGNMENT INVARIANT: NeurID must index the good-unit list
        'NeurID %g out of range (exp %g has %d good units) — good-unit ordering mismatch.', ...
        picks.NeurID(k), r.exp, numel(r.phy));
    picks.Exp(k)   = r.exp;                          % globally-unique experiment ID
    picks.PhyID(k) = r.phy(picks.NeurID(k));         % Phy cluster ID (reverse of plotRaster's ismember)
end

% =========================================================================
% 5. Order, report, return
% =========================================================================
T = sortrows(picks, {'stimulus','distToMedian'});    % group by stimulus, most-representative first

fprintf('\nextractExampleNeurons: %d candidate(s) around the %gth percentile of %s.\n', ...
    height(T), params.percentile, metric);           % header line
for stim = params.stimuli                            % per-stimulus summary
    ts = T(T.stimulus == stim, :);                   % this stimulus's picks
    if isempty(ts), continue; end                    % skip stimuli with no picks
    fprintf('  %s: median %s = %.3f | %d neuron(s): %s\n', stim, metric, ...
        median(ts.(matlab.lang.makeValidName(metric))), height(ts), ...  % median metric across picks
        strjoin(compose('%s-ins%g-exp%g-phy%g', ts.animal, ts.realInsertion, ts.Exp, ts.PhyID), ', ')); % one token per pick
end

% =========================================================================
% 6. (optional) Plot each candidate's raster via its experiment's plotRaster
% =========================================================================
if params.plot                                       % only when plotting was requested
    % Save target: subfolder "neuron candidates_stimA_stimB" under outRoot (paper figs only).
    if params.outRoot == ""                          % no explicit root given
        if ischar(cacheFile) || (isstring(cacheFile) && isscalar(cacheFile)) % a path was passed
            outRoot = string(fileparts(cacheFile));  % default: the cache's own folder
        else                                         % struct/table input has no path
            outRoot = string(pwd);                   % fall back to the current folder
        end
    else
        outRoot = params.outRoot;                    % caller-specified root
    end
    saveDir = fullfile(outRoot, "neuron candidates_" + strjoin(params.stimuli, "_")); % requested subfolder name
    if params.PaperFig && ~exist(saveDir, 'dir'), mkdir(saveDir); end % create only when we will save

    okRows = T(~isnan(T.Exp), :);                     % only picks whose experiment resolved
    gExp = findgroups(okRows.Exp, okRows.stimulus);  % group by experiment x stimulus -> one shared analysis object per group
    for g = 1:max(gExp)                               % each (experiment, stimulus) group
        sub  = okRows(gExp == g, :);                  % the picks in this group
        stim = sub.stimulus(1);                       % stimulus for this group
        exp  = sub.Exp(1);                            % experiment for this group
        key  = pairKey(char(sub.animal(1)), sub.realInsertion(1)); % (animal, realInsertion) key -> cached NP
        NPh  = resolved(key).NP;                       % reuse the NP loaded during resolution

        try
            obj = buildStimObj(NPh, stim);            % this experiment's analysis object (shared across its neurons below)
        catch ME
            warning('Exp %g %s: could not build analysis object (%s) — skipping.', exp, stim, ME.message); % report and continue
            continue
        end

        % One plotRaster call PER NEURON, not batched: plotRaster (both MB's and
        % RG's) closes every figure except the LAST neuron in its exNeuronsPhyID
        % list, so a batched multi-neuron call would silently drop or mislabel
        % every earlier pick's figure(s). Calling per-neuron makes each newFigs
        % capture unambiguously belong to that one phy ID.
        for r = 1:height(sub)                          % each neuron picked for this experiment+stimulus
            phyId = sub.PhyID(r);                       % this neuron's own Phy cluster ID

            before = findall(0, 'Type', 'figure');      % figures open before this neuron's call
            try
                plotPicks(obj, stim, phyId);             % draw this one neuron's raster (RG also opens a raw-data panel)
            catch ME
                warning('Exp %g %s phyID %g: plotRaster failed (%s) — skipping.', exp, stim, phyId, ME.message); % report and continue
                continue
            end
            if ~params.PaperFig, continue; end          % nothing to save unless paper figures requested
            newFigs = setdiff(findall(0, 'Type', 'figure'), before); % figure(s) this neuron's call created (1 for MB, 2 for RG)
            for fi = 1:numel(newFigs)                    % save each figure belonging to this neuron
                fname = sprintf('%s_%s_ins%g_exp%g_U%d.jpg', stim, sub.animal(1), sub.realInsertion(1), exp, phyId); % base name, this neuron's own phy ID
                if numel(newFigs) > 1, fname = strrep(fname, '.jpg', sprintf('_%d.jpg', fi)); end % disambiguate raster vs raw-data panel
                % exportgraphics on a network share (W:) can hit a transient sharing
                % violation when something else (Explorer thumbnail generation, AV/
                % indexer) briefly opens a just-written .jpg — libjpeg surfaces that as
                % a generic "out of disk space?" write error even with ample free space
                % (confirmed: 15.9 TB free on W:, and most neurons in this same batch
                % exported fine). Retry a few times before giving up; a genuine failure
                % (bad path, truly full disk) will keep failing across all attempts.
                nAttempts = 4;
                for attempt = 1:nAttempts
                    try
                        exportgraphics(newFigs(fi), fullfile(saveDir, fname), 'Resolution', 300); % write into the subfolder
                        break                             % succeeded — stop retrying
                    catch ME
                        if attempt == nAttempts           % out of retries — surface it
                            warning('Export failed after %d attempt(s) for %s: %s', nAttempts, fname, ME.message);
                        else                               % transient — brief pause then retry
                            pause(1);
                        end
                    end
                end
            end
            % plotRaster's own figure-closing logic only fires for a multi-neuron batch
            % ("close all but the last"), which never applies now that each call is a
            % single neuron — so nothing else closes these figures. Close them here once
            % saved, or they accumulate for the whole run (up to nPer x numel(stimuli)
            % neurons, 2 figures each for RG) and exhaust graphics/memory resources.
            close(newFigs);
        end
    end
    if params.PaperFig, fprintf('Saved candidate rasters to:\n  %s\n', saveDir); end % tell the user where they went
end

end

% =========================================================================
% LOCAL HELPERS
% =========================================================================
function obj = buildStimObj(NP, stim)
% buildStimObj  Analysis object for a stimulus abbreviation (mirrors plotFalseNegPosNeurons).
    switch stim                                          % dispatch on the stimulus label
        case "RG",   obj = rectGridAnalysis(NP);         % rectified grid (RG/SB)
        case "MB",   obj = linearlyMovingBallAnalysis(NP);% moving ball
        case "MBR",  obj = linearlyMovingBarAnalysis(NP); % moving bar
        case {"SDGs","SDGm"}, obj = StaticDriftingGratingAnalysis(NP); % grating
        case "NI",   obj = imageAnalysis(NP);            % natural image
        case "NV",   obj = movieAnalysis(NP);            % natural movie
        case "FFF",  obj = fullFieldFlashAnalysis(NP);   % full-field flash
        otherwise,   error('Unknown stimulus "%s".', stim); % unrecognised label
    end
end

function plotPicks(obj, stim, phyIds)
% plotPicks  Call the experiment's plotRaster for a set of Phy IDs, per stimulus.
%   Argument sets follow the user's canonical MB / RG example calls, EXCEPT
%   PaperFig=false and overwrite=false: plotRaster's printFig self-saves on EITHER
%   PaperFig=true (to Paper_figs) OR overwrite=true (to visualStimPlotsFolder), and
%   neither of those internal filenames includes the phy ID. Forcing both off makes
%   this function the sole writer, so rasters land only in the subfolder, named
%   with the phy ID that was actually requested.
    switch stim                                          % stimulus-specific plotRaster arguments
        case "MB"                                        % moving ball: multi-speed, needs per-neuron speed pick
            % plotRaster's own 'speed'="max" default picks the HIGHEST-NUMBERED
            % Speed field (a naming coincidence), not the neuron's best-responding
            % speed. TableStimComp's MB responsiveness (what selected this neuron
            % as an example in the first place) instead picks the speed with the
            % LOWEST p-value per neuron (AllExpAnalysis.m's extractStimData). Without
            % this fix, the raster title could show a non-significant Speed2 p-value
            % for a neuron that is actually significant (and was selected) at Speed1.
            bestSpeed = resolveBestSpeed(obj, phyIds);    % per-neuron best-p Speed field, as a "1"/"2"/... string
            obj.plotRaster('exNeuronsPhyID', phyIds, 'overwrite', false, ...
                'MergeNtrials', 3, 'PaperFig', false, 'OneLuminosity', 'white', 'speed', bestSpeed); % MB example call
        case "MBR"                                       % moving bar: AllExpAnalysis always uses Speed1 for this
            % stimulus (extractStimData's 'MBR' case has no multi-speed selection),
            % so plotRaster's "max" default already matches — no per-neuron pick needed.
            obj.plotRaster('exNeuronsPhyID', phyIds, 'overwrite', false, ...
                'MergeNtrials', 3, 'PaperFig', false, 'OneLuminosity', 'white'); % MBR example call
        case "RG"                                        % rectified grid (RG/SB)
            obj.plotRaster('MergeNtrials', 1, 'overwrite', false, 'AllResponsiveNeurons', false, ...
                'exNeuronsPhyID', phyIds, 'selectedLum', 255, 'oneTrial', true, 'PaperFig', false); % RG example call
        otherwise                                        % other stimuli: minimal generic call
            obj.plotRaster('exNeuronsPhyID', phyIds, 'overwrite', false, ...
                'AllResponsiveNeurons', false, 'PaperFig', false); % generic fallback
    end
end

function speedStr = resolveBestSpeed(obj, phyId)
% resolveBestSpeed  Per-neuron best-p Speed field for MB, mirroring
%   AllExpAnalysis.m's extractStimData (argmin of pvalsResponse across Speed1,
%   Speed2, ...), so the raster title's p-value reflects the same speed that
%   drove this neuron's inclusion as a responsive MB example, rather than
%   plotRaster's own "max" default (highest-numbered Speed field).
    Stats = obj.StatisticsPerNeuron();                        % bare call -> returns cached results only, never recomputes
    speedFields = fieldnames(Stats);                          % all top-level fields of the cached stats struct
    speedFields = speedFields(contains(speedFields, 'Speed')); % keep only Speed1, Speed2, ... (drop 'params' etc.)
    if numel(speedFields) <= 1                                % single-speed experiment -> nothing to choose between
        speedStr = "max";                                     % plotRaster's own default already resolves correctly
        return
    end

    p = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder); % phy sorting table (same source plotRaster itself uses)
    phy_IDg = p.phy_ID(string(p.label') == 'good');               % phy IDs of good units, in unit-index order
    unitIdx = find(phy_IDg == phyId, 1);                          % this neuron's unit index within the good-unit list
    assert(~isempty(unitIdx), ...                                 % ALIGNMENT INVARIANT: phy ID must exist among good units
        'resolveBestSpeed: phy ID %d not found among good units for %s.', phyId, obj.dataObj.recordingName);

    allP = nan(1, numel(speedFields));                        % this neuron's p-value at each speed field
    for iS = 1:numel(speedFields)                             % one speed field at a time
        allP(iS) = Stats.(speedFields{iS}).pvalsResponse(unitIdx); % pvalsResponse is what plotRaster's title itself displays
    end
    [~, bestIdx] = min(allP);                                 % lowest p-value wins (matches extractStimData's selection)
    speedStr = extractAfter(speedFields{bestIdx}, "Speed");   % "Speed2" -> "2", the format plotRaster's 'speed' param expects
end

function el = pickExpList(S, argExpList)
% pickExpList  Prefer the cache's saved expList, else the caller-supplied one.
    if isfield(S,'expList') && ~isempty(S.expList)   % cache stored the experiment list
        el = S.expList;                              % use it
    else                                             % cache had none
        el = argExpList;                             % fall back to the params argument
    end
end
