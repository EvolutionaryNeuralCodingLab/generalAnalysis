function D = findFalseNegAndPos(cacheFile, params)
% findFalseNegAndPos  Shortlist candidate FALSE NEGATIVES and FALSE POSITIVES
%   in a union-mode raster cache, then confirm each against the actual
%   StatisticsPerNeuron numbers. Output is labelled by animal, insertion and
%   Phy cluster ID.
%
%   TWO STAGES
%   ----------
%   Stage 1 — RASTER SCREEN (what you see in the figure):
%     For each neuron, take the peak post-onset z and a sustained-bin count
%     in each panel. Calibrate a per-panel "responsive-looking" bar from the
%     "both" group (cells the test calls responsive to both stimuli, hence
%     genuinely responsive in both panels): bar = refPrctile-th percentile
%     of their peak z. Then
%       FALSE NEGATIVE candidate for Y: neuron NOT in Y's group, but its
%          Y-panel peak >= bar(Y) and sustained >= nBinsThresh.
%       FALSE POSITIVE candidate for X: neuron IS in X's group, but its
%          X-panel peak <  bar(X) and sustained <  nBinsThresh.
%
%   Stage 2 — CONFIRMATION (the actual test numbers; confirm=true):
%     For each shortlisted neuron, pull from its StatisticsPerNeuron cache
%     (matched by Phy ID; MB = min p across speeds) the p-value and ZScoreU
%     for the relevant stimulus — the OFF stimulus for an FN, the ON stimulus
%     for an FP — and assign a verdict:
%       FN, alpha<=p<fnMarginBand            -> "FN:marginal"   (near-miss; flips with a matched/lenient test)
%       FN, p>=fnMarginBand                  -> "FN:discrepant" (firm miss but visible: focal/window — inspect)
%       FN, p<alpha                          -> "FN:CACHE-MISMATCH" (raster label disagrees with current stats)
%       FP, fpBorderlineFrac*alpha<=p<alpha  -> "FP:borderline" (barely significant — the genuine FPs)
%       FP, p<fpBorderlineFrac*alpha         -> "FP:firm-sig"   (strong stat, flat raster: response outside the window, e.g. MB>500ms)
%       FP, p>=alpha                         -> "FP:CACHE-MISMATCH"
%
%   Candidates are not verdicts. For a 0-500 ms display the trustworthy
%   directions are RG-FN, RG-FP and MB-FN; MB-FP is mostly the display crop.
%
%   EXAMPLE CALL
%   ------------
%     D = findFalseNegAndPos( ...
%           "W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\Ex_1-97_n44_Raster_MB-RG_union.mat");
%     head(sortrows(D,'Severity','descend'), 20)

arguments
    cacheFile (1,1) string                                  % full path to the *_union.mat cache
    % How "responsive-looking" is judged from the raster (Stage 1):
    params.binPrctile  double = 99    % a post-onset bin is "above noise" if it exceeds this percentile of the baseline (pre-onset) z
    % DURATION-NORMALIZED threshold: a neuron "looks responsive" if the FRACTION
    % of its post-onset bins above noise reaches minFrac. Using a fraction (not a
    % fixed bin count) makes MB (~2.3 s, ~280 bins) and RG (~0.5 s, ~70 bins)
    % judged at equal sensitivity — a fixed count of 10 bins was 3.6% of the MB
    % window but 14% of the RG window, systematically over-flagging short RG
    % responses. minBins is retained only as an absolute floor so very short
    % windows cannot pass on 1-2 chance bins.
    params.minFrac     double = 0.05   % fraction of post-onset bins above noise required to "look responsive"
    params.minBins     double = 3      % absolute floor (bins); the per-panel threshold is max(minBins, round(minFrac*nPost))
    params.confirm         logical = true                    % run Stage 2 (pull p / ZScoreU from StatisticsPerNeuron)
    params.fnMarginBand    double  = 0.20                    % FN with alpha<=p<this is "marginal"
    % An FN candidate must be CORROBORATED by a positive stat effect size
    % (ZScoreU). When the raster's apparent response is a z-score artifact
    % (near-silent baseline inflating the raster z to 10-50 while the stat's
    % effect size sits at ~0 or negative), the neuron is not a genuine miss —
    % the label and the stat's magnitude agree there is no response. Such
    % candidates get verdict "FN:z-inflation" (benign) instead of marginal/
    % discrepant, so they are excluded by plotRaster's verdictFilter. Uses the
    % effect size (ZScoreU), which is independent of the p-value that did the
    % labelling, so this is a magnitude cross-check, not p-on-p circularity.
    params.fnMinZ          double  = 0.5                     % min stat ZScoreU for an FN candidate to count as genuine
    params.fpBorderlineFrac double = 0.5                     % FP with p>=this*alpha is "borderline"
    params.saveCSV         logical = true                    % write the table next to the cache
end

% =========================================================================
% STAGE 1 — RASTER SCREEN
% =========================================================================
S = load(cacheFile);                                        % cached union raster struct
assert(isfield(S,'params'), 'Cache lacks params — not a plotRaster cache?');
stimTypes = string(S.params.stimTypes);                     % stimulus order (e.g. ["MB","RG"])
nStim     = numel(stimTypes);  assert(nStim >= 2, 'Need >= 2 panels.');
preBase   = S.lockedPreBase;                                % baseline length (ms); stim onset sits here
alpha     = S.params.alpha;                                 % classification threshold actually used by the raster
% The raster labels neurons responsive with pValTTest when it was built with
% useTtest=true, otherwise with pvalsResponse. Read the same flag here so
% Stage-2 confirmation compares the IDENTICAL p-value field the labels came
% from — mismatching the field would flag spurious CACHE-MISMATCH verdicts.
useTtest  = isfield(S.params,'useTtest') && S.params.useTtest; % default false for older caches lacking the field

% Per-panel raster + post-onset mask
raster = cell(1,nStim); postIx = cell(1,nStim);
for s = 1:nStim
    f         = matlab.lang.makeValidName(stimTypes(s));    % struct-safe field stem
    raster{s} = S.(f + "_raster");                          % z-scored PSTH [nNeu x nBins]
    postIx{s} = S.lockedEdges{s}(1:end-1) >= preBase;       % post-onset bins
end
f1   = matlab.lang.makeValidName(stimTypes(1));
expV = double(S.(f1 + "_exp"))';                            % experiment ID per neuron  [N x 1]
phyV = double(S.(f1 + "_phy"))';                            % Phy cluster ID per neuron [N x 1]
grpV = string(S.(f1 + "_respGroup"))';                      % responsiveness group      [N x 1]
N    = numel(phyV);

% For each panel, find the noise line from the baseline and count how many
% post-onset bins rise above it. The baseline (pre-onset) period has no
% response, so its bins show how high a bin gets by chance.
isBoth    = grpV == "both";                                % cells the test calls responsive to BOTH stimuli
peakZ     = nan(N,nStim);                                  % strongest post-onset bin (kept for the table; not used to decide)
respCount = zeros(N,nStim);                                % number of post-onset bins above the noise line
looksResp = false(N,nStim);                                % the decision: does the neuron look responsive in this panel?
noiseLine = nan(1,nStim);                                  % per-panel noise level (z), read off the baseline
minBinsVec = zeros(1,nStim);                               % per-panel duration-normalized bin threshold (used again for severity/plot)

for s = 1:nStim
    baseZ        = raster{s}(:, ~postIx{s});               % baseline (pre-onset) z of every neuron = the noise [N x nBaseBins]

    % PER-NEURON noise line: each neuron's own high-percentile baseline bin,
    % instead of one pooled line across all neurons. A pooled line is dominated
    % by high-variance units, so a modest-but-real response in a low-variance
    % unit never clears it — the source of the firm-sig FPs (permutation test
    % firmly responsive, screen calls it flat). Per-neuron mirrors the
    % permutation test's per-unit baseline logic. Compared row-wise below.
    noiseLinePN  = prctile(baseZ, params.binPrctile, 2);   % [N x 1] each neuron's own noise line (NaN for empty/placeholder rows)
    noiseLine(s) = median(noiseLinePN, 'omitnan');         % scalar summary, printout only

    nPost          = sum(postIx{s});                       % number of post-onset bins in this panel
    minBinsVec(s)  = max(params.minBins, round(params.minFrac * nPost)); % proportional threshold, floored at minBins
    postZ          = raster{s}(:, postIx{s});              % post-onset z of every neuron [N x nPost]
    peakZ(:,s)     = max(postZ, [], 2, 'omitnan');         % strongest post-onset bin (reference only)
    respCount(:,s) = sum(postZ > noiseLinePN, 2, 'omitnan');    % bins above THIS neuron's own line (implicit row-wise expansion)
    looksResp(:,s) = respCount(:,s) >= minBinsVec(s);      % responsive if enough bins clear the line (duration-normalized)

    % Check: "both" cells genuinely respond here, so most should clear the
    % line. A low percentage means the noise line is set too high.
    fprintf('  %s: median per-neuron noise line = %.3f z (%gth pctl of baseline), threshold = %d/%d bins (%.1f%%) — captures %.0f%% of "both" cells\n', ...
        stimTypes(s), noiseLine(s), params.binPrctile, minBinsVec(s), nPost, ...
        100*params.minFrac, 100*mean(looksResp(isBoth,s)));
end

% Compare the test's label to how the neuron looks, per panel:
%   false negative = NOT labelled responsive, but looks responsive
%   false positive = labelled responsive, but looks flat
FNfor = strings(N,1); FPfor = strings(N,1); severity = zeros(N,1);
for k = 1:N
    calledResp = (grpV(k) == "both") | (stimTypes == grpV(k));   % stimuli the test labelled this neuron responsive to
    for s = 1:nStim
        if ~calledResp(s) && looksResp(k,s)               % false negative
            FNfor(k) = stimTypes(s);
            severity(k) = max(severity(k), respCount(k,s));               % more bins above noise = stronger miss
        elseif calledResp(s) && ~looksResp(k,s)           % false positive
            FPfor(k) = stimTypes(s);
            severity(k) = max(severity(k), minBinsVec(s) - respCount(k,s));  % further below the panel threshold = more clearly flat
        end
    end
end
isCand = (FNfor ~= "") | (FPfor ~= "");

% =========================================================================
% STAGE 2 — per-experiment resolution: animal/insertion + (optional) stats
% =========================================================================
animal = strings(N,1) + "<unresolved>"; insertion = nan(N,1);
FN_p = nan(N,1); FN_z = nan(N,1); FP_p = nan(N,1); FP_z = nan(N,1);  % confirmation columns
statMap = containers.Map();                                 % key "stim|exp" -> struct(phy,p,z)
statKey = @(st,e) sprintf('%s|%d', st, e);

uex = unique(expV(isCand));
for e = uex(:)'
    try, NP = loadNPclassFromTable(e); catch, continue; end

    % --- animal / insertion (same parse as AllExpAnalysis) ---
    a = string(regexp(NP.recordingName, 'PV\d+', 'match', 'once'));
    if a == "", a = string(regexp(NP.recordingName, 'SA\d+', 'match', 'once')); end
    tok = strsplit(NP.recordingName, '_');
    aThis = a; iThis = str2double(tok{end});               % insertion = trailing recording-name token

    rows = find(isCand & expV == e);
    for k = rows', animal(k) = aThis; insertion(k) = iThis; end

    if ~params.confirm, continue; end

    % --- load stats for each stimulus actually needed for this experiment ---
    needed = unique([FNfor(rows); FPfor(rows)]);            % stimuli referenced by this exp's candidates
    needed = needed(needed ~= "");
    for st = needed'
        try
            obj   = buildObj(NP, st);                       % stimulus-specific analysis object
            % BARE call (no name-value pairs): returns EXACTLY the last-saved
            % StatisticsPerNeuron cache — the same one plotRaster_MultiExp read
            % (also via a bare call) to build the raster labels. Passing explicit
            % params here (e.g. SpatialGridMode/maxCategory) would fail
            % computationParamsMatch against the focused cache AllExpAnalysis
            % wrote, forcing a recompute with a DIFFERENT test AND overwriting
            % the cache on disk — producing spurious CACHE-MISMATCH verdicts and
            % desyncing the stats plotRaster depends on.
            Stats = obj.StatisticsPerNeuron;                % loads cache (identical to the raster's source)
            ps    = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
            phy   = ps.phy_ID(string(ps.label') == 'good'); % good-unit Phy IDs (order matches stats vectors)
            [pv, zv] = extractPZ(Stats, st, useTtest);      % per good-unit p and ZScoreU (MB = min across speeds); p-field matches the raster
            statMap(statKey(st,e)) = struct('phy', phy(:), 'p', pv(:), 'z', zv(:));
        catch ME
            fprintf('  exp %d %s: stats unavailable (%s)\n', e, st, ME.message);
        end
    end
end

% --- fill confirmation p/z per candidate by Phy-ID lookup ---
if params.confirm
    for k = find(isCand)'
        if FNfor(k) ~= ""
            [FN_p(k), FN_z(k)] = lookupPZ(statMap, statKey(FNfor(k), expV(k)), phyV(k));
        end
        if FPfor(k) ~= ""
            [FP_p(k), FP_z(k)] = lookupPZ(statMap, statKey(FPfor(k), expV(k)), phyV(k));
        end
    end
end

% --- verdict ---
Verdict = strings(N,1);
for k = find(isCand)'
    v = strings(1,0);
    if FNfor(k) ~= ""
        if isnan(FN_p(k)),                 v(end+1) = "FN:no-stat";        %#ok<AGROW>
        elseif FN_p(k) <  alpha,           v(end+1) = "FN:CACHE-MISMATCH"; %#ok<AGROW>
        elseif FN_z(k) <= params.fnMinZ,   v(end+1) = "FN:z-inflation";    %#ok<AGROW> effect size does not corroborate the raster's apparent response
        elseif FN_p(k) <  params.fnMarginBand, v(end+1) = "FN:marginal";   %#ok<AGROW>
        else,                              v(end+1) = "FN:discrepant";     %#ok<AGROW>
        end
    end
    if FPfor(k) ~= ""
        if isnan(FP_p(k)),                 v(end+1) = "FP:no-stat";        %#ok<AGROW>
        elseif FP_p(k) >= alpha,           v(end+1) = "FP:CACHE-MISMATCH"; %#ok<AGROW>
        elseif FP_p(k) >= params.fpBorderlineFrac*alpha, v(end+1) = "FP:borderline"; %#ok<AGROW>
        else,                              v(end+1) = "FP:firm-sig";       %#ok<AGROW>
        end
    end
    Verdict(k) = strjoin(v, "; ");
end

% =========================================================================
% Assemble, rank, report
% =========================================================================
D = table(animal, insertion, phyV, expV, grpV, FNfor, FPfor, ...
          peakZ(:,1), peakZ(:,2), FN_p, FN_z, FP_p, FP_z, Verdict, severity, ...
    'VariableNames', {'Animal','Insertion','PhyID','Exp','Group','FN_for','FP_for', ...
        char(stimTypes(1)+"_peakZ"), char(stimTypes(2)+"_peakZ"), ...
        'FN_p','FN_z','FP_p','FP_z','Verdict','Severity'});
D = D(isCand,:);
D = sortrows(D,'Severity','descend');

fprintf('\n%d candidates of %d union neurons.\n', height(D), N);
for s = 1:nStim
    fprintf('  FN for %s: %d   FP for %s: %d\n', ...
        stimTypes(s), sum(D.FN_for==stimTypes(s)), stimTypes(s), sum(D.FP_for==stimTypes(s)));
end
% Genuine vs artifact split for FN: z-inflation candidates have no corroborating
% effect size, so they are display artifacts, not real misses.
if params.confirm
    nZinfl   = sum(contains(D.Verdict,"FN:z-inflation"));                  % raster z-score artifacts (benign)
    nGenuine = sum(contains(D.Verdict,"FN:marginal") | contains(D.Verdict,"FN:discrepant")); % corroborated FN
    fprintf('  FN split: %d genuine (marginal/discrepant), %d z-inflation artifacts (ZScoreU<=%.2g, excluded from genuine).\n', ...
        nGenuine, nZinfl, params.fnMinZ);
end
if params.confirm && any(contains(D.Verdict,"CACHE-MISMATCH"))
    fprintf('  WARNING: %d CACHE-MISMATCH rows — raster labels disagree with current stats (stale cache?).\n', ...
        sum(contains(D.Verdict,"CACHE-MISMATCH")));
end

if params.saveCSV && ~isempty(D)
    [d,b] = fileparts(cacheFile);
    out = fullfile(d, b + "_falseNegPos.csv");
    writetable(D, out); fprintf('Wrote %s\n', out);
end

% Stage-1 picture: how responsive each neuron looks in each panel (bins above
% noise), coloured by the test's label. A point whose colour disagrees with
% where it sits is an FN or FP. Dashed lines mark the per-panel cut-offs
% (now duration-normalized, so the two axes can differ).
figure('Color','w'); hold on
for g = unique(grpV)'
    m = grpV == g;
    scatter(respCount(m,1)+0.15*randn(sum(m),1), respCount(m,2)+0.15*randn(sum(m),1), ...
        16, 'filled', 'DisplayName', char(g));            % small jitter so equal counts don't overlap
end
yl = ylim; xl = xlim;
plot(xl, (minBinsVec(2)-0.5)*[1 1], 'k--', 'HandleVisibility','off');   % stim2 (y-axis) threshold
plot((minBinsVec(1)-0.5)*[1 1], yl, 'k--', 'HandleVisibility','off');   % stim1 (x-axis) threshold
xlabel(sprintf('%s: post-onset bins above noise', stimTypes(1)));
ylabel(sprintf('%s: post-onset bins above noise', stimTypes(2)));
legend('Location','northeastoutside'); box off
title('Colour = test label; position = how responsive it looks; disagreements are FN/FP');

% --- Scatter 2: confirmation — test p vs ZScoreU for each candidate ---
if params.confirm
    figure('Color','w'); hold on
    pAll = [D.FN_p; D.FP_p]; zAll = [D.FN_z; D.FP_z];
    lab  = [repmat("FN",height(D),1); repmat("FP",height(D),1)];
    keep = ~isnan(pAll);
    for L = ["FN","FP"]
        m = keep & lab==L; scatter(max(pAll(m),1e-4), zAll(m), 18, 'filled', 'DisplayName', char(L));
    end
    set(gca,'XScale','log'); xline(alpha,'k--','HandleVisibility','off');
    xlabel('StatisticsPerNeuron p (relevant stimulus)'); ylabel('ZScoreU (effect size)');
    legend('Location','best'); box off
    title('Stage 2: FN should sit at p\geq\alpha with real z; FP at p just <\alpha with low z');
end
end

% =========================================================================
% LOCAL HELPERS
% =========================================================================
function obj = buildObj(NP, stim)
% buildObj  Construct the analysis object for a stimulus abbreviation.
    switch stim
        case "RG",   obj = rectGridAnalysis(NP);
        case "MB",   obj = linearlyMovingBallAnalysis(NP);
        case "MBR",  obj = linearlyMovingBarAnalysis(NP);
        case {"SDGs","SDGm","SDG"}, obj = StaticDriftingGratingAnalysis(NP);
        case "NI",   obj = imageAnalysis(NP);
        case "NV",   obj = movieAnalysis(NP);
        case "FFF",  obj = fullFieldFlashAnalysis(NP);
        otherwise,   error('Unknown stim %s', stim);
    end
end

function [p, z] = extractPZ(Stats, stim, useTtest)
% extractPZ  Per good-unit p-value and ZScoreU for a stimulus.
%   MB/MBR: minimum p across speeds (matches AllExpAnalysis), z at that speed.
%   useTtest selects the p-value FIELD so it matches the one the raster used to
%   label neurons: pValTTest when true, pvalsResponse when false.
    if nargin < 3, useTtest = false; end                    % default: permutation-test p (pvalsResponse)
    pField = 'pvalsResponse';                               % permutation-test p-value (raster default)
    if useTtest, pField = 'pValTTest'; end                  % t-test p-value (raster's useTtest=true path)
    switch stim
        case {"MB","MBR"}
            spF = fieldnames(Stats); spF = spF(startsWith(string(spF),"Speed"));
            if isempty(spF)
                p = Stats.(pField)(:); z = Stats.ZScoreU(:);
            else
                P = nan(numel(Stats.(spF{1}).(pField)), numel(spF)); Z = nan(size(P));
                for k = 1:numel(spF)
                    P(:,k) = Stats.(spF{k}).(pField)(:);
                    Z(:,k) = Stats.(spF{k}).ZScoreU(:);
                end
                [p, idx] = min(P, [], 2);                    % best speed per neuron
                z = Z(sub2ind(size(Z), (1:size(Z,1))', idx));% z at the winning speed
            end
        case "SDGs", p = Stats.Static.(pField)(:); z = Stats.Static.ZScoreU(:);
        case "SDGm", p = Stats.Moving.(pField)(:); z = Stats.Moving.ZScoreU(:);
        otherwise,   p = Stats.(pField)(:); z = Stats.ZScoreU(:);
    end
end

function [p, z] = lookupPZ(statMap, key, phyID)
% lookupPZ  Fetch p and ZScoreU for one Phy ID from a cached stats entry.
    p = NaN; z = NaN;
    if ~isKey(statMap, key); return; end
    s   = statMap(key);
    pos = find(s.phy == phyID, 1);                          % Phy-ID match (robust to positional drift)
    if isempty(pos); return; end
    p = s.p(pos); z = s.z(pos);
end