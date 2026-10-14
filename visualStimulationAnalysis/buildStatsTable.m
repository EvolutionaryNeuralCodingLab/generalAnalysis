function [T, outPaths] = buildStatsTable(spec, params)
% buildStatsTable  Supplementary statistics table across saved AllExpAnalysis runs.
%
%   [T, outPaths] = buildStatsTable(spec, params)
%
%   Recomputes, for every requested comparison, the same statistics that appear
%   on the swarm panels - observed mean and median, hierarchical-bootstrap 95%
%   CI, Cohen's d / d_z with CI, P(>0) with CI, and the exact p-value - and
%   emits them as one table with an explicit exclusion-criteria column.
%
%   The numbers come from the SAME routines that make the figures
%   (bootstrapPairDifference, bootstrapFractionDifference and
%   plotSwarmBootstrapWithComparisons), not from a parallel implementation, so
%   the table cannot drift away from the panels it documents.
%
%   INPUTS
%     spec    Nx3 cell array, one row per comparison:
%               {figureLabel, structPath, pairsSpec}
%             figureLabel - char, e.g. 'Figure 1G'
%             structPath  - full path to an *_analysisStruct.mat written by
%                           AllExpAnalysis with PaperFig=true
%             pairsSpec   - "all" for every pairwise combination, or an Mx2
%                           cell of stimulus-label pairs to restrict to
%     params  Name-value:
%             outDir    - output directory (default: folder of the first struct)
%             outName   - base filename (default 'SupplementaryStatsTable')
%             nBootPval - bootstrap replicates for p-values (default 1e5)
%             pvalSeed  - RNG seed applied before every p-value (default 0)
%             metrics   - metrics to include (default all three)
%             write     - false to compute without writing files
%
%   OUTPUT
%     T         table, one row per (comparison x metric x pair)
%     outPaths  cellstr of files written
%
%   NOTE ON n: it means different things per metric, which is precisely why
%   this table is worth having. For Z-score and SpkR, n counts NEURONS (paired
%   within insertion by NeurID). For FracResponsive, n counts INSERTIONS. The
%   nUnitIs column states which, and nAnimals/nInsertions disambiguate further.
%
%   See also: AllExpAnalysis, bootstrapPairDifference,
%             bootstrapFractionDifference, plotSwarmBootstrapWithComparisons

arguments
    spec cell
    params.outDir      char    = ''
    params.outName     char    = 'SupplementaryStatsTable'
    params.nBootPval   double  = 1e5
    params.pvalSeed    double  = 0
    params.metrics     cell    = {'Z-score','SpkR','FracResponsive'}
    params.write       logical = true
end

assert(size(spec,2) == 3, ...
    'spec must be Nx3: {figureLabel, structPath, pairsSpec}.');

% Smallest resolvable p: the statistic is two-tailed, so with zero counts on
% one side the minimum non-zero value is 2/nBoot (see formatPValue).
pValRes = 2 / params.nBootPval;

rows = table();                                        % accumulator

for si = 1:size(spec,1)

    figLabel   = spec{si,1};
    structPath = spec{si,2};
    pairsSpec  = spec{si,3};

    fprintf('\n[%d/%d] %s\n  %s\n', si, size(spec,1), figLabel, structPath);

    % Fail loudly rather than skipping: a missing struct would otherwise leave
    % a silent hole in the supplementary table.
    assert(exist(structPath, 'file') == 2, ...
        'buildStatsTable: struct not found for %s:\n  %s', figLabel, structPath);

    S = load(structPath);

    % ---- Normalise exactly as AllExpAnalysis does before plotting --------
    TS = S.TableStimComp;
    TS.('Z-score')(isnan(TS.('Z-score'))) = 0;
    TS.SpkR(isnan(TS.SpkR))               = 0;
    TS.animal    = categorical(cellstr(string(TS.animal)));
    TS.insertion = categorical(cellstr(string(TS.insertion)));
    TS.stimulus  = categorical(cellstr(string(TS.stimulus)));

    % ---- Resolve which pairs to report -----------------------------------
    compLabels = cellstr(categories(TS.stimulus));
    if ischar(pairsSpec) || isstring(pairsSpec)
        assert(strcmpi(string(pairsSpec), "all"), ...
            'pairsSpec must be "all" or an Mx2 cell of stimulus pairs.');
        pairsUse = nchoosek(compLabels, 2);
    else
        pairsUse = pairsSpec;
        % Every requested label must exist, or the row would silently vanish
        missing = setdiff(unique(pairsUse(:)), compLabels);
        assert(isempty(missing), ...
            ['buildStatsTable: %s requests stimuli absent from the struct: %s\n' ...
             '  available: %s'], figLabel, strjoin(missing', ', '), ...
            strjoin(compLabels', ', '));
    end

    exclusionText = describeExclusions(S);             % one statement per comparison

    % ---- Unit-level metrics (Z-score, SpkR) ------------------------------
    for mi = 1:numel(params.metrics)
        metric = params.metrics{mi};
        if strcmp(metric, 'FracResponsive'), continue; end    % handled separately

        rng(params.pvalSeed);                          % reproducible, plot-independent
        pv = zeros(1, size(pairsUse,1));
        for pi = 1:size(pairsUse,1)
            pv(pi) = bootstrapPairDifference(TS, pairsUse(pi,:), ...
                params.nBootPval, metric);
        end

        st   = statsFromPanel(TS, pairsUse, pv, {metric}, pValRes);
        rows = [rows; makeRows(figLabel, metric, st, exclusionText, ...
            'neurons', pValRes)];                              %#ok<AGROW>
    end

    % ---- Fraction-responsive ---------------------------------------------
    if any(strcmp(params.metrics, 'FracResponsive')) && isfield(S, 'TableRespNeurs')

        TR = S.TableRespNeurs;

        % Do NOT re-encode a column that is already categorical. Rebuilding the
        % category list renumbers the animal/insertion ids, which changes the
        % order hierBootMatchFreq enumerates its pools and therefore which
        % resample a given seed produces. The distribution is unchanged but the
        % realization is not, so the table would quote a different p from the
        % figure it documents: measured 0.00044 here vs 0.00052 from
        % AllExpAnalysis on the Figure 1 fraction comparison, same data.
        % AllExpAnalysis passes TableRespNeurs through untouched (unlike
        % TableStimComp, which it does re-encode), so this must match.
        if ~iscategorical(TR.animal)
            TR.animal    = categorical(cellstr(string(TR.animal)));
        end
        if ~iscategorical(TR.insertion)
            TR.insertion = categorical(cellstr(string(TR.insertion)));
        end
        if ~iscategorical(TR.stimulus)
            TR.stimulus  = categorical(cellstr(string(TR.stimulus)));
        end

        % AllExpAnalysis SECTION 7 rewrites SDG angle labels from degrees to
        % radians, but only on TableStimComp -- TableRespNeurs keeps degrees.
        % Saved structs from mode 2/3 SDG-angle runs therefore carry two
        % different label conventions internally, and the pairs resolved from
        % TableStimComp would match nothing here. Reconcile before use.
        TR = reconcileAngleLabels(TR, compLabels, figLabel);

        % Keep only insertions holding every compared stimulus, as SECTION 10 does
        fracLabels = cellstr(categories(TR.stimulus));
        [G, ~]     = findgroups(TR.insertion);
        hasAll     = splitapply(@(s) all(ismember(categorical(fracLabels), s)), ...
            TR.stimulus, G);
        tempTable  = TR(hasAll(G) & ismember(TR.stimulus, categorical(fracLabels)), :);

        % Restrict to the same pairs, where both labels exist in this table
        keepPair  = all(ismember(pairsUse, fracLabels), 2);
        pairsFrac = pairsUse(keepPair, :);

        if isempty(pairsFrac)
            warning('buildStatsTable:noFractionPairs', ...
                '%s: no requested pair exists in TableRespNeurs; skipping FracResponsive.', ...
                figLabel);
        else
            rng(params.pvalSeed);
            pvF = zeros(1, size(pairsFrac,1));
            for pi = 1:size(pairsFrac,1)
                pvF(pi) = bootstrapFractionDifference(TR, pairsFrac(pi,:), ...
                    params.nBootPval);
            end

            stF  = statsFromPanel(tempTable, pairsFrac, pvF, ...
                {'respNeur','totalSomaticN'}, pValRes);
            rows = [rows; makeRows(figLabel, 'FracResponsive', stF, ...
                exclusionText, 'insertions', pValRes)];         %#ok<AGROW>
        end
    end
end

T = rows;

% ---- Write -------------------------------------------------------------
outPaths = {};
if params.write && ~isempty(T)
    outDir = params.outDir;
    if isempty(outDir), outDir = fileparts(spec{1,2}); end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    xlsxPath = fullfile(outDir, [params.outName '.xlsx']);
    csvPath  = fullfile(outDir, [params.outName '.csv']);
    writetable(T, xlsxPath);
    writetable(T, csvPath);
    outPaths = {xlsxPath, csvPath};

    fprintf('\nWrote %d rows to:\n  %s\n  %s\n', height(T), xlsxPath, csvPath);
end

end % buildStatsTable


% =========================================================================
% LOCAL: statsFromPanel
% Obtain per-group and per-pair statistics by going through the same entry
% point the figures use, then discard the transient figures. computeSwarmStats
% is a local function of plotSwarmBootstrapWithComparisons and cannot be called
% directly; routing through the public function guarantees these numbers are
% computed identically to the plotted ones.
% =========================================================================
function st = statsFromPanel(tbl, pairsUse, pvals, valueField, pValRes)

isFraction = numel(valueField) == 2;                   % two fields => ratio mode

[figMain, ~, figDiffs, st] = plotSwarmBootstrapWithComparisons( ...
    tbl, pairsUse, pvals, valueField, ...
    fraction = isFraction, showBothAndDiff = false, diff = false, ...
    plotMeanSem = true, pValRes = pValRes);

% Close the transient figures; only the statistics are wanted
if isgraphics(figMain),  close(figMain);  end
if ~isempty(figDiffs) && isgraphics(figDiffs), close(figDiffs); end

end % statsFromPanel


% =========================================================================
% LOCAL: makeRows
% Flatten one stats struct into table rows, one per pair.
% =========================================================================
function R = makeRows(figLabel, metric, st, exclusionText, unitName, pValRes)

D = st.diffs;
R = table();

for i = 1:height(D)

    % Per-group descriptives for the two sides of this pair
    gA = st.groups(st.groups.group == D.stimA(i), :);
    gB = st.groups(st.groups.group == D.stimB(i), :);

    [pNum, pOp] = formatPValue(D.p(i), pValRes);

    r = table( ...
        string(figLabel), string(metric), ...
        D.stimA(i), D.stimB(i), ...
        D.n(i), D.nAnimals(i), D.nInsertions(i), string(unitName), ...
        pick(gA,'mean'), pick(gA,'median'), ...
        pick(gB,'mean'), pick(gB,'median'), ...
        D.mean(i), D.median(i), D.ciLo(i), D.ciHi(i), ...
        D.effect(i), D.effectLo(i), D.effectHi(i), ...
        D.propPos(i), D.propPosLo(i), D.propPosHi(i), ...
        D.p(i), string(['p' pOp pNum]), ...
        string(exclusionText), ...
        'VariableNames', { ...
            'Figure','Metric','StimA','StimB', ...
            'n','nAnimals','nInsertions','nUnitIs', ...
            'meanA','medianA','meanB','medianB', ...
            'diffMean','diffMedian','diffCIlo','diffCIhi', ...
            'CohensD','dCIlo','dCIhi', ...
            'PropAgtB','propCIlo','propCIhi', ...
            'pValue','pFormatted','ExclusionCriteria'});

    R = [R; r];                                        %#ok<AGROW>
end

end % makeRows


% =========================================================================
% LOCAL: reconcileAngleLabels
% Bring TableRespNeurs stimulus labels into the same convention as
% TableStimComp for mode 2/3 SDG-angle runs.
%
% ROOT CAUSE: AllExpAnalysis SECTION 7 converts SDG angle labels from degrees
% to radians (0->1.57, 90->0, 180->4.71, 270->3.14) but applies the rewrite
% ONLY to S.TableStimComp. S.TableRespNeurs is saved with the original degree
% labels, so one file ends up holding two conventions. This is a compatibility
% shim for structs already on disk; the durable fix is to relabel both tables
% in AllExpAnalysis.
%
% The conversion is applied only when it demonstrably helps: if the label sets
% already agree it is a no-op, and if applying the map does not produce a match
% the table is returned untouched with a warning, so this can never silently
% scramble labels that were already correct.
% =========================================================================
function TR = reconcileAngleLabels(TR, compLabels, figLabel)

trLabels = cellstr(categories(TR.stimulus));

% Already consistent with the comparison labels: nothing to do
if isempty(setdiff(trLabels, compLabels))
    return
end

% Only meaningful for SDG angle labels
if ~any(contains(string(trLabels), "SDG")) || ~any(contains(string(trLabels), "ang"))
    return
end

% Apply the SECTION 7 map to the SDG rows, exactly as AllExpAnalysis does
s      = string(TR.stimulus);
idxSDG = contains(s, "SDG");
sNew   = s;
sNew(idxSDG) = replace(s(idxSDG), ...
    ["0",   "90", "180", "270"], ...
    ["1.57","0",  "4.71","3.14"]);

% Accept only if the rewrite actually reconciles the two conventions
if isempty(setdiff(cellstr(unique(sNew)), compLabels))
    TR.stimulus = categorical(cellstr(sNew));
    fprintf(['  %s: TableRespNeurs angle labels were in degrees while ' ...
             'TableStimComp used radians; reconciled for the fraction rows.\n'], figLabel);
else
    warning('buildStatsTable:angleLabelMismatch', ...
        ['%s: TableRespNeurs labels (%s) do not match TableStimComp (%s) and ' ...
         'the degree-to-radian map does not reconcile them. Fraction rows may ' ...
         'be skipped.'], figLabel, strjoin(trLabels', ', '), strjoin(compLabels', ', '));
end

end % reconcileAngleLabels


function v = pick(g, field)
% pick  Scalar field from a one-row group table; NaN when the group is absent.
if height(g) == 1
    v = g.(field);
else
    v = NaN;
end
end


% =========================================================================
% LOCAL: describeExclusions
% Build the inclusion/exclusion statement for one comparison from the
% parameters actually used, so the table documents what was run rather than
% what someone remembers running.
% =========================================================================
function txt = describeExclusions(S)

p     = S.params;
parts = {};

% --- Responsiveness criterion (the OR-mask) ---
if isfield(p,'useTtest') && p.useTtest
    testName = 'paired t-test on the response window';
else
    testName = 'max-statistic sign-flip permutation test, right-tailed';
end
parts{end+1} = sprintf(['Units: somatic units responsive to at least one ' ...
    'compared stimulus (OR-mask), %s, alpha=%g'], testName, p.threshold);

% --- Multiple-comparison handling ---
fdrBits = {};
if isfield(p,'applyFDR') && p.applyFDR
    fdrBits{end+1} = 'BH-FDR inside the statistic';
end
if isfield(p,'useFDR') && p.useFDR
    fdrBits{end+1} = 'BH-FDR post-hoc';
end
if isempty(fdrBits)
    parts{end+1} = ['No FDR correction applied; the max-statistic null ' ...
        'controls family-wise error across categories and positions within ' ...
        'a stimulus, but not across stimulus types'];
else
    parts{end+1} = ['Correction: ' strjoin(fdrBits, ' + ')];
end

% --- Response window ---
winTxt = sprintf('Response window: base %g ms', p.BaseRespWindow);
if isfield(p,'durationWindow') && ~isempty(p.durationWindow)
    winTxt = sprintf('%s, sliding %g ms', winTxt, p.durationWindow);
end
parts{end+1} = winTxt;

if isfield(p,'maxCategory') && p.maxCategory
    parts{end+1} = ['Statistic maximised across categories, with the null ' ...
        'taking the same maximum and nullMean subtracted to remove the ' ...
        'winner-take-max bias'];
end

if isfield(p,'SpatialGridMode') && p.SpatialGridMode
    parts{end+1} = 'Spatial-grid mode: base window forced to 200 ms';
end

% --- Trial focus (Lock) ---
if isfield(p,'Lock') && ~isempty(p.Lock)
    lockTxt = {};
    for i = 1:numel(p.Lock)
        L = p.Lock{i};
        for k = 1:2:numel(L)-1
            lockTxt{end+1} = sprintf('%s=%g', string(L{k}), L{k+1});  %#ok<AGROW>
        end
    end
    if ~isempty(lockTxt)
        parts{end+1} = ['Trials restricted to ' strjoin(unique(lockTxt), ', ')];
    end
end

% --- Insertion-level requirement ---
parts{end+1} = ['Insertions: only those containing every compared stimulus; ' ...
    'fraction-responsive rows additionally require both stimuli within the ' ...
    'same insertion'];

% --- Mode-specific requirements ---
switch S.mode
    case 2
        parts{end+1} = ['Within-stimulus category comparison; recordings with ' ...
            'fewer than 2 levels of the category excluded'];
    case 3
        parts{end+1} = ['Specific-level comparison; recordings lacking any ' ...
            'requested level within a single session excluded'];
end

txt = strjoin(parts, '. ');

end % describeExclusions
