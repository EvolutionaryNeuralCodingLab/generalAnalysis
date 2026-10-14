function [fig, randiColors, figAllDiffs, stats] = plotSwarmBootstrapWithComparisons(tbl, pairs, pValues, valueField, params)
% PLOTSWARMBOOTSTRAPWITHCOMPARISONS  Per-stimulus swarm plot with hierarchical
%   bootstrap central tendency, uncertainty bar, and pairwise significance brackets.
%
%   [fig, randiColors, figAllDiffs, stats] = plotSwarmBootstrapWithComparisons(tbl, ...
%       pairs, pValues, valueField, params)
%
%   tbl        - One row per observation. Required columns: stimulus, animal,
%                insertion (all categorical). Optional: NeurID (numeric, used
%                in diff mode for neuron-level data), zScore (numeric).
%                Two granularities are auto-detected:
%                  * neuron-level   : multiple rows per (insertion,stimulus)
%                  * insertion-level: one row per (insertion,stimulus)
%   pairs      - Nx2 cell array of stimulus name pairs to compare/test.
%   pValues    - N-element vector of pre-computed p-values for those pairs.
%   valueField - 1-cell {field} for raw value, 2-cell {num,den} for ratio.
%
%   Selected params (see arguments block for full list):
%     nBoot          - bootstrap replicates (default 10000)
%     fraction       - true => valueField{1} ./ valueField{2}
%     diff           - plot per-neuron stimA-stimB difference instead of raw
%     showBothAndDiff- two-tile layout: raw on left, difference on right
%     ciMethod       - 'sem' or 'percentile' (95% bootstrap CI, default)
%     bootGroupVars  - cell of column names defining the bootstrap hierarchy.
%                      Default auto-fills from {'animal','insertion'} based on
%                      what exists in the table. Pass {} explicitly to force a
%                      flat bootstrap.
%     rngSeed        - bootstrap RNG seed for reproducibility (default 0)
%
%   Returns the figure handle, random dot-draw permutation, and (when >1 pair)
%   a second figure showing every pairwise difference in separate tiles.
%
%   stats - struct with two summary tables (the numbers behind the plotted bars):
%     stats.groups - one row per stimulus, in x-axis (category) order.
%     stats.diffs  - one row per entry of `pairs`, in the SAME row order as
%                    `pairs`, with stimA/stimB/p appended.
%   Both tables carry the columns:
%     group       - label ('MB' for a group, 'MB-SB' for a difference)
%     n           - number of non-NaN observations (neurons or insertions)
%     nAnimals    - distinct animals contributing to those observations
%     nInsertions - distinct insertions contributing to those observations
%     mean,median - PLAIN pooled sample statistics of the observations actually
%                   plotted: every observation counts equally. `mean` is the
%                   true centre of the dot cloud and IS the value drawn as the
%                   thick black bar.
%     bootMean    - mean of the hierarchical-bootstrap distribution. DIAGNOSTIC
%                   ONLY, not plotted (see note below).
%     ciLo, ciHi  - bootstrap uncertainty bounds on the mean (params.ciMethod:
%                   95% percentile interval of the bootstrap distribution, or
%                   the observed mean +/- the bootstrap SE).
%     effect      - Cohen's d of the plotted values (mean / SD). In
%                   stats.groups this is a one-sample d against zero; in
%                   stats.diffs the values ARE the paired differences, so it is
%                   Cohen's d_z, the standard paired effect size.
%     effectLo/Hi - 95% hierarchical-bootstrap percentile CI on effect.
%     propPos     - fraction of values greater than zero. In stats.diffs this
%                   is the common-language effect size, P(stimA > stimB) per
%                   unit; Cliff's delta = 2*propPos - 1 when there are no ties.
%     propPosLo/Hi- 95% hierarchical-bootstrap percentile CI on propPos.
%   Effect-size CIs are always percentile intervals even when ciMethod='sem',
%   because a symmetric +/-SE interval can leave the valid range of a bounded
%   quantity like propPos.
%   Point estimate vs interval: following the standard bootstrap convention
%   (Efron & Tibshirani), the bar is the OBSERVED statistic and the hierarchical
%   resampling supplies only the interval, which is what carries the
%   cross-animal inference. bootMean is reported so the bootstrap bias is
%   visible: it is NOT an equal-weight-per-animal mean, because
%   hierBootMatchFreq pools all resampled observations and takes one flat mean,
%   making the resampled statistic a ratio with a random denominator. Its
%   expectation drifts off the observed mean when neuron counts are unbalanced
%   across animals. A large (mean - bootMean) gap therefore signals strong
%   count imbalance, and a percentile CI may sit noticeably asymmetric about
%   the bar — that asymmetry is informative, not a bug.
%   hierBootMatchFreq returns the mean of each resample, so no bootstrap CI is
%   available for the median; the reported median is a plain sample median.
%   Groups with n < 3 get mean/median but NaN bootstrap columns (no bar drawn).
%
%   Bootstrap details: uses hierBoot (Saravanan et al. 2020) when grouping
%   levels are present, resampling each level with replacement in turn.

% -------------------------------------------------------------------------
% Argument validation block
% -------------------------------------------------------------------------
arguments
    tbl table                                          % observation table
    pairs cell = {}                                    % stim pairs to test
    pValues double = []                                % p-value per pair
    valueField cell = {}                               % field name(s) of value column
    params.nBoot           (1,1) double  = 10000       % bootstrap replicates
    params.fraction        logical       = false       % ratio mode (num/den)
    params.yLegend         char          = 'value'     % y-axis label
    params.diff            logical       = false       % plot per-pair difference only
    params.Xjitter                       = 'density'   % swarm jitter scheme
    params.dotSize                       = 5           % marker size (pts^2)
    params.yMaxVis                       = 1           % visible y-axis cap
    params.filled          logical       = true        % filled vs open markers
    params.Alpha                         = 0.2         % marker face/edge alpha
    params.plotMeanSem     logical       = true        % overlay mean +/- uncertainty
    params.colorByZScore   logical       = false       % color dots by zScore (else by animal)
    params.showBothAndDiff logical       = true        % two-tile raw + diff layout
    params.drawLines       logical       = false       % connect paired observations
    params.rngSeed         (1,1) double  = 0           % bootstrap reproducibility
    params.ciMethod        char          = 'percentile'% 'sem' | 'percentile'
    params.bootGroupVars   cell          = {'__auto__'}% hierarchical bootstrap levels
    params.markUnits       cell          = {}          % N×4 {NeurID,stimulus,animal,insertion}; all marked with the same symbol
    params.nearSigThresh   (1,1) double  = 0.1         % raw panel: draw a bracket only when p < this (exact p is always the label)
    params.pValRes                       = []          % smallest resolvable p behind pValues (2/nBoot of the test that produced them).
                                                       %   Sets the '<' floor printed when a p is exactly 0. Empty -> conservative 0.001.
    params.swarmWidth      (1,1) double  = 0.6         % swarm jitter width; mean/CI bar drawn wider so it always shows
end

% -------------------------------------------------------------------------
% Up-front input validation
% -------------------------------------------------------------------------

% Either raw mode (1 field) or fraction mode (2 fields)
if params.fraction
    assert(numel(valueField) == 2, 'Fraction mode requires two valueField entries.');
else
    assert(~isempty(valueField), 'valueField must contain at least one column name.');
end

% colorByZScore requires the column to exist; downgrade with warning if absent
if params.colorByZScore && ~ismember('zScore', tbl.Properties.VariableNames)
    warning('colorByZScore=true but tbl has no zScore column; falling back to animal coloring.');
    params.colorByZScore = false;
end

% showBothAndDiff places diff in its own tile; it overrides params.diff
if params.showBothAndDiff && params.diff
    warning('showBothAndDiff=true overrides params.diff; diff appears in the right tile only.');
    params.diff = false;
end

% Seed RNG once so bootstraps and dot-draw orders are deterministic
rng(params.rngSeed);

% -------------------------------------------------------------------------
% Resolve bootstrap grouping variables
% -------------------------------------------------------------------------
if isscalar(params.bootGroupVars) && strcmp(params.bootGroupVars{1}, '__auto__')
    cands                = {'animal','insertion'};                          % candidate hierarchy columns
    params.bootGroupVars = cands(ismember(cands, tbl.Properties.VariableNames));
else
    missing = ~ismember(params.bootGroupVars, tbl.Properties.VariableNames);
    assert(~any(missing), 'bootGroupVars contains missing columns: %s', ...
        strjoin(params.bootGroupVars(missing), ', '));
end

% -------------------------------------------------------------------------
% Detect data granularity
% -------------------------------------------------------------------------
% If every (insertion, stimulus) pair appears at most once, the table is
% insertion-level; otherwise neuron-level.
isInsertionLevel = height(tbl) == height(unique(tbl(:, {'insertion','stimulus'})));

% -------------------------------------------------------------------------
% Padding / spacing constants derived from the y-axis cap
% -------------------------------------------------------------------------
yMaxVis    = params.yMaxVis;
bracketPad = yMaxVis * 0.05;                           % gap between data and first bracket
stackPad   = yMaxVis * 0.05;                           % vertical stacking between brackets

% -------------------------------------------------------------------------
% Pre-process tbl: rename legacy labels, reorder categories, compute value
% -------------------------------------------------------------------------
tbl   = renameStimulusLabels(tbl);                     % RG->SB, SDGs->SG, SDGm->MG
pairs = renamePairLabels(pairs);                       % apply same rename to pair labels
tbl   = reorderStimulusByLevel(tbl);                   % sort categories by trailing number

% Normalise the marker specs and apply the SAME label renaming as the table,
% so a user-supplied 'RG' lands in the renamed 'SB' column.
params.markUnitsResolved = resolveMarkUnits(params.markUnits);

if params.fraction
    tbl.value = tbl.(valueField{1}) ./ tbl.(valueField{2});  % element-wise ratio
else
    tbl.value = tbl.(valueField{1});                   % raw value
end

% Drop unused categorical levels so colormaps and counts are accurate
tbl.stimulus  = removecats(tbl.stimulus);
tbl.animal    = removecats(tbl.animal);
tbl.insertion = removecats(tbl.insertion);

% -------------------------------------------------------------------------
% Build figure: single axes or 1x2 tiledlayout
% -------------------------------------------------------------------------
fig = figure;
set(fig, 'Color', 'w');                               % white background for publication

% Summary accumulators. Each panel returns the stats it computed; whatever a
% given plotting mode does not produce is filled in afterwards so the returned
% struct has the same schema in every mode.
statsGroups  = statsSchema();                         % per-stimulus summary (raw panel)
statsDiffOne = statsSchema();                         % single-pair summary (diff panel)

if params.showBothAndDiff
    % Left tile: every stimulus shown raw; right tile: most-significant diff
    tl     = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    axRaw  = nexttile(tl, 1);
    axDiff = nexttile(tl, 2);

    [randiColors, statsGroups] = plotRawSwarm(axRaw, tbl, pairs, pValues, params, ...
        yMaxVis, bracketPad, stackPad, isInsertionLevel);

    % Pick the most significant pair for the diff tile
    if ~isempty(pValues)
        [~, sigIdx] = min(pValues);
    else
        sigIdx = 1;
    end
    pairForDiff = pairs(sigIdx, :);
    pValForDiff = pValues(sigIdx);

    tblDiff = buildDiffTable(tbl, pairForDiff, params, isInsertionLevel);
    [~, statsDiffOne] = plotDiffSwarm(axDiff, tblDiff, pairForDiff, pValForDiff, params, ...
        yMaxVis, bracketPad);
else
    % Single-axes mode: either the raw swarm or the difference
    ax = axes(fig);
    hold(ax, 'on');
    set(ax, 'Clipping', 'off');

    if params.diff
        tblDiff     = buildDiffTable(tbl, pairs, params, isInsertionLevel);
        [randiColors, statsDiffOne] = plotDiffSwarm(ax, tblDiff, pairs, pValues, params, ...
            yMaxVis, bracketPad);
    else
        [randiColors, statsGroups] = plotRawSwarm(ax, tbl, pairs, pValues, params, ...
            yMaxVis, bracketPad, stackPad, isInsertionLevel);
    end
end

% -------------------------------------------------------------------------
% Additional figure: one tile per pairwise difference (only when >1 pair)
% -------------------------------------------------------------------------
if size(pairs, 1) > 1
    [figAllDiffs, statsDiffs] = plotAllPairDiffs(tbl, pairs, pValues, params, ...
        isInsertionLevel, yMaxVis, bracketPad);
else
    figAllDiffs = [];
    statsDiffs  = statsDiffOne;                       % 0 or 1 pair: whatever the diff panel produced
end

% -------------------------------------------------------------------------
% Assemble the summary struct
% -------------------------------------------------------------------------

% Pure-difference mode never draws the raw panel, so the per-stimulus summary
% is computed here instead (same routine, so identical numbers either way).
if isempty(statsGroups)
    statsGroups = computeSwarmStats(tbl, categories(tbl.stimulus), params);
end

% Single pair with no difference panel drawn (raw-only mode): build it now so
% the caller still gets the pairwise difference summary.
if size(pairs, 1) == 1 && isempty(statsDiffs)
    tblDiff1   = buildDiffTable(tbl, pairs, params, isInsertionLevel);
    statsDiffs = computeSwarmStats(tblDiff1, categories(tblDiff1.stimulus), params);
end

% Attach the pair identity and its p-value. plotAllPairDiffs emits exactly one
% row per pair, in pair order (empty tiles emit an all-NaN placeholder row), so
% this is a positional join guarded by a hard assertion.
if ~isempty(pairs) && ~isempty(statsDiffs)
    assert(height(statsDiffs) == size(pairs, 1), ...
        ['plotSwarmBootstrapWithComparisons: %d difference-summary rows for %d ', ...
         'pairs; row alignment broken.'], height(statsDiffs), size(pairs, 1));
    statsDiffs.stimA = string(pairs(:, 1));            % first stimulus of the pair
    statsDiffs.stimB = string(pairs(:, 2));            % second stimulus (subtracted)
    if numel(pValues) == size(pairs, 1)
        statsDiffs.p = pValues(:);                     % caller-supplied p-value per pair
    end
end

stats = struct('groups', statsGroups, 'diffs', statsDiffs);

end % main function




% =========================================================================
% LOCAL FUNCTION: plotRawSwarm
% =========================================================================
function [randiColors, statsGroups] = plotRawSwarm(ax, tbl, pairs, pValues, params, ...
    yMaxVis, bracketPad, stackPad, isInsertionLevel)

hold(ax, 'on');
set(ax, 'Clipping', 'off');                           % brackets may sit above yMaxVis

stimuli = categories(tbl.stimulus);                    % ordered category list
tblPlot = tbl;

% Random permutation for dot draw order (seeded in main)
randiColors = randperm(height(tblPlot));

% Choose dot color source
if params.colorByZScore
    colorData = tblPlot.zScore(randiColors);
else
    colorData = tblPlot.animal(randiColors);
end

% Draw swarm
if params.filled
    s = swarmchart(ax, tblPlot.stimulus(randiColors), tblPlot.value(randiColors), ...
        params.dotSize, colorData, 'filled', ...
        'MarkerFaceAlpha', params.Alpha);
else
    s = swarmchart(ax, tblPlot.stimulus(randiColors), tblPlot.value(randiColors), ...
        params.dotSize, colorData, ...
        'MarkerEdgeAlpha', params.Alpha, 'LineWidth', 1, 'SizeData', 30);
end
s.XJitter = params.Xjitter;
s.XJitterWidth = params.swarmWidth;                    % bound dot spread so the mean/CI bar can always exceed it

% Build short tick labels from encoded category names
Str = string(stimuli);
out = buildTickLabels(Str);

% Apply custom tick labels only if categories contain embedded numbers
numStr = regexp(Str, '-?\d+\.?\d*', 'match', 'once');
hold(ax, 'on');
if any(~cellfun(@isempty, numStr))
    xticklabels(ax, out);
end

% Configure colormap
if params.colorByZScore
    colormap(ax, buildRdBuColormap(256));
    maxZ = max(abs(tblPlot.zScore), [], 'omitnan');
    if isempty(maxZ) || maxZ == 0, maxZ = 1; end
    clim(ax, [-maxZ maxZ]);
    cb = colorbar(ax);
    cb.Label.String = 'Z-score';
else
    colormap(ax, lines(numel(categories(tblPlot.animal))));
end

% -------------------------------------------------------------------------
% Optional connecting lines between paired observations
% -------------------------------------------------------------------------
if params.drawLines && numel(stimuli) <= 2
    if isInsertionLevel
        unitIDvar = 'insertion';                       % insertion IS the unit
    elseif ismember('NeurID', tblPlot.Properties.VariableNames)
        unitIDvar = 'NeurID';                          % neuron is the unit
    else
        unitIDvar = '';
        warning(['drawLines=true on neuron-level data without NeurID; ', ...
                 'skipping connecting lines.']);
    end

    if ~isempty(unitIDvar)
        cats = categories(tblPlot.stimulus);
        xMap = containers.Map(cats, 1:numel(cats));    % stimulus -> x-position
        xNum = arrayfun(@(i) xMap(char(tblPlot.stimulus(i))), 1:height(tblPlot));

        unitIDs = unique(tblPlot.(unitIDvar));
        for u = 1:numel(unitIDs)
            idx = tblPlot.(unitIDvar) == unitIDs(u);
            if nnz(idx) < 2, continue; end
            line(ax, xNum(idx), tblPlot.value(idx), ...
                'Color', [0 0 0 0.1], 'LineWidth', 0.1);
        end
    end
end

ylabel(ax, params.yLegend);
ax.Box   = 'off';
ax.Layer = 'top';

% Hierarchical bootstrap mean +/- SE (or 95% CI). The numbers are computed
% unconditionally so the caller always receives the summary table; only the
% drawing is gated on plotMeanSem.
statsGroups = computeSwarmStats(tblPlot, stimuli, params);   % one row per stimulus, in category order
if params.plotMeanSem
    drawMeanSemBars(ax, statsGroups, params);
end

% Significance brackets. When >5 pairs, only bracket adjacent groups to
% avoid visual clutter; the full set is in figAllDiffs.
maxBracketY = yMaxVis;                                 % default top when no brackets are drawn
if ~isempty(pairs) && numel(pValues) == size(pairs, 1)
    adjacentOnly = size(pairs, 1) > 6;
    maxBracketY = plotBrackets(ax, tblPlot, stimuli, pairs, pValues, ...
        yMaxVis, bracketPad, stackPad, adjacentOnly, params.nearSigThresh, ...
        params.pValRes);
end

% Overlay user-specified example-unit markers (same symbol for every unit)
overlayMarkedUnits(ax, tblPlot, params, true);

% Cap the visible y-range, but extend the top to fit lifted bracket labels so
% they aren't clipped (text near the top of a tiled layout can be cut even
% though the brackets themselves use Clipping=off).
ylim(ax, [ax.YLim(1), max(yMaxVis, maxBracketY)]);

end % plotRawSwarm


% =========================================================================
% LOCAL FUNCTION: plotDiffSwarm
% Draws a swarm plot of pairwise differences for one condition pair,
% overlays a bootstrap mean ± SE bar, a zero reference line, and
% (if significant) a star annotation.
% =========================================================================
function [randiColors, statsDiff] = plotDiffSwarm(ax, tblDiff, pairs, pValues, params, ...
    yMaxVis, bracketPad)

% Hold axes so subsequent drawing commands add to the same panel
hold(ax, 'on');

% Allow markers and text to render outside the strict axes box;
% required so significance stars above yMaxVis remain visible
set(ax, 'Clipping', 'off');

% ---- Guard: empty diff table --------------------------------------------
% This should already be caught by the caller, but defend here too
if height(tblDiff) == 0
    randiColors = [];
    statsDiff   = emptyStatsRow(sprintf('%s-%s', pairs{1,1}, pairs{1,2}));  % NaN placeholder keeps row/pair alignment
    text(ax, 0.5, 0.5, 'No paired data', ...
        'Units',              'normalized', ...
        'HorizontalAlignment','center', ...
        'FontSize',           8, ...
        'Color',              [0.6 0.6 0.6]);
    return
end

% ---- Reproducible draw order --------------------------------------------
% Fix the RNG seed immediately before randperm so that the colour-to-unit
% mapping is identical on every re-render (required for publication).
% Use a deterministic seed derived from the pair names so different panels
% get different (but stable) permutations.
rng(sum(double(char(strjoin({pairs{1,1}, pairs{1,2}}, '')))));
randiColors = randperm(height(tblDiff));   % shuffled row index for plotting order

% ---- Colour source ------------------------------------------------------
% Choose per-dot colour based on z-score (diverging) or animal identity
if params.colorByZScore && ismember('zScore', tblDiff.Properties.VariableNames)
    % Numeric z-score maps onto a diverging colormap set below
    colorData = tblDiff.zScore(randiColors);
else
    % Animal identity is categorical; convert to numeric so swarmchart
    % accepts it as a valid colour index vector
    colorData = double(tblDiff.animal(randiColors));
end

% ---- Swarm chart ---------------------------------------------------------
if params.filled
    % Filled circles; face transparency controlled by params.Alpha
    s = swarmchart(ax, ...
        tblDiff.stimulus(randiColors), ...   % x: condition label (categorical)
        tblDiff.value(randiColors), ...      % y: difference value
        params.dotSize, ...                  % marker area (pt²)
        colorData, ...                       % colour data (numeric)
        'filled', ...
        'MarkerFaceAlpha', params.Alpha);
else
    % Open circles; use params.dotSize for size consistency with filled branch
    s = swarmchart(ax, ...
        tblDiff.stimulus(randiColors), ...
        tblDiff.value(randiColors), ...
        params.dotSize, ...                  % fixed: was hardcoded 30, ignoring params
        colorData, ...
        'MarkerEdgeAlpha', params.Alpha, ...
        'LineWidth',       1);
end

% Horizontal jitter style for the swarm (e.g. 'rand', 'none')
s.XJitter = params.Xjitter;
s.XJitterWidth = params.swarmWidth;                    % bound dot spread so the mean/CI bar can always exceed it

% ---- X tick label --------------------------------------------------------
% Combine the two condition names into a single readable label,
% e.g. "MB 1.57 - MG 90"
pairLabels = buildTickLabels(string({pairs{1,1}, pairs{1,2}}));
xticklabels(ax, join(pairLabels, " - "));

% ---- Colormap & colour bar -----------------------------------------------
if params.colorByZScore && ismember('zScore', tblDiff.Properties.VariableNames)
    % Symmetric diverging colormap centred on zero
    colormap(ax, buildRdBuColormap(256));
    maxZ = max(abs(tblDiff.zScore), [], 'omitnan');   % symmetric colour limit
    if isempty(maxZ) || maxZ == 0, maxZ = 1; end      % avoid degenerate clim
    clim(ax, [-maxZ maxZ]);
    cb = colorbar(ax);
    cb.Label.String = 'Z-score';
else
    % One colour per unique animal; categories() requires a categorical column
    colormap(ax, lines(numel(categories(tblDiff.animal))));
end

% ---- Reference line at zero ---------------------------------------------
% Makes it immediately clear which differences are positive vs negative
yline(ax, 0, 'LineWidth', 1, 'Alpha', 0.7);

% Y-axis label (e.g. 'Δ firing rate (spk/s)')
ylabel(ax, params.yLegend);

% Clean axes appearance for publication
ax.Box   = 'off';
ax.Layer = 'top';   % draw axes on top of data so tick marks are not obscured

% ---- Bootstrap mean ± SE bars -------------------------------------------
% Stats are computed unconditionally (the caller collects them); only drawing
% is gated on plotMeanSem.
stimuli   = categories(tblDiff.stimulus);        % ordered condition list (single 'A-B' entry)
statsDiff = computeSwarmStats(tblDiff, stimuli, params);
if params.plotMeanSem
    drawMeanSemBars(ax, statsDiff, params);
end

% ---- Significance annotation --------------------------------------------
% Capture the current auto-scaled Y limits AFTER all data have been drawn
% so the lower bound reflects the true data range
ylims     = ylim(ax);
topNeeded = yMaxVis;   % how high the axis must reach to fit the label (default: the cap)

if ~isempty(pValues) && numel(pValues) >= 1
    pVal = pValues(1);   % scalar p-value for this difference panel

    % Highest data point within the visible window (anchors the label just
    % above the top of the visible data, not above data clipped by yMaxVis).
    vals       = tblDiff.value;
    maxVisible = min(max(vals(:), [], 'omitnan'), yMaxVis);
    if isempty(maxVisible) || isnan(maxVisible)
        maxVisible = yMaxVis;
    end

    % Label sits bracketPad above the highest visible point
    yText = maxVisible + bracketPad;

    % The exact p-value is always printed on the single difference panel (no
    % stars): thresholded symbols flip between categories on Monte Carlo noise
    % alone when p sits near 0.001 or 0.05. NaN p prints nothing.
    [pNum, pOp] = formatPValue(pVal, params.pValRes);            % plain decimal, no exponent
    if ~isempty(pNum)
        txt = ['p' pOp pNum];                                    % 'p=0.012' or 'p<0.001'
    else
        txt = '';                                                % NaN p: draw nothing
    end

    if ~isempty(txt)
        % Bottom-anchored so the label grows upward, clear of the data
        text(ax, 1, yText, txt, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'bottom', ...
            'FontSize',            7, ...
            'Clipping',            'off');
        topNeeded = yText + bracketPad;   % headroom above the label
    end
end


% Cap the upper Y limit, extending it to fit the significance label so it isn't
% clipped. The lower bound stays auto-scaled (per-panel data range); for the
% all-pairs figure, plotAllPairDiffs takes the union across panels via linkaxes.
ylim(ax, [ylims(1), max(yMaxVis, topNeeded)]);


end % plotDiffSwarm


% =========================================================================
% LOCAL FUNCTION: buildDiffTable
% Per-unit (stimA - stimB) within each insertion.
% Pairing strategy depends on data granularity:
%   insertion-level: direct subtraction (one row per insertion)
%   neuron-level + NeurID: match by NeurID (intersect)
%   neuron-level without NeurID: row-order fallback with warning
% =========================================================================
function tblDiff = buildDiffTable(tbl, pairs, params, isInsertionLevel)

assert(~isempty(pairs) && size(pairs, 1) >= 1, ...
    'diff mode requires at least one stimulus pair.');

stimA = strtrim(pairs{1,1});                           % trim whitespace for safety
stimB = strtrim(pairs{1,2});

hasNeurID = ismember('NeurID', tbl.Properties.VariableNames);

if ~hasNeurID && ~isInsertionLevel
    warning(['buildDiffTable: NeurID column absent for neuron-level data. ', ...
             'Pairing by row order — fragile if rows are reordered.']);
end

ins      = categories(tbl.insertion);                  % unique insertion labels
diffVals = [];                                         % accumulator: paired differences
animals  = categorical.empty(0, 1);                    % accumulator: animal per diff row
insers   = categorical.empty(0, 1);                    % accumulator: insertion per diff row
zScores  = [];                                         % accumulator: zScore (if colorByZScore)
neurIDs  = [];                                          % accumulator: NeurID per diff row (NaN where unavailable)
hasRealIns = ismember('realInsertion', tbl.Properties.VariableNames);
realIns    = categorical.empty(0, 1);                  % accumulator: real insertion per diff row
useZ = params.colorByZScore && ismember('zScore', tbl.Properties.VariableNames);

for i = 1:numel(ins)
    idxA = tbl.insertion == ins{i} & tbl.stimulus == stimA;
    idxB = tbl.insertion == ins{i} & tbl.stimulus == stimB;

    % Skip insertions where either stimulus is absent
    if ~any(idxA) || ~any(idxB)
        continue
    end

    if isInsertionLevel
        % One row per side; direct subtraction
        vA = tbl.value(idxA);
        vB = tbl.value(idxB);
        an = tbl.animal(idxA);
        insCat = tbl.insertion(idxA);                  % preserve original categorical
        nidPiece = NaN;                                % insertion-level: no neuron identity
        if useZ
            zPair = (tbl.zScore(idxA) + tbl.zScore(idxB)) / 2;
        end

    elseif hasNeurID
        % Neuron-level with explicit IDs: safest matching via intersect
        tA = tbl(idxA, :);
        tB = tbl(idxB, :);
        [~, iA, iB] = intersect(tA.NeurID, tB.NeurID, 'stable');
        if isempty(iA)
            continue
        end

        vA = tA.value(iA);
        vB = tB.value(iB);
        an = tA.animal(iA);
        insCat = tA.insertion(iA);
        nidPiece = tA.NeurID(iA);                      % keep the IDs of the matched neurons
        if useZ
            zPair = (tA.zScore(iA) + tB.zScore(iB)) / 2;
        end

    else
        % Row-order fallback (fragile)
        vA = tbl.value(idxA);
        vB = tbl.value(idxB);
        if numel(vA) ~= numel(vB)
            warning('Insertion %s: %d stimA rows vs %d stimB rows; skipping.', ...
                ins{i}, numel(vA), numel(vB));
            continue
        end
        an     = repmat(tbl.animal(find(idxA, 1)), numel(vA), 1);
        insCat = repmat(tbl.insertion(find(idxA, 1)), numel(vA), 1);
        nidPiece = NaN(numel(vA), 1);                  % row-order fallback: IDs not reliable
        if useZ
            zPair = (tbl.zScore(idxA) + tbl.zScore(idxB)) / 2;
        end
    end

    diffVals = [diffVals; vA - vB];                                        %#ok<AGROW>
    animals  = [animals;  an];                                             %#ok<AGROW>
    insers   = [insers;   insCat];                                         %#ok<AGROW>
    neurIDs  = [neurIDs;  nidPiece(:)];
    if hasRealIns
        riThis  = tbl.realInsertion(find(idxA, 1));        % constant within this insertion
        realIns = [realIns; repmat(riThis, numel(vA), 1)]; %#ok<AGROW>
    end
    if useZ
        zScores = [zScores; zPair];                                        %#ok<AGROW>
    end
end

% Drop NaN differences (e.g. from zero-denominator fractions)
valid    = ~isnan(diffVals);
stimName = sprintf('%s-%s', stimA, stimB);

tblDiff           = table();
tblDiff.insertion = insers(valid);                     % categorical insertion labels
tblDiff.stimulus  = categorical(repmat({stimName}, sum(valid), 1));
tblDiff.animal    = animals(valid);
tblDiff.value     = diffVals(valid);
tblDiff.NeurID    = neurIDs(valid);
if hasRealIns
    tblDiff.realInsertion = realIns(valid);
end

if useZ
    tblDiff.zScore = zScores(valid);
end

end % buildDiffTable


% -------------------------------------------------------------------------
% statsSchema, emptyStatsRow, computeSwarmStats and hierBootEffect used to be
% local functions here. They are now public functions in "general functions/"
% so that analyses which do not plot a swarm (e.g. the depth and
% stripe-direction comparisons in analyzeStripeNeurons) can produce
% schema-compatible stats blocks from the same code, and so buildStatsTable can
% compute statistics without creating a figure just to harvest them.
% Behaviour is unchanged; this file calls them exactly as before.
% -------------------------------------------------------------------------



% =========================================================================
% LOCAL FUNCTION: drawMeanSemBars
% Draw the bootstrap central tendency and uncertainty computed by
% computeSwarmStats. Row i of statsTbl corresponds to x-position i, which is
% the invariant computeSwarmStats guarantees by emitting one row per category.
% =========================================================================
function drawMeanSemBars(ax, statsTbl, params)

for i = 1:height(statsTbl)

    % Rows with no bootstrap (absent group, or n < 3) get no bar
    if isnan(statsTbl.bootMean(i))
        continue
    end

    % Point estimate = the OBSERVED mean of the plotted points, so the bar
    % always sits at the centre of its own swarm. The bootstrap contributes
    % only the interval (see computeSwarmStats).
    mu  = statsTbl.mean(i);                            % observed pooled mean (drawn)
    yLo = statsTbl.ciLo(i);                            % lower uncertainty bound
    yHi = statsTbl.ciHi(i);                            % upper uncertainty bound

    xb    = i;                                         % centre the bar on the swarm column
    halfW = params.swarmWidth/2 + 0.2;                % extend past the (now bounded) dot cloud so the bar always shows

    % SEM/CI as a light-gray patch spanning [yLo yHi], centred on the column.
    pat = patch(ax, xb + halfW*[-1 1 1 -1], [yLo yLo yHi yHi], ...  % rectangle covering the uncertainty range
        [0.85 0.85 0.85], 'EdgeColor', 'none');        % light gray fill, no border

    % Mean as a wide, thick horizontal bar at mu.
    mLine = line(ax, xb + halfW*[-1 1], [mu mu], ...   % full-width line at the mean
        'Color', 'k', 'LineWidth', 2.5);               % thicker so it reads clearly over the dots

    % Send both behind the swarm dots (patch furthest back, then the mean line).
    uistack(mLine, 'bottom');                          % mean line drops just below the dots
    uistack(pat,   'bottom');                          % patch drops below the mean line (backmost)
end

end % drawMeanSemBars


% =========================================================================
% LOCAL FUNCTION: plotBrackets
% Pairwise significance brackets. When adjacentOnly is true, only brackets
% between groups at positions i and i+1 are drawn (prevents visual clutter
% with many comparisons). All pairs are always reported to plotAllPairDiffs.
% =========================================================================
function maxY  = plotBrackets(ax, tblPlot, stimuli, pairs, pValues, ...
    yMaxVis, bracketPad, stackPad, adjacentOnly,nearSigThresh, pValRes)

% Track y-positions of placed brackets to prevent overlap
usedHeights = zeros(size(pairs, 1), 1);
maxY = yMaxVis;   % highest point drawn, so the caller can extend ylim to fit labels
for k = 1:size(pairs, 1)

    % Skip NaN and clearly non-significant pairs. Pairs with
    % 0.05 <= p < nearSigThresh are kept and get the exact p above the bracket.
    if isnan(pValues(k)) || pValues(k) >= nearSigThresh
        continue
    end

    % Find x-positions for both stimuli in this pair
    x1 = find(strcmp(stimuli, pairs{k,1}));
    x2 = find(strcmp(stimuli, pairs{k,2}));
    if isempty(x1) || isempty(x2), continue; end

    % --- ADJACENT-ONLY FILTER ---
    % When adjacentOnly is true, skip any pair where the two groups are
    % not next to each other on the x-axis.  This prevents 22+ overlapping
    % brackets when there are many significant comparisons.  The full set
    % of pairwise differences is shown in the separate figAllDiffs figure.
    if adjacentOnly && abs(x1 - x2) > 1
        continue
    end

    % Cap individual values at yMaxVis so the bracket sits at the visible edge
    vals1 = tblPlot.value(tblPlot.stimulus == pairs{k,1});
    vals2 = tblPlot.value(tblPlot.stimulus == pairs{k,2});
    maxVisible = max(min([vals1; vals2], yMaxVis));
    yBase      = maxVisible + bracketPad;

    % Vertical stacking: nudge up if a previous bracket is too close
    y = yBase;
    while any(abs(usedHeights(1:k-1) - y) < stackPad)
        y = y + stackPad;
    end
    usedHeights(k) = y;

    % Horizontal bracket + two short vertical ticks
    line(ax, [x1 x2], [y y],                 'Color', 'k', 'LineWidth', 1.2, 'Clipping', 'off');
    line(ax, [x1 x1], [y - yMaxVis*0.01, y], 'Color', 'k', 'LineWidth', 1.2, 'Clipping', 'off');
    line(ax, [x2 x2], [y - yMaxVis*0.01, y], 'Color', 'k', 'LineWidth', 1.2, 'Clipping', 'off');

    % Exact p above every bracket that is drawn (no stars). Which brackets get
    % drawn is still filtered by nearSigThresh above, so the panel does not fill
    % with labels for clearly null pairs; the complete set of p-values is in the
    % returned stats.diffs table and in the all-pairs figure.
    [pNum, pOp] = formatPValue(pValues(k), pValRes);   % plain decimal, no exponent
    txt = ['p' pOp pNum];                              % 'p=0.012' or 'p<0.001'

    % Numeric labels are taller than stars were, so lift them clear of the
    % bracket and anchor by their bottom edge so they grow upward.
    yTxt = y + bracketPad;                             % lift for numeric labels
    text(ax, mean([x1 x2]), yTxt, txt, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...           % text grows up, off the bracket
        'FontSize', 8, 'Clipping', 'off');

    maxY = max(maxY, yTxt + bracketPad);               % headroom above the lifted label
end

end % plotBrackets


% =========================================================================
% LOCAL FUNCTION: plotAllPairDiffs
% Stand-alone figure with one tile per pairwise difference.
% =========================================================================
function [figAll, statsDiffs] = plotAllPairDiffs(tbl, pairs, pValues, params, ...
    isInsertionLevel, yMaxVis, bracketPad)

% Total number of pairwise comparisons to tile
nPairs = size(pairs, 1);

% One summary row per pair, appended in pair order (including placeholder rows
% for pairs with no paired data) so the caller can join it positionally.
statsDiffs = statsSchema();

% Guard: nothing to plot
if nPairs == 0
    figAll = figure;
    return
end

% Create figure with white background for publication
figAll = figure;
set(figAll, 'Color', 'w');

% Determine grid dimensions: cap columns at 7 to keep tiles legible;
% add as many rows as needed to accommodate all pairs
nCols = min(nPairs, 7);
nRows = ceil(nPairs / nCols);

% Build the tiled layout with tight spacing for a compact multi-panel figure
tl = tiledlayout(figAll, nRows, nCols, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

% Overall figure title
title(tl, 'All pairwise differences');

% Pre-allocate axes handle array so we can link Y axes after the loop.
% Initialise with gobjects so empty-data tiles can be excluded cleanly.
axHandles = gobjects(nPairs, 1);

for k = 1:nPairs

    % Open the next tile in reading order
    ax = nexttile(tl);

    % Extract the condition pair (1×2 cell) and its omnibus p-value for tile k
    pairK  = pairs(k, :);
    pValK  = pValues(k);

    % Build the difference table for this pair (one row per unit/insertion)
    tblDiff = buildDiffTable(tbl, pairK, params, isInsertionLevel);

    % --- Empty-data path ---------------------------------------------------
    if height(tblDiff) == 0
        % Placeholder summary row so statsDiffs keeps one row per pair
        statsDiffs = [statsDiffs; emptyStatsRow(sprintf('%s-%s', pairK{1}, pairK{2}))];  %#ok<AGROW>

        % Suppress all axes decorations so the blank tile is invisible in print
        axis(ax, 'off');

        % Build a human-readable label from the two condition names
        pairLabels = buildTickLabels(string({pairK{1}, pairK{2}}));

        % Place a low-contrast annotation so the missing pair is still legible
        % when inspecting the figure interactively
        text(ax, 0.5, 0.5, join(pairLabels, " - ") + " (no data)", ...
            'Units',              'normalized', ...
            'HorizontalAlignment','center', ...
            'FontSize',           6, ...
            'Color',              [0.6 0.6 0.6]);

        % Do NOT store this handle — we exclude empty tiles from Y-axis linking
        continue
    end
    % -----------------------------------------------------------------------

    % Draw the swarm + significance bracket into ax, keeping its summary row
    [~, statsK] = plotDiffSwarm(ax, tblDiff, pairK, pValK, params, ...
        yMaxVis, bracketPad);

    % Each diff table holds exactly one condition ('A-B'), so one row per pair
    assert(height(statsK) == 1, ...
        'plotAllPairDiffs: pair %d produced %d summary rows; expected 1.', ...
        k, height(statsK));
    statsDiffs = [statsDiffs; statsK];                 %#ok<AGROW>

    % Record the handle for this successfully plotted tile
    axHandles(k) = ax;

end % for k

% --- Shared Y axis ---------------------------------------------------------
% Retain only the axes handles that were actually plotted (non-null gobjects)
validAx = axHandles(isvalid(axHandles));

if numel(validAx) > 1
    % Couple all valid tile axes so interactive zoom/pan stays synchronised
    linkaxes(validAx, 'y');

    % Compute the union of all individual Y ranges so no data are clipped.
    % We do this explicitly rather than relying on linkaxes auto-scaling,
    % because plotDiffSwarm may have expanded ylim for brackets/text and
    % we want the *tightest* common range that honours those expansions.
    allYLims = cell2mat(arrayfun(@(a) ylim(a), validAx, ...
        'UniformOutput', false));      % numel(validAx) × 2 matrix
    sharedYLim = [min(allYLims(:,1)),  max(allYLims(:,2))];

    % Apply the shared limits; linkaxes propagates this to every coupled axis
    ylim(validAx(1), sharedYLim);
end
% ---------------------------------------------------------------------------

end % plotAllPairDiffs


% =========================================================================
% LOCAL FUNCTION: renameStimulusLabels
% =========================================================================
function tbl = renameStimulusLabels(tbl)
s = string(tbl.stimulus);
s = replace(s, "RG",   "SB");
s = replace(s, "SDGs", "SG");
s = replace(s, "SDGm", "MG");
tbl.stimulus = categorical(s);
end


% =========================================================================
% LOCAL FUNCTION: renamePairLabels
% =========================================================================
function pairs = renamePairLabels(pairs)
if isempty(pairs), return; end
for i = 1:numel(pairs)
    p = string(pairs{i});
    p = replace(p, "RG",   "SB");
    p = replace(p, "SDGs", "SG");
    p = replace(p, "SDGm", "MG");
    pairs{i} = char(p);
end
end


% =========================================================================
% LOCAL FUNCTION: buildTickLabels
% Decode category names into short human-readable tick labels.
% e.g. 'MB_dir_1p57' -> 'MB 1.57', 'MG_ang_90' -> 'MG 90'
% =========================================================================
% =========================================================================
function out = buildTickLabels(Str)
% buildTickLabels  Convert coded condition strings (e.g. "MB0p05") into
%   reader-friendly tick labels (e.g. "MB 0.05"), stripping redundant
%   prefixes when only one prefix exists in the set, and silently
%   discarding any entries whose numeric payload is NaN.
%
%   Input:  Str  – string array of coded condition names
%   Output: out  – string array of formatted labels (same size as Str,
%                   with NaN entries replaced by "")

% ---- Step 1: extract the leading uppercase prefix from each entry ----
% Matches one or more capital letters at the start of each string.
% Returns <missing> for entries without an uppercase prefix.
prefixes = strings(size(Str));                    % pre-allocate output
for i = 1:numel(Str)
    prefixes(i) = regexp(Str(i), '^[A-Z]+', 'match', 'once');
end

% ---- Step 2: decide whether the prefix carries information ----
% If every entry shares the same prefix (e.g. all start with "TF"),
% the prefix is redundant and can be dropped to declutter the axis.
uniquePrefixes = unique(prefixes(~ismissing(prefixes)));  % ignore <missing> when counting
removePrefix   = numel(uniquePrefixes) <= 1;              % true → prefixes are uninformative

% ---- Step 3: format each entry ----
out = strings(size(Str));                         % pre-allocate output

for i = 1:numel(Str)
    s      = Str(i);                              % current coded string
    prefix = prefixes(i);                         % its uppercase prefix (or <missing>)

    % -- Extract the numeric part of the string --
    % Matches an optional minus sign, one or more digits, and an optional
    % decimal part delimited by 'p' or '.'.  (Inside a character class
    % the dot is already literal, so no backslash is needed.)
    numStr = regexp(s, '-?\d+(?:[p.]\d+)?', 'match', 'once');

    % -- Guard: skip entries with no numeric payload --
    if ismissing(numStr)
        out(i) = s;                               % preserve non-numeric labels as-is
        continue
    end

    % -- Convert the 'p' decimal convention to a real decimal point --
    numStr = replace(numStr, 'p', '.');

    % -- Convert to a numeric value for formatting decisions --
    numVal = str2double(numStr);

    % ================================================================
    % BUG FIX: upstream code occasionally produces NaN-valued entries
    % (e.g. "MB NaN").  The digit regex above will not match "NaN",
    % but it *will* match if the string is e.g. "MBNaN" with stray
    % digits nearby.  Either way, any NaN that survives to this point
    % must be caught or it will appear as "MB NaN" on the axis.
    % ================================================================
    if isnan(numVal)
        out(i) = "";                              % blank label; caller should
        continue                                  %   remove or handle these
    end

    % -- Choose a readable numeric format --
    % Use scientific notation for very small non-zero magnitudes;
    % otherwise use fixed-point with two decimal places.
    if abs(numVal) < 0.01 && numVal ~= 0
        numFormatted = compose("%.2e", numVal);   % e.g. "1.50e-03"
    else
        numFormatted = compose("%.2f", numVal);   % e.g. "0.25"
    end

    % -- Strip unnecessary trailing zeros from the mantissa --
    % For scientific notation the zeros sit before the 'e', so we
    % first handle the mantissa, then reassemble with the exponent.
    if contains(numFormatted, 'e')
        parts = split(numFormatted, 'e');         % {"1.50", "-03"}
        parts(1) = regexprep(parts(1), '\.?0+$', '');  % "1.50" → "1.5"
        numFormatted = join(parts, 'e');           % "1.5e-03"
    else
        numFormatted = regexprep(numFormatted, '\.?0+$', '');  % "0.250" → "0.25"
    end

    % ---- Reassemble: prefix + number (or number alone) ----
    if ~ismissing(prefix) && ~removePrefix
        out(i) = prefix + " " + numFormatted;     % e.g. "MB 0.25"
    else
        out(i) = numFormatted;                     % e.g. "0.25"
    end
end
end

% =========================================================================
% LOCAL FUNCTION: reorderStimulusByLevel
% Reorder tbl.stimulus categories ascending by the trailing numeric token
% of each label (e.g. 'MB_dir_0' < 'MB_dir_30').  Reverses the encoding
% used by AllExpAnalysis: 'p' -> '.', 'neg' -> '-'.
% No-op if fewer than 2 labels have a numeric trailing token.
% =========================================================================
function tbl = reorderStimulusByLevel(tbl)

cats = categories(tbl.stimulus);
nums = nan(numel(cats), 1);

for i = 1:numel(cats)
    parts = strsplit(cats{i}, '_');
    if numel(parts) < 2, continue; end                 % no underscore => no level token

    last = parts{end};                                 % trailing token
    last = strrep(last, 'p',   '.');                   % decode decimal
    last = strrep(last, 'neg', '-');                    % decode negative

    v = str2double(last);
    if ~isnan(v), nums(i) = v; end
end

% Need at least 2 numeric tokens to reorder; otherwise leave alphabetical
if sum(~isnan(nums)) < 2
    stimOrder = unique(string(tbl.stimulus), 'stable');
    tbl.stimulus = reordercats(tbl.stimulus, cellstr(stimOrder));
    return
end

% Two-step stable sort: primary numeric ascending, secondary alphabetical
[catsAlpha, idxAlpha] = sort(cats);
numsAlpha             = nums(idxAlpha);
[~, idxNum]           = sort(numsAlpha, 'ascend', 'MissingPlacement', 'last');
catsFinal             = catsAlpha(idxNum);

tbl.stimulus = reordercats(tbl.stimulus, catsFinal);

end % reorderStimulusByLevel


% =========================================================================
% LOCAL FUNCTION: buildRdBuColormap
% =========================================================================
function cmap = buildRdBuColormap(n)
half = floor(n/2);
blueToWhite = [linspace(0.02, 1, half)', ...
               linspace(0.44, 1, half)', ...
               linspace(0.69, 1, half)'];
whiteToRed  = [linspace(1, 0.70, half)', ...
               linspace(1, 0.09, half)', ...
               linspace(1, 0.09, half)'];
cmap = [blueToWhite; whiteToRed];
end

% =========================================================================
% LOCAL FUNCTION: resolveMarkUnits
% Normalise the markUnits cell into a struct array, applying the same stimulus
% renaming the table receives, and assigning a marker symbol per row.
% =========================================================================
function specs = resolveMarkUnits(markUnits)


% Single marker shared by every example unit (same symbol for all)
markerSymbol = 'o';

% Empty / unset input: return a 0x0 struct with the right fields (overlay no-ops)
if isempty(markUnits)
    specs = struct('id', {}, 'stim', {}, 'animal', {}, 'insertion', {}, 'symbol', {});
    return
end

% Shape check: must be an N×4 cell {NeurID, stimulus, animal, insertion}
assert(iscell(markUnits) && size(markUnits, 2) == 4, ...
    'markUnits must be an N×4 cell: {NeurID, stimulus, animal, insertion}.');

n = size(markUnits, 1);                                % number of requested units (all share one symbol; no cap)

% Pre-allocate the output struct array
specs = struct('id', cell(1, n), 'stim', cell(1, n), ...
               'animal', cell(1, n), 'insertion', cell(1, n), 'symbol', cell(1, n));

for i = 1:n
    % NeurID -> double, tolerating numeric or char/categorical input
    rawID = markUnits{i, 1};                           % raw NeurID entry
    if isnumeric(rawID)
        specs(i).id = double(rawID);                   % already numeric
    else
        specs(i).id = str2double(string(rawID));       % char/categorical -> numeric (not double())
    end

    % Stimulus: apply the SAME replacements as renameStimulusLabels, same order
    stimStr = string(markUnits{i, 2});                 % raw stimulus label
    stimStr = replace(stimStr, "RG",   "SB");          % moving/static ball grid rename
    stimStr = replace(stimStr, "SDGs", "SG");          % static grating rename
    stimStr = replace(stimStr, "SDGm", "MG");          % moving grating rename
    specs(i).stim = char(stimStr);                     % store as char for equality tests
    
    % Animal / insertion: store as char so string() comparison is unambiguous
    specs(i).animal    = char(string(markUnits{i, 3}));         % animal label (matched as string)
    specs(i).insertion = str2double(string(markUnits{i, 4}));   % REAL insertion number (recordingName trailing token)

    % Same symbol for every unit
    specs(i).symbol = markerSymbol;                    % all example units use the same marker
end

end % resolveMarkUnits


% =========================================================================
% LOCAL FUNCTION: overlayMarkedUnits
% Draw the example-unit symbols onto an already-rendered swarm. A neuron is
% identified only by NeurID x animal x insertion (NeurID repeats across both),
% plus stimulus in the raw panels. Markers sit at the category centre and the
% unit's exact value; the jittered x is not retrievable on a categorical ruler.
% =========================================================================
function overlayMarkedUnits(ax, tblPlot, params, matchByStimulus)

% Pull the resolved specs; nothing to do if none were given
specs = params.markUnitsResolved;                      % struct array (possibly empty)
if isempty(specs)
    return
end

% Marking is keyed on NeurID, so the column must exist
if ~ismember('NeurID', tblPlot.Properties.VariableNames)
    warning('overlayMarkedUnits: tbl has no NeurID column; cannot mark units.');
    return
end

if ~ismember('realInsertion', tblPlot.Properties.VariableNames)
    warning('overlayMarkedUnits: tbl has no realInsertion column; cannot mark units.');
    return
end

% Numeric NeurID vector (categorical loaded from .mat needs str2double(string()))
if iscategorical(tblPlot.NeurID)
    nid = str2double(string(tblPlot.NeurID));          % category labels -> numeric
else
    nid = double(tblPlot.NeurID);                      % already numeric
end

% String views of the grouping columns for robust equality
animalStr = string(tblPlot.animal);                    % per-row animal
stimStr   = string(tblPlot.stimulus);                  % per-row stimulus (renamed already)
realInsTbl = str2double(string(tblPlot.realInsertion));% per-row REAL insertion number

% Loop over the requested example units
for m = 1:numel(specs)
    spec = specs(m);                                   % current spec

    % Rows for the requested animal
    isAnimal = animalStr == spec.animal;
    if ~any(isAnimal)
        warning('overlayMarkedUnits: no rows for animal %s; marker skipped.', spec.animal);
        continue
    end

    % Rows for the requested animal
    isAnimal = animalStr == spec.animal;
    if ~any(isAnimal)
        warning('overlayMarkedUnits: no rows for animal %s; marker skipped.', spec.animal);
        continue
    end

    % Match by the REAL insertion number. The number you type IS the recording's
    % own insertion label, so this is a direct value match — no offset, gap-safe,
    % and unaffected by a dropped first insertion.
    if ~any(realInsTbl(isAnimal) == spec.insertion)
        warning(['overlayMarkedUnits: insertion %g is not present in animal %s ', ...
            '(likely no somatic/responsive units there); marker skipped.'], ...
            spec.insertion, spec.animal);
        continue
    end

    % Identity mask: animal + real insertion number + neuron
    mask = isAnimal ...
        & realInsTbl == spec.insertion ...            % real insertion, by value
        & nid        == spec.id;
    if matchByStimulus
        mask = mask & stimStr == spec.stim;            % raw panel: restrict to the stimulus column
    end

    rows = find(mask);

    % Insertion exists but this neuron/stimulus combination doesn't
    if isempty(rows)
        if matchByStimulus
            warning(['overlayMarkedUnits: insertion %g of animal %s exists, but no ', ...
                     'NeurID %g under stimulus %s; marker skipped.'], ...
                     spec.insertion, spec.animal, spec.id, spec.stim);
        else
            fprintf(['overlayMarkedUnits: NeurID %g (animal %s, insertion %g) ', ...
                     'absent from this difference panel; skipped.\n'], ...
                     spec.id, spec.animal, spec.insertion);
        end
        continue
    end

    if numel(rows) > 1
        warning(['overlayMarkedUnits: %d rows match NeurID %g (animal %s, ', ...
            'insertion %g); marking all of them.'], ...
            numel(rows), spec.id, spec.animal, spec.insertion);
    end

    % Draw each matched unit (normally exactly one)
    for r = rows'
        xMark = tblPlot.stimulus(r);                   % categorical category centre (axis-safe)
        yMark = tblPlot.value(r);                      % the unit's exact value

        % White halo first, then symbol on top; both small and at alpha 0.6 so
        % the underlying swarm dot stays visible through the mark.
        % plot(ax, xMark, yMark, spec.symbol, ...
        %     'Color', [1 1 1 0.6], 'MarkerSize', 6, 'LineWidth', 2.0, 'Clipping', 'off');  % halo
        plot(ax, xMark, yMark, spec.symbol, ...
            'Color', [0 0 0 0.6], 'MarkerSize', 4, 'LineWidth', 1.2, 'Clipping', 'off');  % symbol
    end
end

end % overlayMarkedUnits