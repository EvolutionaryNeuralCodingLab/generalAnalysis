function [fig, randiColors, figAllDiffs] = plotSwarmBootstrapWithComparisons(tbl, pairs, pValues, valueField, params)
% PLOTSWARMBOOTSTRAPWITHCOMPARISONS  Per-stimulus swarm plot with hierarchical
%   bootstrap central tendency, uncertainty bar, and pairwise significance brackets.
%
%   [fig, randiColors, figAllDiffs] = plotSwarmBootstrapWithComparisons(tbl, ...
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
textPad    = yMaxVis * 0.01;                           % gap between bracket and star text
semAlpha   = 0.6;                                      % alpha for bootstrap error bars

% -------------------------------------------------------------------------
% Pre-process tbl: rename legacy labels, reorder categories, compute value
% -------------------------------------------------------------------------
tbl   = renameStimulusLabels(tbl);                     % RG->SB, SDGs->SG, SDGm->MG
pairs = renamePairLabels(pairs);                       % apply same rename to pair labels
tbl   = reorderStimulusByLevel(tbl);                   % sort categories by trailing number

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

if params.showBothAndDiff
    % Left tile: every stimulus shown raw; right tile: most-significant diff
    tl     = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    axRaw  = nexttile(tl, 1);
    axDiff = nexttile(tl, 2);

    randiColors = plotRawSwarm(axRaw, tbl, pairs, pValues, params, ...
        yMaxVis, bracketPad, stackPad, textPad, semAlpha, isInsertionLevel);

    % Pick the most significant pair for the diff tile
    if ~isempty(pValues)
        [~, sigIdx] = min(pValues);
    else
        sigIdx = 1;
    end
    pairForDiff = pairs(sigIdx, :);
    pValForDiff = pValues(sigIdx);

    tblDiff = buildDiffTable(tbl, pairForDiff, params, isInsertionLevel);
    plotDiffSwarm(axDiff, tblDiff, pairForDiff, pValForDiff, params, ...
        yMaxVis, bracketPad, textPad);
else
    % Single-axes mode: either the raw swarm or the difference
    ax = axes(fig);
    hold(ax, 'on');
    set(ax, 'Clipping', 'off');

    if params.diff
        tblDiff     = buildDiffTable(tbl, pairs, params, isInsertionLevel);
        randiColors = plotDiffSwarm(ax, tblDiff, pairs, pValues, params, ...
            yMaxVis, bracketPad, textPad);
    else
        randiColors = plotRawSwarm(ax, tbl, pairs, pValues, params, ...
            yMaxVis, bracketPad, stackPad, textPad, semAlpha, isInsertionLevel);
    end
end

% -------------------------------------------------------------------------
% Additional figure: one tile per pairwise difference (only when >1 pair)
% -------------------------------------------------------------------------
if size(pairs, 1) > 1
    figAllDiffs = plotAllPairDiffs(tbl, pairs, pValues, params, ...
        isInsertionLevel, yMaxVis, bracketPad, textPad);
else
    figAllDiffs = [];
end

end % main function




% =========================================================================
% LOCAL FUNCTION: plotRawSwarm
% =========================================================================
function randiColors = plotRawSwarm(ax, tbl, pairs, pValues, params, ...
    yMaxVis, bracketPad, stackPad, textPad, semAlpha, isInsertionLevel)

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

% Hierarchical bootstrap mean +/- SE (or 95% CI)
if params.plotMeanSem
    plotMeanSemBars(ax, tblPlot, stimuli, params, semAlpha);
end

% Significance brackets. When >5 pairs, only bracket adjacent groups to
% avoid visual clutter; the full set is in figAllDiffs.
if ~isempty(pairs) && numel(pValues) == size(pairs, 1)
    adjacentOnly = size(pairs, 1) > 6;
    plotBrackets(ax, tblPlot, stimuli, pairs, pValues, ...
        yMaxVis, bracketPad, stackPad, textPad, adjacentOnly);
end

% Cap visible y-range (brackets use Clipping=off so they remain visible)
ylim(ax, [ax.YLim(1) yMaxVis]);

end % plotRawSwarm


% =========================================================================
% LOCAL FUNCTION: plotDiffSwarm
% Draws a swarm plot of pairwise differences for one condition pair,
% overlays a bootstrap mean ± SE bar, a zero reference line, and
% (if significant) a star annotation.
% =========================================================================
function randiColors = plotDiffSwarm(ax, tblDiff, pairs, pValues, params, ...
    yMaxVis, bracketPad, textPad)

% Hold axes so subsequent drawing commands add to the same panel
hold(ax, 'on');

% Allow markers and text to render outside the strict axes box;
% required so significance stars above yMaxVis remain visible
set(ax, 'Clipping', 'off');

% ---- Guard: empty diff table --------------------------------------------
% This should already be caught by the caller, but defend here too
if height(tblDiff) == 0
    randiColors = [];
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
if params.plotMeanSem
    stimuli = categories(tblDiff.stimulus);      % ordered condition list
    plotMeanSemBars(ax, tblDiff, stimuli, params, 0.6);
end

% ---- Significance annotation --------------------------------------------
% Capture the current auto-scaled Y limits AFTER all data have been drawn
% so the lower bound reflects the true data range
ylims = ylim(ax);

if ~isempty(pValues) && numel(pValues) >= 1
    pVal = pValues(1);   % scalar p-value for this difference panel

    % Compute the highest data point that falls within the visible window.
    % This anchors the star text just above the top of the visible data,
    % not above data that has been clipped by yMaxVis.
    % Equivalent to min(max(vals), yMaxVis) but written for clarity.
    vals       = tblDiff.value;
    maxVisible = min(max(vals(:), [], 'omitnan'), yMaxVis);  % fixed: was max(min(...)) with no 'omitnan'

    % Fallback if vals was empty or all-NaN (isempty catches []; isnan catches NaN)
    if isempty(maxVisible) || isnan(maxVisible)
        maxVisible = yMaxVis;
    end

    % Text sits bracketPad above the highest visible data point
    yText = maxVisible + bracketPad;

    % Draw significance stars only when the result clears the α = 0.05 threshold
    if ~isnan(pVal) && pVal < 0.05
        if     pVal < 0.001, txt = '***';
        elseif pVal < 0.01,  txt = '**';
        else,                txt = '*';
        end
        % x = 1 centres the label over the single swarm column;
        % 'Clipping','off' keeps it visible even when yText > yMaxVis
        text(ax, 1, yText, txt, ...
            'HorizontalAlignment', 'center', ...
            'FontSize',            7, ...
            'Clipping',            'off');
    end
end

% Fix the upper Y limit to yMaxVis so all panels in plotAllPairDiffs
% start from the same ceiling before linkaxes applies the shared lower bound.
% The lower bound is left as auto-scaled so it reflects each panel's data range;
% plotAllPairDiffs will then extend it to the global minimum via linkaxes.
ylim(ax, [ylims(1), yMaxVis]);

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
        if useZ
            zPair = (tbl.zScore(idxA) + tbl.zScore(idxB)) / 2;
        end
    end

    diffVals = [diffVals; vA - vB];                                        %#ok<AGROW>
    animals  = [animals;  an];                                             %#ok<AGROW>
    insers   = [insers;   insCat];                                         %#ok<AGROW>
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

if useZ
    tblDiff.zScore = zScores(valid);
end

end % buildDiffTable


% =========================================================================
% LOCAL FUNCTION: plotMeanSemBars
% Hierarchical-bootstrap central tendency and uncertainty per stimulus.
% Uses hierBoot (Saravanan et al. 2020) for hierarchical resampling.
% =========================================================================
function plotMeanSemBars(ax, tblPlot, stimuli, params, semAlpha)

for i = 1:numel(stimuli)
    idx = tblPlot.stimulus == stimuli{i};
    if ~any(idx), continue; end

    % Pull values and drop NaNs
    vals = tblPlot.value(idx);
    keep = ~isnan(vals);
    vals = vals(keep);

    n = numel(vals);
    if n < 3
        fprintf('Stimulus %s: n=%d < 3; skipping mean/SE bar.\n', char(stimuli{i}), n);
        continue
    end

    % Pull each grouping column aligned with NaN drop
    groupVars = params.bootGroupVars;
    groupVals = cell(1, numel(groupVars));
    for g = 1:numel(groupVars)
        col = tblPlot.(groupVars{g})(idx);
        col = col(keep);
        if iscategorical(col)
            col = double(col);                         % hierBoot requires numeric
        end
        groupVals{g} = col;
    end

    % Hierarchical or flat bootstrap
    if isempty(groupVars)
        bootMean = bootstrp(params.nBoot, @mean, vals);
    else
        bootMean = hierBoot(vals, params.nBoot, groupVals{:});
    end

    % Point estimate: mean of the bootstrap distribution
    % (weights animals/insertions equally, matching the mixed model)
    mu = mean(bootMean);

    % Uncertainty bar
    switch lower(params.ciMethod)
        case 'sem'
            se  = std(bootMean);
            yLo = mu - se;
            yHi = mu + se;
        case 'percentile'
            yLo = prctile(bootMean,  2.5);
            yHi = prctile(bootMean, 97.5);
        otherwise
            error('Unknown ciMethod: %s. Use ''sem'' or ''percentile''.', params.ciMethod);
    end

    % Vertical uncertainty line
    line(ax, [i i], [yLo yHi], ...
        'Color', [0 0 0 semAlpha], 'LineWidth', 2);
    % End caps
    capW = 0.1;
    line(ax, [i-capW i+capW], [yHi yHi], ...
        'Color', [0 0 0 semAlpha], 'LineWidth', 2);
    line(ax, [i-capW i+capW], [yLo yLo], ...
        'Color', [0 0 0 semAlpha], 'LineWidth', 2);
    % Mean line
    dx = 0.15;
    plot(ax, [i-dx i+dx], [mu mu], 'k-', 'LineWidth', 1.2);
end

end % plotMeanSemBars


% =========================================================================
% LOCAL FUNCTION: plotBrackets
% Pairwise significance brackets. When adjacentOnly is true, only brackets
% between groups at positions i and i+1 are drawn (prevents visual clutter
% with many comparisons). All pairs are always reported to plotAllPairDiffs.
% =========================================================================
function plotBrackets(ax, tblPlot, stimuli, pairs, pValues, ...
    yMaxVis, bracketPad, stackPad, textPad, adjacentOnly)

% Track y-positions of placed brackets to prevent overlap
usedHeights = zeros(size(pairs, 1), 1);

for k = 1:size(pairs, 1)

    % Skip non-significant (no bracket, no text)
    if isnan(pValues(k)) || pValues(k) >= 0.05
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

    % Star annotation
    if     pValues(k) < 0.001, txt = '***';
    elseif pValues(k) < 0.01,  txt = '**';
    elseif pValues(k) < 0.05,  txt = '*';
    end
    text(ax, mean([x1 x2]), y + textPad, txt, ...
        'HorizontalAlignment', 'center', 'FontSize', 7, 'Clipping', 'off');
end

end % plotBrackets


% =========================================================================
% LOCAL FUNCTION: plotAllPairDiffs
% Stand-alone figure with one tile per pairwise difference.
% =========================================================================
function figAll = plotAllPairDiffs(tbl, pairs, pValues, params, ...
    isInsertionLevel, yMaxVis, bracketPad, textPad)

% Total number of pairwise comparisons to tile
nPairs = size(pairs, 1);

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

    % Draw the swarm + significance bracket into ax
    plotDiffSwarm(ax, tblDiff, pairK, pValK, params, ...
        yMaxVis, bracketPad, textPad);

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