function [fig, randiColors] = plotSwarmBootstrapWithComparisons(tbl, pairs, pValues, valueField, params)

arguments
    tbl table
    pairs cell = {}
    pValues double = []
    valueField cell = {}
    params.nBoot        (1,1) double  = 10000
    params.fraction     logical       = false
    params.yLegend      char          = 'value'
    params.diff         logical       = false   % compute difference between first pair
    params.Xjitter                    = 'density'
    params.dotSize                    = 7
    params.yMaxVis                    = 1
    params.filled       logical       = true
    params.Alpha                      = 0.2
    params.plotMeanSem  logical       = true
    params.colorByZScore logical      = false   % color dots by zScore column in tbl instead of by animal
    params.showBothAndDiff logical    = true   % show raw (both stim types) in left tile AND difference in right tile
    params.drawLines logical          = false %draw lines between points
end

% -------------------------------------------------------------------------
% Validate Z-score coloring request
% -------------------------------------------------------------------------
if params.colorByZScore && ~ismember('zScore', tbl.Properties.VariableNames)
    warning('colorByZScore=true but tbl has no zScore column — falling back to animal coloring.');
    params.colorByZScore = false;
end

% showBothAndDiff implicitly requires diff=false at the top level, because
% it manages the diff internally in the right tile
if params.showBothAndDiff && params.diff
    warning('showBothAndDiff=true overrides params.diff — diff will be shown in the right tile only.');
    params.diff = false;
end

% -------------------------------------------------------------------------
% Shared parameters derived from yMaxVis
% -------------------------------------------------------------------------
yMaxVis    = params.yMaxVis;
bracketPad = yMaxVis * 0.05;
stackPad   = yMaxVis * 0.05;
textPad    = yMaxVis * 0.01;
semAlpha   = 0.6;

% -------------------------------------------------------------------------
% Pre-process tbl: rename stim labels and compute value column
% -------------------------------------------------------------------------
tbl = renameStimulusLabels(tbl);           % RG->SB, SDGs->SG, SDGm->MG
pairs = renamePairLabels(pairs);           % same substitution in pairs

if params.fraction
    assert(numel(valueField) == 2, 'Fraction mode requires two valueField entries.');
    tbl.value = tbl.(valueField{1}) ./ tbl.(valueField{2});
else
    tbl.value = tbl.(valueField{1});
end

tbl.stimulus  = removecats(tbl.stimulus);
tbl.animal    = removecats(tbl.animal);
tbl.insertion = removecats(tbl.insertion);

% -------------------------------------------------------------------------
% Build figure: single axes OR tiledlayout(1,2) for showBothAndDiff
% -------------------------------------------------------------------------
fig = figure;
set(fig, 'Color', 'w');

if params.showBothAndDiff
    % Two tiles: [raw both stim types | difference]
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    axRaw  = nexttile(tl, 1);  % left tile: raw swarm
    axDiff = nexttile(tl, 2);  % right tile: difference swarm

    % --- LEFT TILE: raw (both stim types) ---
    randiColors = plotRawSwarm(axRaw, tbl, pairs, pValues, params, ...
        yMaxVis, bracketPad, stackPad, textPad, semAlpha);

    % --- RIGHT TILE: difference ---
    tblDiff = buildDiffTable(tbl, pairs, params);
    plotDiffSwarm(axDiff, tblDiff, pairs, pValues, params, ...
        yMaxVis, bracketPad, textPad);

else
    % Single axes mode
    ax = axes(fig);  %#ok<LAXES>
    hold(ax, 'on');
    set(ax, 'Clipping', 'off');

    if params.diff
        % Difference mode only
        tblDiff     = buildDiffTable(tbl, pairs, params);
        randiColors = plotDiffSwarm(ax, tblDiff, pairs, pValues, params, ...
            yMaxVis, bracketPad, textPad);
    else
        % Raw mode only
        randiColors = plotRawSwarm(ax, tbl, pairs, pValues, params, ...
            yMaxVis, bracketPad, stackPad, textPad, semAlpha);
    end
end

end % main function


% =========================================================================
% LOCAL FUNCTION: plotRawSwarm
% Plots raw values for all stim types with lines between paired neurons.
% Returns randiColors (permuted indices used for dot ordering).
% =========================================================================
function randiColors = plotRawSwarm(ax, tbl, pairs, pValues, params, ...
    yMaxVis, bracketPad, stackPad, textPad, semAlpha)

hold(ax, 'on');
set(ax, 'Clipping', 'off');

stimuli    = categories(tbl.stimulus);
tblPlot    = tbl;

% Randomize dot draw order so overlapping colors are visible
randiColors = randperm(height(tblPlot));

% Determine dot color data: Z-score (diverging) or animal (categorical)
if params.colorByZScore
    % Use Z-score values — colormap set to diverging RdBu below
    colorData = tblPlot.zScore(randiColors);
else
    % Use animal labels — colormap set to lines() below
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

% Set colormap and colorbar
if params.colorByZScore
    % Diverging red-blue colormap centred at 0
    colormap(ax, buildRdBuColormap(256));
    maxZ = max(abs(tblPlot.zScore), [], 'omitnan');
    if maxZ == 0, maxZ = 1; end
    clim(ax, [-maxZ maxZ]);
    cb = colorbar(ax);
    cb.Label.String = 'Z-score';
else
    colormap(ax, lines(numel(categories(tblPlot.animal))));
end

% Draw lines between paired neurons across stim types
if params.drawLines
    cats = categories(tblPlot.stimulus);
    xMap = containers.Map(cats, 1:numel(cats));
    xNum = arrayfun(@(i) xMap(char(tblPlot.stimulus(i))), 1:height(tblPlot));


    try
        tblPlot.NeurID;
        UI = 'NeurID';
    catch
        UI = 'insertion';
    end

    for i = 1:numel(unique(tblPlot.(UI)))
        idx = double(tblPlot.(UI)) == i;
        if sum(idx) < 2, continue; end
        line(ax, xNum(idx), tblPlot.value(idx), ...
            'Color', [0 0 0 0.1], 'LineWidth', 0.1);
    end
end

ylabel(ax, params.yLegend);
ax.Box   = 'off';
ax.Layer = 'top';

% Bootstrap mean + SEM per stimulus
if params.plotMeanSem
    plotMeanSemBars(ax, tblPlot, stimuli, params, semAlpha);
end

% Pairwise significance brackets
if ~isempty(pairs) && numel(pValues) == size(pairs,1)
    plotBrackets(ax, tblPlot, stimuli, pairs, pValues, yMaxVis, bracketPad, stackPad, textPad);
end

ylim(ax, [ax.YLim(1) yMaxVis]);

end % plotRawSwarm


% =========================================================================
% LOCAL FUNCTION: plotDiffSwarm
% Plots the per-neuron difference (stimA - stimB) as a single swarm column.
% =========================================================================
function randiColors = plotDiffSwarm(ax, tblDiff, pairs, pValues, params, ...
    yMaxVis, bracketPad, textPad)

hold(ax, 'on');
set(ax, 'Clipping', 'off');

randiColors = randperm(height(tblDiff));

% Determine dot color data
if params.colorByZScore && ismember('zScore', tblDiff.Properties.VariableNames)
    colorData = tblDiff.zScore(randiColors);
else
    colorData = tblDiff.animal(randiColors);
end

% Draw swarm
if params.filled
    s = swarmchart(ax, tblDiff.stimulus(randiColors), tblDiff.value(randiColors), ...
        params.dotSize, colorData, 'filled', ...
        'MarkerFaceAlpha', params.Alpha);
else
    s = swarmchart(ax, tblDiff.stimulus(randiColors), tblDiff.value(randiColors), ...
        params.dotSize, colorData, ...
        'MarkerEdgeAlpha', params.Alpha, 'LineWidth', 1, 'SizeData', 30);
end
s.XJitter = params.Xjitter;

% Set colormap
if params.colorByZScore && ismember('zScore', tblDiff.Properties.VariableNames)
    colormap(ax, buildRdBuColormap(256));
    maxZ = max(abs(tblDiff.zScore), [], 'omitnan');
    if maxZ == 0, maxZ = 1; end
    clim(ax, [-maxZ maxZ]);
    cb = colorbar(ax);
    cb.Label.String = 'Z-score';
else
    colormap(ax, lines(numel(categories(tblDiff.animal))));
end

% Zero reference line
yline(ax, 0, 'LineWidth', 1, 'Alpha', 0.7);

ylabel(ax, params.yLegend);
ax.Box   = 'off';
ax.Layer = 'top';

% Bootstrap mean + SEM
if params.plotMeanSem
    stimuli = categories(tblDiff.stimulus);
    plotMeanSemBars(ax, tblDiff, stimuli, params, 0.6);
end

% Significance annotation for diff mode
ylims = ylim(ax);
if ~isempty(pValues) && numel(pValues) >= 1

    fprintf('=== DIFF MODE SIGNIFICANCE ===\n');
    fprintf('p-value: %.4e\n', pValues(1));

    vals = tblDiff.value;
    if isempty(vals)
        fprintf('No values to annotate.\n');
        ylim(ax, [ylims(1) yMaxVis]);
        return
    end

    % Guard against empty vals producing empty maxVisible
    maxVisible = max(min(vals(:), yMaxVis(1)));
    if isempty(maxVisible), maxVisible = yMaxVis; end
    yText = maxVisible + bracketPad;

    if pValues(1) < 1e-3
        txt = '***';
        if pValues(1) == 0, txt = '****'; end

        text(ax, 1, yText, txt, ...
            'HorizontalAlignment', 'center', 'FontSize', 7, 'Clipping', 'off');

        stimA     = pairs{1,1};
        stimB     = pairs{1,2};
        compText  = sprintf('%s > %s', stimA, stimB);
        yCompText = yText + textPad * 10;

        text(ax, 1, yCompText, compText, ...
            'HorizontalAlignment', 'center', 'FontSize', 10, 'Clipping', 'off');

        requiredHeight = yCompText + textPad * 10;
        if requiredHeight > yMaxVis
            ylim(ax, [ylims(1) requiredHeight]);
        else
            ylim(ax, [ylims(1) yMaxVis]);
        end
    else
        ylim(ax, [ylims(1) yMaxVis]);
    end
else
    ylim(ax, [ylims(1) yMaxVis]);
end

fprintf('=== END DIFF MODE SIGNIFICANCE ===\n');

end % plotDiffSwarm


% =========================================================================
% LOCAL FUNCTION: buildDiffTable
% Computes per-neuron difference (stimA - stimB) within each insertion.
% If colorByZScore, also computes mean zScore across the pair for each neuron.
% =========================================================================
function tblDiff = buildDiffTable(tbl, pairs, params)

assert(~isempty(pairs) && size(pairs,1) >= 1, ...
    'diff mode requires at least one stimulus pair.');

stimA = pairs{1,1};
stimB = pairs{1,2};

ins      = categories(tbl.insertion);
diffVals = [];
animals  = [];
insers   = [];
zScores  = [];  % mean zScore across pair, used only if colorByZScore

for i = 1:numel(ins)
    idxA = tbl.insertion == ins{i} & tbl.stimulus == stimA;
    idxB = tbl.insertion == ins{i} & tbl.stimulus == stimB;

    if ~any(idxA) || ~any(idxB), continue; end

    if params.fraction
        vA = tbl.(valueField{1})(idxA) ./ tbl.(valueField{2})(idxA);
        vB = tbl.(valueField{1})(idxB) ./ tbl.(valueField{2})(idxB);
    else
        vA = tbl.value(idxA);
        vB = tbl.value(idxB);
    end

    diffVals = [diffVals; vA - vB];                                       %#ok<AGROW>
    animals  = [animals;  repmat(tbl.animal(find(idxA,1)), length(vA), 1)]; %#ok<AGROW>
    insers   = [insers;   repmat(i, length(vA), 1)];                      %#ok<AGROW>

    % For Z-score coloring: average zScore across both stim types per neuron
    if params.colorByZScore && ismember('zScore', tbl.Properties.VariableNames)
        zA = tbl.zScore(idxA);
        zB = tbl.zScore(idxB);
        zScores = [zScores; (zA + zB) / 2];                               %#ok<AGROW>
    end
end

valid = ~isnan(diffVals);
stimName = sprintf('%s-%s', stimA, stimB);

tblDiff          = table();
tblDiff.insertion = categorical(insers(valid));
tblDiff.stimulus  = categorical(repmat({stimName}, sum(valid), 1));
tblDiff.animal    = animals(valid);
tblDiff.value     = diffVals(valid);

if params.colorByZScore && ~isempty(zScores)
    tblDiff.zScore = zScores(valid);
end

end % buildDiffTable


% =========================================================================
% LOCAL FUNCTION: plotMeanSemBars
% Draws bootstrapped mean ± SEM bars for each stimulus group.
% =========================================================================
function plotMeanSemBars(ax, tblPlot, stimuli, params, semAlpha)

for i = 1:numel(stimuli)
    idx  = tblPlot.stimulus == stimuli{i};
    if ~any(idx), continue; end

    vals = tblPlot.value(idx);
    if numel(vals) < 3
        fprintf('Number of values to bootstrap is less than 3\n');
        continue
    end

    if height(tblPlot) < 500
        bootMean = bootstrp(params.nBoot, @mean, vals);
        mu  = mean(bootMean);
        sem = std(bootMean);
    else
        mu  = mean(vals, 'omitnan');
        sem = std(vals,  'omitnan') / sqrt(sum(~isnan(vals)));
    end

    % SEM error bar
    line(ax, [i i], mu + [-1 1]*sem, ...
        'Color', [0 0 0 semAlpha], 'LineWidth', 2);

    capW = 0.1;
    line(ax, [i-capW i+capW], [mu+sem mu+sem], ...
        'Color', [0 0 0 semAlpha], 'LineWidth', 2);
    line(ax, [i-capW i+capW], [mu-sem mu-sem], ...
        'Color', [0 0 0 semAlpha], 'LineWidth', 2);

    % Mean line
    dx = 0.15;
    plot(ax, [i-dx i+dx], [mu mu], 'k-', 'LineWidth', 1.2);
end

end % plotMeanSemBars


% =========================================================================
% LOCAL FUNCTION: plotBrackets
% Draws pairwise significance brackets above the swarm.
% =========================================================================
function plotBrackets(ax, tblPlot, stimuli, pairs, pValues, ...
    yMaxVis, bracketPad, stackPad, textPad)

fprintf('=== DEBUGGING BRACKETS ===\n');
fprintf('Number of pairs: %d\n', size(pairs,1));
fprintf('Number of pValues: %d\n', numel(pValues));
fprintf('Stimuli in plot: %s\n', strjoin(cellstr(stimuli), ', '));

usedHeights = zeros(size(pairs,1), 1);

for k = 1:size(pairs,1)

    fprintf('\n--- Pair %d: %s vs %s ---\n', k, pairs{k,1}, pairs{k,2});

    x1 = find(strcmp(stimuli, pairs{k,1}));
    x2 = find(strcmp(stimuli, pairs{k,2}));

    fprintf('x1 index: %d, x2 index: %d\n', x1, x2);

    if isempty(x1) || isempty(x2)
        fprintf('SKIPPING: One or both stimuli not found!\n');
        continue
    end

    vals1 = tblPlot.value(tblPlot.stimulus == pairs{k,1});
    vals2 = tblPlot.value(tblPlot.stimulus == pairs{k,2});

    fprintf('vals1 count: %d, vals2 count: %d\n', numel(vals1), numel(vals2));

    maxVisible = max(min([vals1; vals2], yMaxVis));
    yBase = maxVisible + bracketPad;

    y = yBase;
    while any(abs(usedHeights(1:k-1) - y) < stackPad)
        y = y + stackPad;
    end
    usedHeights(k) = y;

    fprintf('Bracket y position: %.3f\n', y);
    fprintf('p-value: %.4e\n', pValues(k));

    % Bracket lines
    line(ax, [x1 x2], [y y], 'Color', 'k', 'LineWidth', 1.2, 'Clipping', 'off');
    line(ax, [x1 x1], [y - yMaxVis*0.01, y], 'Color', 'k', 'LineWidth', 1.2, 'Clipping', 'off');
    line(ax, [x2 x2], [y - yMaxVis*0.01, y], 'Color', 'k', 'LineWidth', 1.2, 'Clipping', 'off');

    % Significance stars
    if pValues(k) < 1e-3
        txt = '***';
        if pValues(k) == 0, txt = '****'; end
        fprintf('Drawing text: %s\n', txt);
        text(ax, mean([x1 x2]), y + textPad, txt, ...
            'HorizontalAlignment', 'center', 'FontSize', 7, 'Clipping', 'off');
    else
        fprintf('p-value not significant enough (>= 1e-3)\n');
    end
end

end % plotBrackets


% =========================================================================
% LOCAL FUNCTION: renameStimulusLabels
% Replaces legacy stimulus abbreviations in tbl.stimulus.
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
% Applies same label substitutions to the pairs cell array.
% =========================================================================
function pairs = renamePairLabels(pairs)
if isempty(pairs), return; end
for i = 1:numel(pairs)
    if strcmp(pairs{i}, 'RG'),   pairs{i} = 'SB'; end
    if strcmp(pairs{i}, 'SDGm'), pairs{i} = 'MG'; end
    if strcmp(pairs{i}, 'SDGs'), pairs{i} = 'SG'; end
end
end


% =========================================================================
% LOCAL FUNCTION: buildRdBuColormap
% Returns an n-row diverging Red-Blue colormap centred at 0.
% Blue = negative (low Z), White = zero, Red = positive (high Z).
% =========================================================================
function cmap = buildRdBuColormap(n)
half = floor(n/2);

% Blue -> White (low to mid)
blueToWhite = [linspace(0.02, 1, half)', ...
               linspace(0.44, 1, half)', ...
               linspace(0.69, 1, half)'];

% White -> Red (mid to high)
whiteToRed  = [linspace(1, 0.70, half)', ...
               linspace(1, 0.09, half)', ...
               linspace(1, 0.09, half)'];

cmap = [blueToWhite; whiteToRed];
end