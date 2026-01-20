function plotSwarmBootstrapWithComparisons(tbl, pairs, pValues,valueField, params)
arguments
    tbl table
    pairs cell = {} %Pair names (abreviations)
    pValues double = [] %P values of comparison pairs
    valueField cell = {} %Valuues of points to plot
    params.nBoot (1,1) double = 10000
    params.fraction = false %If there are 2 value fields, then a fraction between the 1st and the 2nd is computed
    params.yLegend = 'value'
end
%% ----------------- PARAMETERS -----------------
yMaxVis     = 1;   % visualization limit
bracketPad  = yMaxVis*0.05;     % space above visible data
stackPad    = yMaxVis*0.05;     % vertical spacing between brackets (increased)
textPad     = yMaxVis*0.01;   % distance between p-value text and bracket
semAlpha    = 0.6;

%% ----------------- PARAMETERS -----------------
if params.fraction

    values = tbl.(valueField{1})./tbl.(valueField{2});
else
    values = tbl.(valueField{1});
end


%% ----------------- SWARM -----------------
figure
tbl.stimulus = removecats(tbl.stimulus);

tbl.animal = removecats(tbl.animal);

swarmchart(tbl.stimulus, values, 18, tbl.animal, ...
    'filled', 'MarkerFaceAlpha', 0.6);
hold on
colormap(lines(numel(unique(tbl.animal))))
ylabel(params.yLegend)
%xlabel('Stimulus')
ax = gca;
ax.Box   = 'off';    % remove top/right axes
ax.Layer = 'top';
stimuli     = unique(tbl.stimulus);
experiments = unique(tbl.insertion);

% Replace 'RG' with 'SB' in x-axis labels
stimuliLabels = stimuli;
rgIdx = ismember(stimuliLabels, categorical({'RG'}));
if any(rgIdx)
    stimuliLabels(rgIdx) = 'SB';
end

rgIdx = ismember(stimuliLabels, categorical({'SDGm'}));
if any(rgIdx)
    stimuliLabels(rgIdx) = 'MG';
end

rgIdx = ismember(stimuliLabels, categorical({'SDGs'}));
if any(rgIdx)
    stimuliLabels(rgIdx) = 'SG';
end

ax.XTickLabel = stimuliLabels;

%% ----------------- BOOTSTRAP SEM -----------------
for i = 1:numel(stimuli)
    % one value per experiment
    vals = [];
    for e = 1:numel(experiments)
        idx = tbl.insertion == experiments(e) & ...
              tbl.stimulus   == stimuli(i);
        if any(idx)
            vals(end+1,1) = values(idx); %#ok<AGROW>
        end
    end
    % require >= 3 experiments
    if numel(vals) < 3
        continue
    end
    bootMean = bootstrp(params.nBoot, @mean, vals);
    mu  = mean(bootMean);
    sem = std(bootMean);
    % SEM vertical line
    line([i i], mu + [-1 1]*sem, ...
        'Color', [0 0 0 semAlpha], 'LineWidth', 2);
    % SEM caps
    capWidth = 0.1;
    line([i-capWidth i+capWidth], [mu+sem mu+sem], ...
        'Color', [0 0 0 semAlpha], 'LineWidth', 2);
    line([i-capWidth i+capWidth], [mu-sem mu-sem], ...
        'Color', [0 0 0 semAlpha], 'LineWidth', 2);
    % Mean
    dx = 0.15;
    plot([i-dx i+dx], [mu mu], 'k-', 'LineWidth', 1.2);
end
%% ----------------- PAIRWISE COMPARISONS -----------------
% Safety check: ensure pValues matches pairs
if size(pairs,1) ~= numel(pValues)
    warning('Number of pairs (%d) does not match number of pValues (%d). Skipping comparisons.', ...
        size(pairs,1), numel(pValues));
    ylim([0 yMaxVis])
    hold off
    return
end

usedHeights = zeros(size(pairs,1),1);
for k = 1:size(pairs,1)
    x1 = find(ismember(stimuli, pairs{k,1}));
    x2 = find(ismember(stimuli, pairs{k,2}));
    
    % Check if both stimuli exist
    if isempty(x1) || isempty(x2)
        warning('Stimulus pair {%s, %s} not found in data. Skipping.', ...
            pairs{k,1}, pairs{k,2});
        continue
    end
    
    % Find max visible data for ONLY these two stimuli
    idx1 = tbl.stimulus == pairs{k,1};
    idx2 = tbl.stimulus == pairs{k,2};
    
    % Get values below yMaxVis for each stimulus
    vals1 = values(idx1);
    vals2 = values(idx2);
    vals1_visible = vals1(vals1 <= yMaxVis);
    vals2_visible = vals2(vals2 <= yMaxVis);
    
    % Find the maximum visible value (closest to 100)
    if isempty(vals1_visible)
        max1 = yMaxVis;
    else
        max1 = max(vals1_visible);
    end
    
    if isempty(vals2_visible)
        max2 = yMaxVis;
    else
        max2 = max(vals2_visible);
    end
    
    maxVisible = max(max1, max2);
    
    yBase = maxVisible + bracketPad;
    
    % auto-stack brackets
    y = yBase;
    while any(abs(usedHeights(1:k-1) - y) < stackPad)
        y = y + stackPad;
    end
    usedHeights(k) = y;
    
    % bracket
    line([x1 x2], [y y], 'Color','k', 'LineWidth',1.2);
    line([x1 x1], [y-yMaxVis*0.01 y], 'Color','k', 'LineWidth',1.2);
    line([x2 x2], [y-yMaxVis*0.01 y], 'Color','k', 'LineWidth',1.2);
    
    % probability label (only show if p < 1e-3)
    if pValues(k) < 1e-3
        if pValues(k) == 0
            %pTxt = sprintf('p < %.0e', 1/params.nBoot);
            pTxt = '****';
        else
            %pTxt = sprintf('p = %.4f', pValues(k));
            pTxt = '***';
        end
        text(mean([x1 x2]), y + textPad, pTxt, ...
            'HorizontalAlignment','center', ...
            'FontSize', 7);
    end
end
%% ----------------- OVERFLOW MARKERS -----------------
for i = 1:numel(stimuli)
    nAbove = sum(tbl.stimulus == stimuli(i) & values > yMaxVis);
    if nAbove == 0
        continue
    end
    % highest bracket involving this stimulus
    involved = any(strcmp(pairs, stimuli(i)),2);
    if any(involved)
        yBracket = max(usedHeights(involved));
        % Position text above the bracket area
        yText = yBracket + 4;
    else
        yText = yMaxVis + 2;
    end
    
    % Only show text label
    text(i, yText, sprintf('n=%d', nAbove), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize', 10, ...
        'Color','r', ...
        'Clipping','off');
end
ylim([0 yMaxVis])
hold off
end