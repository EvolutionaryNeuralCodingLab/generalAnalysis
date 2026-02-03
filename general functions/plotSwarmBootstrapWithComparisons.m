function plotSwarmBootstrapWithComparisons(tbl, pairs, pValues, valueField, params)

arguments
    tbl table
    pairs cell = {}
    pValues double = []
    valueField cell = {}
    params.nBoot (1,1) double = 10000
    params.fraction logical = false
    params.yLegend char = 'value'
    params.diff logical = false   % compute difference between first pair
    params.Xjitter = 'density'
    params.dotSize = 7
    params.yMaxVis = 1
    params.filled = true
    params.Alpha = 0.4
    params.plotMeanSem = true
end

%% ----------------- PARAMETERS -----------------
yMaxVis     = params.yMaxVis;
bracketPad  = yMaxVis * 0.05;
stackPad    = yMaxVis * 0.05;
textPad     = yMaxVis * 0.01;
semAlpha    = 0.6;

%% ----------------- DIFF MODE -----------------
if params.diff

    assert(~isempty(pairs) && size(pairs,1) >= 1, ...
        'params.diff=true requires at least one stimulus pair.');

    stimA = pairs{1,1};
    stimB = pairs{1,2};

    if params.fraction
        assert(numel(valueField) == 2, ...
            'Fraction mode requires two valueField entries.');
    end


    ins = categories(tbl.insertion);

    diffVals = [];
    animals = [];
    insers = [];

    for i = 1:numel(ins)
        idxA = tbl.insertion == ins{i} & tbl.stimulus == stimA;
        idxB = tbl.insertion == ins{i} & tbl.stimulus == stimB;

        if any(idxA) && any(idxB)

            if params.fraction
                vA = tbl.(valueField{1})(idxA) ./ tbl.(valueField{2})(idxA);
                vB = tbl.(valueField{1})(idxB) ./ tbl.(valueField{2})(idxB);
            else
                vA = tbl.(valueField{1})(idxA);
                vB = tbl.(valueField{1})(idxB);
            end

            diffVals = [diffVals; vA - vB];
            animals  = [animals; repmat(tbl.animal(find(idxA,1)),length(vA),1)]; 
            insers = [insers; repmat(i,length(vA),1)];
        end
    end

    valid = ~isnan(diffVals);

    stimName = sprintf('%s-%s', stimA, stimB);

    tblPlot = table();
    tblPlot.insertion = categorical(insers(valid));
    tblPlot.stimulus  = categorical(repmat({stimName}, sum(valid), 1));
    tblPlot.animal    = animals(valid);
    tblPlot.value     = diffVals(valid);

    plotLinesBetween = false;

else
    %% ----------------- RAW MODE -----------------
    if params.fraction
        assert(numel(valueField) == 2, ...
            'Fraction mode requires two valueField entries.');
        tbl.value = tbl.(valueField{1}) ./ tbl.(valueField{2});
    else
        tbl.value = tbl.(valueField{1});
    end

    tblPlot = tbl;
    plotLinesBetween = true;
end

%% ----------------- CHANGE STIM NAMES -----------------
stimuli = unique(tblPlot.stimulus);

% tbl.stimulus = removecats(tbl.stimulus);
% tbl.animal = removecats(tbl.animal);
% tbl.insertion = removecats(tbl.insertion);
%Replace 'RG' with 'SB'

% Convert to string
s = string(tblPlot.stimulus);

% Replace substring wherever it appears
s = replace(s, "RG", "SB");
s = replace(s, "SDGs", "SG");
s = replace(s, "SDGm", "MG");


% Convert back to categorical
tblPlot.stimulus = categorical(s);


%% ----------------- RENAME PAIRS TO MATCH -----------------
if ~isempty(pairs)
    for i = 1:numel(pairs)
        if strcmp(pairs{i}, 'RG')
            pairs{i} = 'SB';
        elseif strcmp(pairs{i}, 'SDGm')
            pairs{i} = 'MG';
        elseif strcmp(pairs{i}, 'SDGs')
            pairs{i} = 'SG';
        end
    end
end


%% ----------------- CLEAN CATEGORIES -----------------
tblPlot.stimulus   = removecats(tblPlot.stimulus);
tblPlot.animal     = removecats(tblPlot.animal);
tblPlot.insertion  = removecats(tblPlot.insertion);

stimuli     = categories(tblPlot.stimulus);
insertions  = categories(tblPlot.insertion);

%% ----------------- FIGURE -----------------
figure;
ax = axes;
hold(ax, 'on');
set(ax, 'Clipping', 'off');  % <-- ADD THIS LINE

% 1) Swarm first
if params.filled
    s=swarmchart(ax, tblPlot.stimulus, tblPlot.value, ...
        params.dotSize, tblPlot.animal, 'filled', ...
        'MarkerFaceAlpha', params.Alpha);
else
    s=swarmchart(ax, tblPlot.stimulus, tblPlot.value, ...
        params.dotSize, tblPlot.animal, ...
        'MarkerEdgeAlpha',params.Alpha);
end

s.XJitter = params.Xjitter;
%s.XJitterWidth = 0.1;

if plotLinesBetween
    % 2) Get numeric x positions of categories
    cats = categories(tblPlot.stimulus);
    xMap = containers.Map(cats, 1:numel(cats));

    xNum = arrayfun(@(i) xMap(char(tblPlot.stimulus(i))), ...
        1:height(tblPlot));


    % 3) Plot lines AFTER swarm
    for i = 1:numel(unique(tblPlot.NeurID))
        idx = double(tblPlot.NeurID) == i;
        if sum(idx) < 2
            continue
        end

        line(ax, ...
            xNum(idx), tblPlot.value(idx), ...
            'Color', [0 0 0 0.1], ...
            'LineWidth', 0.1),'lin';
    end
else
    yline(0,LineWidth=1,Alpha=0.7)
end
colormap(lines(numel(categories(tblPlot.animal))))
ylabel(params.yLegend)

ax = gca;
ax.Box   = 'off';
ax.Layer = 'top';

%% ----------------- BOOTSTRAP MEAN + SEM -----------------

if params.plotMeanSem

    for i = 1:numel(stimuli)

        idx = tblPlot.stimulus  == stimuli{i};

        if any(idx) && params.fraction
            vals = tblPlot.value(idx);
            insers = tblPlot.insertion(idx);
            animals = tblPlot.animal(idx);
        elseif any(idx)
            vals = tblPlot.value(idx);
            insers = tblPlot.insertion(idx);
            animals = tblPlot.animal(idx);
        end

        if numel(vals) < 3
            fprintf('Number of values to bootstrap is less than 3\n')
            continue
        end

        if size(tblPlot,1) < 500
            bootMean = bootstrp(params.nBoot, @mean, vals);
            mu  = mean(bootMean);
            sem = std(bootMean);
        else

            mu = mean(vals);
            sem = std(vals,'omitnan') / sqrt(numel(vals));
        end

        % SEM
        line([i i], mu + [-1 1]*sem, ...
            'Color', [0 0 0 semAlpha], 'LineWidth', 2);

        capW = 0.1;
        line([i-capW i+capW], [mu+sem mu+sem], ...
            'Color', [0 0 0 semAlpha], 'LineWidth', 2);
        line([i-capW i+capW], [mu-sem mu-sem], ...
            'Color', [0 0 0 semAlpha], 'LineWidth', 2);


        % mean
        dx = 0.15;
        plot([i-dx i+dx], [mu mu], 'k-', 'LineWidth', 1.2);
    end
end


%% ----------------- PAIRWISE COMPARISONS -----------------
if ~params.diff && ~isempty(pairs) && numel(pValues) == size(pairs,1)

    fprintf('=== DEBUGGING BRACKETS ===\n');
    fprintf('Number of pairs: %d\n', size(pairs,1));
    fprintf('Number of pValues: %d\n', numel(pValues));
    fprintf('Stimuli in plot: %s\n', strjoin(cellstr(stimuli), ', '));
    
    usedHeights = zeros(size(pairs,1),1);

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
        fprintf('Will draw bracket: YES\n');

        % bracket
        line([x1 x2], [y y], 'Color','k', 'LineWidth',1.2, 'Clipping', 'off');
        line([x1 x1], [y-yMaxVis*0.01 y], 'Color','k', 'LineWidth',1.2, 'Clipping', 'off');
        line([x2 x2], [y-yMaxVis*0.01 y], 'Color','k', 'LineWidth',1.2, 'Clipping', 'off');

        % significance
        if pValues(k) < 1e-3
            txt = '***';
            if pValues(k) == 0
                txt = '****';
            end
            fprintf('Drawing text: %s\n', txt);
            text(mean([x1 x2]), y + textPad, txt, ...
                'HorizontalAlignment','center', ...
                'FontSize', 7, 'Clipping', 'off');
        else
            fprintf('p-value not significant enough (>= 1e-3)\n');
        end
    end
    
end

%% ----------------- SIGNIFICANCE FOR DIFF MODE -----------------
%% ----------------- SIGNIFICANCE FOR DIFF MODE -----------------
if params.diff && ~isempty(pValues) && numel(pValues) >= 1
    
    fprintf('=== DIFF MODE SIGNIFICANCE ===\n');
    fprintf('p-value: %.4e\n', pValues(1));
    
    % There's only one "stimulus" (the difference)
    x1 = 1;  % Position of the single difference bar
    
    % Get all values for this difference
    vals = tblPlot.value;
    
    % Find the maximum visible value
    maxVisible = max(min(vals, yMaxVis));
    yText = maxVisible + bracketPad;
    
    fprintf('Text y position: %.3f\n', yText);
    
    % Draw significance stars
    if pValues(1) < 1e-3
        txt = '***';
        if pValues(1) == 0
            txt = '****';
        end
        fprintf('Drawing significance: %s\n', txt);
        
        % Draw the asterisks
        text(x1, yText, txt, ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 7, ...
            'Clipping', 'off');
        
        % Draw the comparison text above the asterisks
        stimA = pairs{1,1};
        stimB = pairs{1,2};
        compText = sprintf('%s > %s', stimA, stimB);
        
        yCompText = yText + textPad*10;
        
        text(x1, yCompText, compText, ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 10, ...
            'Clipping', 'off');
        
        fprintf('Drawing comparison text: %s\n', compText);

        ylims = ylim;
        
        % Adjust ylim if needed to fit both texts
        requiredHeight = yCompText + textPad*10;  % Extra padding above comparison text
        if requiredHeight > yMaxVis
            ylim([ylims(1) requiredHeight]);
            fprintf('Adjusted ylim to [0 %.3f] to fit text\n', requiredHeight);
        else
            ylim([ylims(1) yMaxVis]);
        end
    else
        fprintf('p-value not significant enough (>= 1e-3)\n');
        ylim([0 yMaxVis]);
    end
    
    fprintf('=== END DIFF MODE SIGNIFICANCE ===\n');
else
    ylim([ylims(1) yMaxVis]);
end


end