function plotPairwiseComparisons(stimuli, pairs, pValues, yMax)

nPairs = size(pairs,1);
offset = 5;   % vertical spacing above yMax for the first line
step   = 5;   % spacing between multiple comparison lines

for k = 1:nPairs
    % find x positions
    x1 = find(strcmp(stimuli, pairs{k,1}));
    x2 = find(strcmp(stimuli, pairs{k,2}));

    % y position for this comparison
    y = yMax + offset + (k-1)*step;

    % horizontal line connecting the two stimuli
    plot([x1 x2], [y y], 'k-', 'LineWidth',1.5)

    % short vertical lines down to x1 and x2
    plot([x1 x1], [y-1 y], 'k-', 'LineWidth',1.5)
    plot([x2 x2], [y-1 y], 'k-', 'LineWidth',1.5)

    % mark significance / probability
    if pValues(k) < 0.05
        text(mean([x1 x2]), y+1, '*', ...
            'HorizontalAlignment','center', 'FontSize',14)
    else
        text(mean([x1 x2]), y+1, sprintf('p=%.2f', pValues(k)), ...
            'HorizontalAlignment','center', 'FontSize',10)
    end
end
end
