function plotZScoreReceptiveField(obj, neuronIdx, direction, speedIdx)
% plotZScoreReceptiveField - Visualise grid-based z-score receptive field.
%
% Inputs:
%   neuronIdx : index of neuron to plot
%   direction : 'best' (default) or a specific direction index
%   speedIdx  : 1 or 2 (default 1)

    if nargin < 3, direction = 'best'; end
    if nargin < 4, speedIdx = 1; end

    % Load cached results
    if ~isfile(obj.getAnalysisFileName)
        error('Run StatisticsPerNeuronSpatialGrid first.');
    end
    S = load(obj.getAnalysisFileName);
    fieldName = sprintf('Speed%d', speedIdx);
    if ~isfield(S, fieldName)
        error('Speed%d not present in results.', speedIdx);
    end

    data = S.(fieldName);

    % Pick direction
    if strcmp(direction, 'best')
        dirIdx = data.prefDirection(neuronIdx);
    else
        dirIdx = direction;
    end

    % Extract mean Diff per cell for this neuron at this direction
    % Use pooled across all factor levels — average the ZScorePerGrid across levels
    factorFields = fieldnames(data.ZScorePerGrid);
    fName        = factorFields{1};  % use first factor's array (all have same dirs)
    zArr         = data.ZScorePerGrid.(fName);   % [nCells × nDirs × nLevels × nNeurons]

    zCell = squeeze(mean(zArr(:, dirIdx, :, neuronIdx), 3, 'omitnan'));  % [nCells × 1]

    % Reshape to grid
    gridSize = data.gridSize;
    zGrid    = reshape(zCell, gridSize, gridSize)';   % transpose because (gy-1)*nGrid + gx

    % Plot
    figure;
    imagesc(zGrid);
    axis image;
    colorbar;
    title(sprintf('Neuron %d — Direction %d — Speed %d', neuronIdx, dirIdx, speedIdx));
    xlabel('grid x');
    ylabel('grid y');
    colormap(parula);
end