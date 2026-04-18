function plotPSTH_MultiExp(exList, params)

arguments
    exList double
    params.stimTypes  (1,:) string  = ["rectGrid", "linearlyMovingBall"]
    params.binWidth   double        = 10
    params.smooth     double        = 0      % smoothing window in ms (0 = no smoothing)
    params.statType   string        = "maxPermuteTest"
    params.speed      string        = "max"
    params.alpha      double        = 0.05
    params.shadeSTD   logical       = true
    params.postStim   double        = 500
    params.preBase    double        = 200
    params.overwrite  logical       = false
    params.TakeTopPercentTrials double = 0.3
    params.zScore     logical       = false
    params.PaperFig   logical       = false
    params.byDepth    logical       = false  % plot 3 depth bins per stim type
end

% -------------------------------------------------------------------------
% Load depth info from saved file (only if byDepth is requested)
% -------------------------------------------------------------------------
if params.byDepth
    depthFile = 'W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\NeuronDepths.mat';
    if ~exist(depthFile, 'file')
        error('NeuronDepths.mat not found. Run getNeuronDepths() first.');
    end
    D = load(depthFile);
    depthTable    = D.depthTable;
    depthBinEdges = D.depthBinEdges;
    nDepthBins    = 3;
    fprintf('Depth bins loaded:\n');
    fprintf('  Bin 1 (shallow): %.0f - %.0f um\n', depthBinEdges(1), depthBinEdges(2));
    fprintf('  Bin 2 (middle) : %.0f - %.0f um\n', depthBinEdges(2), depthBinEdges(3));
    fprintf('  Bin 3 (deep)   : %.0f - %.0f um\n', depthBinEdges(3), depthBinEdges(4));
else
    nDepthBins = 1;
end

% -------------------------------------------------------------------------
% Build save path
% -------------------------------------------------------------------------
NP_first = loadNPclassFromTable(exList(1));
vs_first  = linearlyMovingBallAnalysis(NP_first);

p = extractBefore(vs_first.getAnalysisFileName, 'lizards');
p = [p 'lizards'];
if ~exist([p '\Combined_lizard_analysis'], 'dir')
    cd(p)
    mkdir Combined_lizard_analysis
end
saveDir = [p '\Combined_lizard_analysis'];

stimLabel   = strjoin(params.stimTypes, '-');
depthSuffix = '';
if params.byDepth; depthSuffix = '_byDepth'; end
nameOfFile  = sprintf('\\Ex_%d-%d_Combined_PSTHs_%s%s.mat', ...
    exList(1), exList(end), stimLabel, depthSuffix);

% -------------------------------------------------------------------------
% Decide whether to recompute or load
% -------------------------------------------------------------------------
if exist([saveDir nameOfFile], 'file') == 2 && ~params.overwrite
    S = load([saveDir nameOfFile]);
    if isequal(S.expList, exList)
        fprintf('Loading saved PSTHs from:\n  %s\n', [saveDir nameOfFile]);
        forloop = false;
    else
        fprintf('Experiment list mismatch — recomputing.\n');
        forloop = true;
    end
else
    forloop = true;
end

% =========================================================================
% EXPERIMENT LOOP
% =========================================================================
if forloop

    nStim = numel(params.stimTypes);
    nExp  = numel(exList);

    % psthAll{s,b} — s = stim type, b = depth bin (1 if byDepth is off)
    psthAll = cell(nStim, nDepthBins);

    lockedPreBase = [];
    lockedNBins   = [];
    lockedEdges   = [];

    for ei = 1:nExp

        ex = exList(ei);
        fprintf('\n=== Experiment %d ===\n', ex);

        try
            NP = loadNPclassFromTable(ex);
        catch ME
            warning('Could not load experiment %d: %s', ex, ME.message);
            for s = 1:nStim
                for b = 1:nDepthBins
                    if ~isempty(psthAll{s,b})
                        psthAll{s,b} = [psthAll{s,b}; NaN(1, lockedNBins)];
                    end
                end
            end
            continue
        end

        for s = 1:nStim

            stimType = params.stimTypes(s);

            try
                switch stimType
                    case "rectGrid"
                        obj = rectGridAnalysis(NP);
                    case "linearlyMovingBall"
                        obj = linearlyMovingBallAnalysis(NP);
                    case 'StaticGrating'
                        obj = StaticDriftingGratingAnalysis(NP);
                    case 'MovingGrating'
                        obj = StaticDriftingGratingAnalysis(NP);
                    otherwise
                        error('Unknown stimType: %s', stimType);
                end
            catch ME
                warning('Could not build %s for exp %d: %s', stimType, ex, ME.message);
                for b = 1:nDepthBins
                    if ~isempty(psthAll{s,b})
                        psthAll{s,b} = [psthAll{s,b}; NaN(1, lockedNBins)];
                    end
                end
                continue
            end

            NeuronResp = obj.ResponseWindow;

            if params.statType == "BootstrapPerNeuron"
                Stats = obj.BootstrapPerNeuron;
            elseif  params.statType == "maxPermuteTest"
                Stats = obj.StatisticsPerNeuron;
            else
                 Stats = obj.ShufflingAnalysis;
            end

            if params.speed ~= "max" && isequal(obj.stimName,'linearlyMovingBall')
                fieldName = 'Speed2'; startStim = 0;
            elseif isequal(obj.stimName,'linearlyMovingBall')
                fieldName = 'Speed1'; startStim = 0;
            elseif isequal(params.stimTypes,'StaticGrating')
                fieldName = 'Static'; startStim = 0;
            elseif isequal(params.stimTypes,'MovingGrating')
                startStim = obj.VST.static_time*1000; fieldName = 'Moving';
            else
                startStim = 0;
            end

            p_sort = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
            label  = string(p_sort.label');
            goodU  = p_sort.ic(:, label == 'good');

            try
                pvals = Stats.(fieldName).pvalsResponse;
            catch
                pvals = Stats.pvalsResponse;
            end

            try
                C = NeuronResp.(fieldName).C;
            catch
                C = NeuronResp.C;
            end
            directimesSorted = C(:, 1)' + startStim;

            preBase     = params.preBase;
            windowTotal = preBase + params.postStim;

            if isempty(lockedPreBase)
                lockedPreBase = preBase;
                lockedEdges   = 0 : params.binWidth : windowTotal;
                lockedNBins   = numel(lockedEdges) - 1;
                tAxis         = lockedEdges(1:end-1);
                fprintf('  Locked window: preBase=%d ms, postStim=%d ms, nBins=%d\n', ...
                    lockedPreBase, params.postStim, lockedNBins);
            end

            eNeurons = find(pvals < params.alpha);

            if isempty(eNeurons)
                fprintf('  [%s] No responsive neurons in exp %d.\n', stimType, ex);
                for b = 1:nDepthBins
                    if ~isempty(psthAll{s,b})
                        psthAll{s,b} = [psthAll{s,b}; NaN(1, lockedNBins)];
                    end
                end
                continue
            end

            fprintf('  [%s] %d responsive neuron(s) in exp %d.\n', stimType, ex, numel(eNeurons));

            % ----------------------------------------------------------
            % Build PSTH per neuron
            % ----------------------------------------------------------
            psthRateNeurons = zeros(numel(eNeurons), lockedNBins);
            neuronBinIdx    = zeros(numel(eNeurons), 1);

            for ni = 1:numel(eNeurons)
                u = eNeurons(ni);

                % Assign depth bin
                if params.byDepth
                    depthRow = depthTable.Experiment == ex & depthTable.Unit == u;
                    if ~any(depthRow)
                        neuronBinIdx(ni) = 0;  % unknown depth — skip
                        continue
                    end
                    unitDepth = depthTable.Depth_um(depthRow);
                    if unitDepth <= depthBinEdges(2)
                        neuronBinIdx(ni) = 1;
                    elseif unitDepth <= depthBinEdges(3)
                        neuronBinIdx(ni) = 2;
                    else
                        neuronBinIdx(ni) = 3;
                    end
                else
                    neuronBinIdx(ni) = 1;  % all neurons in single bin
                end

                MRhist = BuildBurstMatrix( ...
                    goodU(:, u), ...
                    round(p_sort.t), ...
                    round(directimesSorted - lockedPreBase), ...
                    round(windowTotal));
                MRhist = squeeze(MRhist);

                if ~isempty(params.TakeTopPercentTrials)
                    MeanTrial  = mean(MRhist, 2);
                    [~, ind]   = sort(MeanTrial, 'descend');
                    takeTrials = ind(1:round(numel(MeanTrial)*params.TakeTopPercentTrials));
                    MRhist     = MRhist(takeTrials, :);
                end

                nTrials    = size(MRhist, 1);
                spikeTimes = repmat((1:size(MRhist,2)), nTrials, 1);
                spikeTimes = spikeTimes(logical(MRhist));
                counts     = histcounts(spikeTimes, lockedEdges);
                psthRateNeurons(ni, :) = (counts / (params.binWidth * nTrials)) * 1000;
            end

            % ----------------------------------------------------------
            % Average per depth bin and append
            % ----------------------------------------------------------
            for b = 1:nDepthBins
                binNeurons = neuronBinIdx == b;
                if ~any(binNeurons)
                    fprintf('  [%s] No neurons in depth bin %d for exp %d.\n', stimType, b, ex);
                    if ~isempty(psthAll{s,b})
                        psthAll{s,b} = [psthAll{s,b}; NaN(1, lockedNBins)];
                    end
                    continue
                end

                psthExp = mean(psthRateNeurons(binNeurons, :), 1, 'omitnan');

                if params.zScore
                    baselineBins = tAxis < lockedPreBase;
                    baselineMean = mean(psthExp(baselineBins));
                    baselineStd  = std(psthExp(baselineBins));
                    if baselineStd > 0
                        psthExp = (psthExp - baselineMean) / baselineStd;
                    else
                        warning('  [%s] Bin %d: baseline std is zero for exp %d.', stimType, b, ex);
                        if ~isempty(psthAll{s,b})
                            psthAll{s,b} = [psthAll{s,b}; NaN(1, lockedNBins)];
                        end
                        continue
                    end
                end

                psthAll{s,b} = [psthAll{s,b}; psthExp(:)'];
                fprintf('  [%s] Bin %d: %d neuron(s) in exp %d.\n', stimType, b, sum(binNeurons), ex);
            end

        end % stim loop
    end % experiment loop

    % ------------------------------------------------------------------
    % Save
    % ------------------------------------------------------------------
    S.expList       = exList;
    S.lockedEdges   = lockedEdges;
    S.lockedPreBase = lockedPreBase;
    S.params        = params;

    for s = 1:numel(params.stimTypes)
        stimField = matlab.lang.makeValidName(params.stimTypes(s));
        for b = 1:nDepthBins
            S.(sprintf('%s_bin%d', stimField, b)) = psthAll{s,b};
        end
    end

    save([saveDir nameOfFile], '-struct', 'S');
    fprintf('\nSaved PSTHs to:\n  %s\n', [saveDir nameOfFile]);

else
    % Load psthAll from disk
    lockedEdges   = S.lockedEdges;
    lockedPreBase = S.lockedPreBase;

    psthAll = cell(numel(params.stimTypes), nDepthBins);
    for s = 1:numel(params.stimTypes)
        stimField = matlab.lang.makeValidName(params.stimTypes(s));
        for b = 1:nDepthBins
            fieldKey = sprintf('%s_bin%d', stimField, b);
            if isfield(S, fieldKey)
                psthAll{s,b} = S.(fieldKey);
            else
                warning('Field "%s" not found in saved file.', fieldKey);
                psthAll{s,b} = [];
            end
        end
    end
end

% =========================================================================
% PLOT
% =========================================================================

tAxis     = lockedEdges(1:end-1);
tAxisPlot = tAxis - lockedPreBase;

baseColors  = lines(numel(params.stimTypes));
depthShades = [0.05, 0.45, 0.78];  % light → dark for shallow → deep
binLabels   = {'shallow', 'middle', 'deep'};

stimLegendMap = containers.Map(...
    {'linearlyMovingBall', 'rectGrid', 'MovingGrating', 'StaticGrating'}, ...
    {'MB',                 'SB',       'MG',            'SG'});

% ------------------------------------------------------------------
% First pass: global ylim
% ------------------------------------------------------------------
yMax = 0;
yMin = inf;

meanAll = cell(numel(params.stimTypes), nDepthBins);
semAll  = cell(numel(params.stimTypes), nDepthBins);

for s = 1:numel(params.stimTypes)
    for b = 1:nDepthBins
        data = psthAll{s,b};
        if isempty(data); continue; end
        validRows = ~all(isnan(data), 2);
        data      = data(validRows, :);
        if isempty(data); continue; end
        meanAll{s,b} = mean(data, 1, 'omitnan');
        semAll{s,b}  = std(data, 0, 1, 'omitnan') / sqrt(sum(~isnan(data(:,1))));
        yMax = max(yMax, max(meanAll{s,b} + semAll{s,b}));
        yMin = min(yMin, min(meanAll{s,b} - semAll{s,b}));
    end
end

yPad = (yMax - yMin) * 0.1;
if params.zScore
    yLims = [yMin - yPad, yMax + yPad];
else
    yLims = [max(0, yMin - yPad), yMax + yPad];
end

% ------------------------------------------------------------------
% Plot
% ------------------------------------------------------------------
fig = figure;
set(fig, 'Units', 'centimeters', 'Position', [5 5 9 10]);
ax = axes(fig);
hold(ax, 'on');

legendHandles = [];
legendLabels  = {};

for s = 1:numel(params.stimTypes)

    stimKey = char(params.stimTypes(s));
    if isKey(stimLegendMap, stimKey)
        shortName = stimLegendMap(stimKey);
    else
        shortName = stimKey;
    end

    for b = 1:nDepthBins

        data = psthAll{s,b};
        if isempty(data); continue; end
        validRows = ~all(isnan(data), 2);
        data      = data(validRows, :);
        if isempty(data); continue; end

        meanPSTH = meanAll{s,b};
        semPSTH  = semAll{s,b};

        % Smooth if requested
        if params.smooth > 0
            smoothBins = round(params.smooth / params.binWidth);  % convert ms to bins
            meanPSTH   = smoothdata(meanPSTH, 'gaussian', smoothBins);
            semPSTH    = smoothdata(semPSTH,  'gaussian', smoothBins);
        end

        % Color and label depend on mode
        if params.byDepth
            lineColor   = baseColors(s,:) * (1 - depthShades(b));
            legendLabel = sprintf('%s %s (%.0f-%.0f um)', ...
                shortName, binLabels{b}, depthBinEdges(b), depthBinEdges(b+1));
        else
            lineColor   = baseColors(s,:);
            legendLabel = shortName;
        end

        % SEM shading
        if params.shadeSTD && size(data,1) > 1
            upper = meanPSTH + semPSTH;
            lower = meanPSTH - semPSTH;
            xFill = [tAxisPlot(:)', fliplr(tAxisPlot(:)')];
            yFill = [upper(:)',     fliplr(lower(:)')    ];
            fill(ax, xFill, yFill, lineColor, 'FaceAlpha', 0.08, 'EdgeColor', 'none');
        end

        % Mean line
        h = plot(ax, tAxisPlot(:)', meanPSTH(:)', ...
            'Color', lineColor, 'LineWidth', 1.5);

        legendHandles(end+1) = h; %#ok<AGROW>
        legendLabels{end+1}  = legendLabel; %#ok<AGROW>

        fprintf('  [%s] n=%d experiments in plot.\n', legendLabel, sum(validRows));
    end
end

xline(ax, 0,               'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(ax, params.postStim, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

if params.zScore; yLabel = 'Z-score'; else; yLabel = '[spk/s]'; end

xlabel(ax, 'Time re. stim onset [ms]', 'FontName', 'helvetica', 'FontSize', 8);
ylabel(ax, yLabel,                      'FontName', 'helvetica', 'FontSize', 8);
xlim(ax, [tAxisPlot(1) tAxisPlot(end)]);
ylim(ax, yLims);

legend(legendHandles, legendLabels, 'Location', 'northwest', ...
    'FontName', 'helvetica', 'FontSize', 7);

ax.FontName       = 'helvetica';
ax.FontSize       = 8;
ax.YAxis.FontSize = 8;
ax.XAxis.FontSize = 8;
hold(ax, 'off');

sgtitle(sprintf('N = %d', numel(exList)), 'FontName', 'helvetica', 'FontSize', 11);
set(fig, 'Units', 'centimeters', 'Position', [20 20 7 4]);

if params.PaperFig
    vs_first.printFig(fig, sprintf('PSTH-depth-%s-%s', ...
        params.stimTypes(1), params.stimTypes(2)), PaperFig = params.PaperFig)
end

end