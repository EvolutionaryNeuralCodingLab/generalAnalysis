function plotPSTH_MultiExpV1(exList, params)

arguments
    exList double
    params.stimTypes  (1,:) string  = ["rectGrid", "linearlyMovingBall"]
    params.bin        double        = 30
    params.binWidth   double        = 10
    params.statType   string        = "BootstrapPerNeuron"
    params.speed      string        = "max"
    params.alpha      double        = 0.05
    params.shadeSTD   logical       = true
    params.postStim   double        = 500    % ms after stim onset to include
    params.preBase    double        = 200    % ms of baseline before stim onset
    params.overwrite  logical       = false  % force recompute even if file exists
    params.TakeTopPercentTrials double = 0.3 %Percentage of highest spiking rate trials to take to calculate PSTHs
    params.zScore     logical       = false  % normalize firing rate to z-score using baseline
     params.PaperFig  logical       = false  %Is this going to be used in the paper? 
end

% -------------------------------------------------------------------------
% Build save path using first experiment to get the analysis folder
% This mirrors the convention used in PlotZScoreComparison
% -------------------------------------------------------------------------

% Load first experiment just to get the folder path
NP_first = loadNPclassFromTable(exList(1));
vs_first  = linearlyMovingBallAnalysis(NP_first);  % used only for path

% Build the save directory path
p = extractBefore(vs_first.getAnalysisFileName, 'lizards');
p = [p 'lizards'];
if ~exist([p '\Combined_lizard_analysis'], 'dir')
    cd(p)
    mkdir Combined_lizard_analysis
end
saveDir = [p '\Combined_lizard_analysis'];

% Build filename — includes stim types so different comparisons don't clash
stimLabel   = strjoin(params.stimTypes, '-');  % e.g. "rectGrid-linearlyMovingBall"
nameOfFile  = sprintf('\\Ex_%d-%d_Combined_PSTHs_%s.mat', ...
    exList(1), exList(end), stimLabel);

% -------------------------------------------------------------------------
% Decide whether to run the experiment loop or load from disk
% forloop = true  → compute PSTHs from scratch
% forloop = false → load saved struct and skip to plotting
% -------------------------------------------------------------------------
if exist([saveDir nameOfFile], 'file') == 2 && ~params.overwrite
    % File exists and overwrite is off — check if expList matches
    S = load([saveDir nameOfFile]);
    if isequal(S.expList, exList)
        fprintf('Loading saved PSTHs from:\n  %s\n', [saveDir nameOfFile]);
        forloop = false;  % skip computation, go straight to plot
    else
        fprintf('Experiment list mismatch — recomputing.\n');
        forloop = true;   % expList changed, recompute
    end
else
    forloop = true;  % file doesn't exist or overwrite requested
end

% =========================================================================
% EXPERIMENT LOOP — only runs if forloop is true
% =========================================================================
if forloop

    nStim = numel(params.stimTypes);
    nExp  = numel(exList);

    % One cell per stim type, grows one row per experiment
    psthAll = cell(1, nStim);
    for s = 1:nStim
        psthAll{s} = [];
    end

    % Locked time window — set from first valid experiment
    lockedPreBase  = [];
    lockedNBins    = [];
    lockedEdges    = [];

    % ------------------------------------------------------------------
    % LOOP OVER EXPERIMENTS
    % ------------------------------------------------------------------
    for ei = 1:nExp

        ex = exList(ei);
        fprintf('\n=== Experiment %d ===\n', ex);

        % Load NP data for this experiment
        try
            NP = loadNPclassFromTable(ex);
        catch ME
            warning('Could not load experiment %d: %s', ex, ME.message);
            % Add NaN placeholder row if window is already locked
            for s = 1:nStim
                if ~isempty(psthAll{s})
                    psthAll{s} = [psthAll{s}; NaN(1, lockedNBins)];
                end
            end
            continue
        end

        % --------------------------------------------------------------
        % LOOP OVER STIMULUS TYPES
        % --------------------------------------------------------------
        for s = 1:nStim

            stimType = params.stimTypes(s);

            % Build analysis object for this stim type
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
                if ~isempty(psthAll{s})
                    psthAll{s} = [psthAll{s}; NaN(1, lockedNBins)];
                end
                continue
            end

            % ----------------------------------------------------------
            % Extract data structures
            % ----------------------------------------------------------

            % ResponseWindow holds trial timing and spike data
            NeuronResp = obj.ResponseWindow;

            % Stats struct for p-values
            if params.statType == "BootstrapPerNeuron"
                Stats = obj.BootstrapPerNeuron;
            else
                Stats = obj.ShufflingAnalysis;
            end

            % Resolve speed field name
            if params.speed ~= "max" && isequal(obj.stimName,'linearlyMovingBall')
                fieldName = 'Speed2';
                startStim = 0;
            elseif isequal(obj.stimName,'linearlyMovingBall')
                fieldName = 'Speed1';
                startStim = 0;
            elseif isequal(params.stimTypes,'StaticGrating')
                fieldName = 'Static';
                startStim = 0;

            elseif isequal(params.stimTypes,'MovingGrating')
                startStim = obj.VST.static_time*1000;
                fieldName = 'Moving';
            else
                startStim = 0;
            end

            % Spike trains of somatic (good) units
            p_sort = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
            label  = string(p_sort.label');
            goodU  = p_sort.ic(:, label == 'good');

            % P-values for each unit
            try
                pvals   = Stats.(fieldName).pvalsResponse;
            catch
                pvals   = Stats.pvalsResponse;
            end

            % Trial onset times in ms
            try
            C  = NeuronResp.(fieldName).C;
            catch
                C = NeuronResp.C;
            end
            directimesSorted = C(:, 1)' + startStim;

            % Use params.preBase directly — no formula needed
            preBase = params.preBase;

            % Total trial window = baseline + post-stim period
            windowTotal = preBase + params.postStim;

            % Lock in time window from first valid experiment
            if isempty(lockedPreBase)
                lockedPreBase = preBase;
                lockedEdges   = 0 : params.binWidth : windowTotal;
                lockedNBins   = numel(lockedEdges) - 1;
                tAxis         = lockedEdges(1:end-1);
                fprintf('  Locked window: preBase=%d ms, postStim=%d ms, nBins=%d\n', ...
                    lockedPreBase, params.postStim, lockedNBins);
            end

            % ----------------------------------------------------------
            % Find responsive neurons
            % ----------------------------------------------------------
            eNeurons = find(pvals < params.alpha);

            if isempty(eNeurons)
                fprintf('  [%s] No responsive neurons in exp %d.\n', stimType, ex);
                if ~isempty(psthAll{s})
                    psthAll{s} = [psthAll{s}; NaN(1, lockedNBins)];
                end
                continue
            end

            fprintf('  [%s] %d responsive neuron(s) in exp %d.\n', ...
                stimType, ex, numel(eNeurons));

            % ----------------------------------------------------------
            % Build PSTH for each responsive neuron
            % BuildBurstMatrix returns nTrials x 1 x nTimeBins
            % Window: from (trialOnset - preBase) for windowTotal ms
            % ----------------------------------------------------------
            psthRateNeurons = zeros(numel(eNeurons), lockedNBins);

            for ni = 1:numel(eNeurons)
                u = eNeurons(ni);

                % Spike matrix: rows = trials, cols = time bins (1ms each)
                MRhist = BuildBurstMatrix( ...
                    goodU(:, u), ...
                    round(p_sort.t), ...
                    round(directimesSorted - lockedPreBase), ...
                    round(windowTotal));

               

                % Remove singleton dimensions → nTrials x nTimeBins
                MRhist  = squeeze(MRhist);

                 if ~isempty(params.TakeTopPercentTrials)
                     MeanTrial = mean(MRhist,2);
                     [~, ind] = sort(MeanTrial,'descend');

                     takeTrials = ind(1:round(numel(MeanTrial)*params.TakeTopPercentTrials));

                     MRhist = MRhist(takeTrials,:);

                end
                nTrials = size(MRhist, 1);

                % Convert to spike times in ms
                spikeTimes = repmat((1:size(MRhist, 2)), nTrials, 1);
                spikeTimes = spikeTimes(logical(MRhist));

                % Bin into locked edges and convert to spk/s
                counts = histcounts(spikeTimes, lockedEdges);
                psthRateNeurons(ni, :) = (counts / (params.binWidth * nTrials)) * 1000;
            end

            % Average across responsive neurons → 1 x lockedNBins
            psthExp = mean(psthRateNeurons, 1, 'omitnan');

            if params.zScore
                baselineBins = tAxis < lockedPreBase;
                baselineMean = mean(psthExp(baselineBins));
                baselineStd  = std(psthExp(baselineBins));
                if baselineStd > 0
                    psthExp = (psthExp - baselineMean) / baselineStd;
                else
                    warning('  [%s] Baseline std is zero for exp %d — skipping experiment.', stimType, ex);
                    if ~isempty(psthAll{s})
                        psthAll{s} = [psthAll{s}; NaN(1, lockedNBins)];
                    end
                    continue  % skip to next experiment, do not append raw rates
                end
            end

            % Append as new row — guaranteed lockedNBins wide
            psthAll{s} = [psthAll{s}; psthExp(:)'];

        end % end stim loop
    end % end experiment loop

    % ------------------------------------------------------------------
    % Save results to struct
    % ------------------------------------------------------------------
    S.expList      = exList;           % experiment list for future matching
    S.lockedEdges  = lockedEdges;      % bin edges used (ms from trial start)
    S.lockedPreBase = lockedPreBase;   % baseline duration in ms
    S.params       = params;           % all parameters used

    % Save one field per stim type, named by stim e.g. S.rectGrid
    for s = 1:numel(params.stimTypes)
        stimField      = matlab.lang.makeValidName(params.stimTypes(s)); % safe field name
        S.(stimField)  = psthAll{s};   % nExp x nBins PSTH matrix
    end

    save([saveDir nameOfFile], '-struct', 'S');
    fprintf('\nSaved PSTHs to:\n  %s\n', [saveDir nameOfFile]);

else
    % ------------------------------------------------------------------
    % Load psthAll from saved struct
    % ------------------------------------------------------------------
    lockedEdges   = S.lockedEdges;
    lockedPreBase = S.lockedPreBase;

    psthAll = cell(1, numel(params.stimTypes));
    for s = 1:numel(params.stimTypes)
        stimField   = matlab.lang.makeValidName(params.stimTypes(s));
        if isfield(S, stimField)
            psthAll{s} = S.(stimField);  % load the nExp x nBins matrix
        else
            % Stim type not found in saved file — warn and leave empty
            warning('Stim type "%s" not found in saved file.', params.stimTypes(s));
            psthAll{s} = [];
        end
    end

end % end forloop

% =========================================================================
% PLOT
% =========================================================================

tAxis     = lockedEdges(1:end-1);
tAxisPlot = tAxis - lockedPreBase;

colors = lines(numel(params.stimTypes));

fig = figure;
set(fig, 'Units', 'centimeters', 'Position', [5 5 9 10]);  % single axis now

% ------------------------------------------------------------------
% Map stimulus type names to short legend labels
% ------------------------------------------------------------------
stimLegendMap = containers.Map(...
    {'linearlyMovingBall', 'rectGrid', 'MovingGrating', 'StaticGrating'}, ...
    {'MB',                 'SB',       'MG',            'SG'});

% ------------------------------------------------------------------
% First pass: compute mean/sem for all stim types and find global ylim
% ------------------------------------------------------------------
meanAll = cell(1, numel(params.stimTypes));
semAll  = cell(1, numel(params.stimTypes));
yMax    = 0;
yMin    = inf;

for s = 1:numel(params.stimTypes)
    data = psthAll{s};
    if isempty(data)
        continue
    end
    validRows  = ~all(isnan(data), 2);
    data       = data(validRows, :);
    if isempty(data)
        continue
    end
    meanAll{s} = mean(data, 1, 'omitnan');
    semAll{s}  = std(data, 0, 1, 'omitnan') / sqrt(sum(~isnan(data(:,1))));
    yMax = max(yMax, max(meanAll{s} + semAll{s}));
    yMin = min(yMin, min(meanAll{s} - semAll{s}));
end

% Y limits with 10% padding
yPad = (yMax - yMin) * 0.1;
if params.zScore
    yLims = [yMin - yPad, yMax + yPad];
else
    yLims = [max(0, yMin - yPad), yMax + yPad];
end

% ------------------------------------------------------------------
% Single axis plot — all stim types overlaid
% ------------------------------------------------------------------
ax = axes(fig);
hold(ax, 'on');

legendHandles = gobjects(numel(params.stimTypes), 1);  % store line handles for legend

for s = 1:numel(params.stimTypes)

    data = psthAll{s};
    if isempty(data)
        continue
    end
    validRows = ~all(isnan(data), 2);
    data      = data(validRows, :);
    if isempty(data)
        continue
    end

    meanPSTH = meanAll{s};
    semPSTH  = semAll{s};

    % Get short legend label for this stim type
    stimKey = char(params.stimTypes(s));
    if isKey(stimLegendMap, stimKey)
        legendLabel = stimLegendMap(stimKey);
    else
        legendLabel = stimKey;  % fallback to full name if not in map
    end

    % Shade ±SEM band
    if params.shadeSTD && size(data, 1) > 1
        upper = meanPSTH + semPSTH;
        lower = meanPSTH - semPSTH;
        xFill = [tAxisPlot(:)', fliplr(tAxisPlot(:)')];
        yFill = [upper(:)',     fliplr(lower(:)')    ];
        fill(ax, xFill, yFill, colors(s,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end

    % Mean PSTH line — store handle for legend
    legendHandles(s) = plot(ax, tAxisPlot(:)', meanPSTH(:)', ...
        'Color', colors(s,:), 'LineWidth', 1.5, 'DisplayName', legendLabel);

    % Number of contributing experiments as text
    nValid = sum(validRows);
    fprintf('  [%s] n=%d experiments in plot.\n', legendLabel, nValid);

end

% Stim onset and end of post-stim window
xline(ax, 0,               'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(ax, params.postStim, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% Y label
if params.zScore
    yLabel = 'Z-score';
else
    yLabel = '[spk/s]';
end

xlabel(ax, 'Time re. stim onset [ms]', 'FontName', 'helvetica', 'FontSize', 8);
ylabel(ax, yLabel,                      'FontName', 'helvetica', 'FontSize', 8);
xlim(ax, [tAxisPlot(1) tAxisPlot(end)]);
ylim(ax, yLims);

% Legend — only show valid handles (skip stim types with no data)
validHandles = legendHandles(isgraphics(legendHandles));
legend(validHandles, 'Location', 'northeast', 'FontName', 'helvetica', 'FontSize', 8);

ax.FontName = 'helvetica';
ax.FontSize  = 8;
hold(ax, 'off');

sgtitle(sprintf('N = %d', numel(exList)), 'FontName', 'helvetica', 'FontSize', 11);

ax = gca;
ax.YAxis.FontSize = 8;
ax.YAxis.FontName = 'helvetica';

ax = gca;
ax.XAxis.FontSize = 8;
ax.XAxis.FontName = 'helvetica';

set(fig, 'Units', 'centimeters');
set(fig, 'Position', [20 20 5 6]);

if params.PaperFig
    vs_first.printFig(fig, sprintf('PSTH-comparison-%s-%s', ...
        params.stimTypes(1), params.stimTypes(2)), PaperFig = params.PaperFig)
end

end