function plotRaster(obj, params)
% plotRaster  Combined static + drifting raster, PSTH, and raw trace.
%
% Plots both grating phases in a single raster, divided only by angle.
% TF/SF sub-divisions are intentionally omitted because the stimulus design
% does not guarantee equal representation of all (TF, SF) combinations across
% angles (see diagnostic output printed at startup).
%
% DESIGN IMBALANCE — OPTIONS FOR THE PAPER:
%   A) Common-subset: restrict quantitative analyses (tuning curves, DSI, OSI)
%      to the (TF, SF) combinations that appear in every angle. The raster
%      can still show all trials.
%   B) Velocity grouping: TF/SF (deg/s) may have been the intended constant.
%      Group by velocity rather than TF and SF independently.
%   C) Per-condition tuning: compute one orientation tuning curve per (TF,SF)
%      combination that has sufficient trials across all angles.
%
% Trial timing (from ResponseWindow):
%   |-- preBase --|-- staticDur --|-- movingDur --|-- preBase --|
%                 ^               ^               ^
%              stim on       phase change      stim off
%
% ResponseWindow fields:
%   Onsets(:,1)  = static onset  per trial (ms, absolute)
%   Onsets(:,2)  = moving onset  per trial (ms, absolute)
%   Offsets(:,2) = moving offset per trial (ms, absolute)
%   C columns    = [stimOnTime, angle, TF, SF]

% --------------------------------------------------------------------------
% BUG FIXES (all carried forward from earlier versions):
%  1. best_row uninitialized when SelectedWindow=false  -> default = 1.
%  2. [nT,nN,nB]=size(Mr2) on a 2-D matrix             -> removed.
%  3. Floating-point filter comparisons                 -> abs(a-b)<tol.
%  4. yline LineWidth=0.05 (invisible)                  -> 0.5.
%  5. ur/u dual-index confusion                         -> documented.
%  6. Overlapping dividing lines at angle boundaries    -> loop from row 2.
%  7. Wrong y-tick formula from moving-ball code        -> midpoint formula.
%  8. Angle-block boundaries now robust to imbalanced designs.
%  9. Diagnostic condition-balance table printed at startup.
%
% FIXED IN THIS VERSION:
% 10. Red patch x-offset bug: 'start' is already in full-window coordinates
%     (0 = first bin of the preBase buffer), so adding another preBase
%     shifted the patch preBase/params.bin bins to the right of the actual
%     best window. Fixed: patch drawn at [start, start+window]/params.bin.
% 11. Raw-trace xlines: the inherited 'xline(-start/1000)' marked a position
%     before the visible window. Replaced with correctly computed event
%     markers (stim on, phase transition, stim off) that are only drawn when
%     they fall inside the 500 ms trace window.
% --------------------------------------------------------------------------

arguments (Input)
    obj
    params.overwrite            logical = false            % Overwrite saved figures
    params.analysisTime                 = datetime('now')  % Provenance timestamp
    params.inputParams          logical = false            % Print params and exit
    params.preBase                      = 200              % Pre/post-stimulus baseline (ms)
    params.bin                          = 30               % Raster bin size (ms/bin)
    params.exNeurons                    = []               % Neuron index into good-unit list
    params.exNeuronsPhyID       double  = []               % Override: select by phy cluster ID
    params.AllSomaticNeurons    logical = false            % Plot all good units
    params.AllResponsiveNeurons logical = true             % Plot units with min(pS,pM) < 0.05
    params.SelectedWindow       logical = true             % Auto-detect best response window
    params.MergeNtrials                 = 1                % Trials averaged per raster row
    params.GaussianLength               = 10               % Gaussian smoothing kernel (bins)
    params.Gaussian             logical = false            % Apply Gaussian smoothing
    params.MaxVal_1             logical = true             % Clamp raster colormap to [0 1]
    params.OneAngle             string  = "all"            % Restrict to one angle, e.g. "90"
    params.OneTF                string  = "all"            % Restrict to one TF (Hz)
    params.OneSF                string  = "all"            % Restrict to one SF (c/deg)
    params.PaperFig             logical = false            % High-quality figure export
    params.statType             string  = "maxPermuteTest" % "maxPermuteTest"|"BootstrapPerNeuron"
    params.plotRaw               logical = true                
end

if params.inputParams, disp(params); return; end

% ==========================================================================
% 1. LOAD PRE-COMPUTED RESULTS
% ==========================================================================

NeuronResp = obj.ResponseWindow;

if params.statType == "BootstrapPerNeuron"
    Stats = obj.BootstrapPerNeuron;
else
    Stats = obj.StatisticsPerNeuron;
end

pvalsS   = Stats.Static.pvalsResponse;
pvalsM   = Stats.Moving.pvalsResponse;
pvalsMin = min(pvalsS, pvalsM);

% ==========================================================================
% 2. SPIKE SORTING AND STIMULUS TIMING
% ==========================================================================

p       = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
phy_IDg = p.phy_ID(string(p.label') == 'good');
label   = string(p.label');
goodU   = p.ic(:, label == 'good');

staticDur    = round(mean(NeuronResp.Onsets(:,2)  - NeuronResp.Onsets(:,1))); % (ms)
totalStimDur = round(mean(NeuronResp.Offsets(:,2) - NeuronResp.Onsets(:,1))); % (ms)
movingDur    = totalStimDur - staticDur;                                        % (ms)
stimInter    = NeuronResp.stimInter;
preBase      = round(stimInter - stimInter / 4);   % 75% of ITI (ms)

% ==========================================================================
% 3. CONDITION MATRIX C AND OPTIONAL FILTERING
% ==========================================================================

C = NeuronResp.C;   % columns: [stimOnTime, angle, TF, SF]

if params.OneAngle ~= "all"
    angleVal = str2double(params.OneAngle);
    C = C(abs(C(:,2) - angleVal) < 1e-3, :);
    if isempty(C), error('No trials found for OneAngle = "%s"', params.OneAngle); end
end
if params.OneTF ~= "all"
    tfVal = str2double(params.OneTF);
    C = C(abs(C(:,3) - tfVal) < 1e-4, :);
    if isempty(C), error('No trials found for OneTF = "%s"', params.OneTF); end
end
if params.OneSF ~= "all"
    sfVal = str2double(params.OneSF);
    C = C(abs(C(:,4) - sfVal) < 1e-4, :);
    if isempty(C), error('No trials found for OneSF = "%s"', params.OneSF); end
end

% Sort angle primary, TF secondary, SF tertiary so that similar conditions
% remain adjacent within each angle block even without TF/SF boundary lines
[C, ~] = sortrows(C, [2 3 4]);

% Static-onset time for each sorted trial — the single reference time for
% all BuildBurstMatrix calls (both phases are captured relative to this)
directimesSorted = C(:,1)';

% ==========================================================================
% 4. CONDITION COUNTS AND DIAGNOSTIC TABLE
% ==========================================================================

uAngle = unique(C(:,2));
uTF    = unique(C(:,3));
uSF    = unique(C(:,4));
angleN = numel(uAngle);
tfN    = numel(uTF);
sfN    = numel(uSF);
nT     = size(C, 1);

fprintf('\n=== Condition balance check ===\n');
fprintf('%-8s  %-8s  %-8s  %-8s  %s\n', 'Angle', 'TF [Hz]', 'SF [c/d]', 'Vel[d/s]', 'nTrials');
fprintf('%s\n', repmat('-', 1, 52));
for a = 1:angleN
    for t = 1:tfN
        for s = 1:sfN
            mask = abs(C(:,2)-uAngle(a))<1e-3 & ...
                   abs(C(:,3)-uTF(t))<1e-4    & ...
                   abs(C(:,4)-uSF(s))<1e-4;
            if sum(mask) > 0
                fprintf('%-8.0f  %-8.2f  %-8.4f  %-8.1f  %d\n', ...
                    uAngle(a), uTF(t), uSF(s), uTF(t)/uSF(s), sum(mask));
            end
        end
    end
end
fprintf('%s\n\n', repmat('=', 1, 52));

% --- Angle-block boundaries (robust to imbalanced designs) ---
% Scanning sorted C for angle transitions avoids the assumption that every
% angle has the same number of trials.
%
% angleChangeIdx: (angleN+1)-vector where:
%   angleChangeIdx(a)   = first row of angle block a  (1-based into C)
%   angleChangeIdx(a+1) = first row of angle block a+1, or nT+1 for the last
angleChangeIdx  = [find([true; diff(C(:,2)) ~= 0]); nT + 1];
nTrialsPerAngle = diff(angleChangeIdx);                             % [angleN x 1]
angleMidpoints  = angleChangeIdx(1:end-1) + (nTrialsPerAngle-1)/2; % midpoint rows

% ==========================================================================
% 5. NEURON SELECTION
% ==========================================================================

if ~isempty(params.exNeuronsPhyID)
    [found, neuronIdx] = ismember(params.exNeuronsPhyID, phy_IDg);
    if any(~found)
        warning('Phy IDs not found in good units (skipped): %s', ...
            num2str(params.exNeuronsPhyID(~found)));
    end
    params.exNeurons = neuronIdx(found);
    fprintf('  Converted phy IDs [%s] -> unit indices [%s]\n', ...
        num2str(params.exNeuronsPhyID(found)), num2str(params.exNeurons));
end

if params.AllSomaticNeurons
    eNeuron = 1:size(goodU, 2);
elseif params.AllResponsiveNeurons
    eNeuron = find(pvalsMin < 0.05);
    if isempty(eNeuron)
        fprintf('No responsive neurons (min p < 0.05 across both phases).\n'); return
    end
else
    eNeuron = params.exNeurons;
end

% ==========================================================================
% 6. BUILD RASTER MATRIX
% ==========================================================================

% Mr: [nTrials x nSelectedNeurons x nBins]
% Window = preBase + staticDur + movingDur + preBase, in params.bin ms bins.
% Column 1 of Mr = preBase ms before static onset (= beginning of window).
% Column preBase/params.bin of Mr = static onset.
% Column (preBase+staticDur)/params.bin of Mr = moving onset.
Mr = BuildBurstMatrix( ...
    goodU(:, eNeuron), ...
    round(p.t / params.bin), ...
    round((directimesSorted - preBase) / params.bin), ...
    round((totalStimDur + preBase*2) / params.bin));

if params.Gaussian
    Mr = ConvBurstMatrix(Mr, fspecial('gaussian', [1 params.GaussianLength], 3), 'same');
end

channels      = goodU(1, eNeuron);
[~, ~, nBins] = size(Mr);

% ==========================================================================
% 7. PER-NEURON FIGURE LOOP
% ==========================================================================
%
% INDEXING:
%   u  = absolute index into goodU   (goodU(:,u), phy_IDg(u), pvalsS(u))
%   ur = relative index into eNeuron (Mr(:,ur,:), channels(ur))
%   ur must be incremented in every code path.

ur = 1;

for u = eNeuron

    fig = figure;

    % ------------------------------------------------------------------
    % 7a. 2-D merged raster [nTrials x nBins]
    % ------------------------------------------------------------------

    mergeTrials = params.MergeNtrials;
    Mr2         = zeros(nT, nBins);

    if mergeTrials > 1
        for i = 1:mergeTrials:nT
            meanb = mean(squeeze(Mr(i:min(i+mergeTrials-1, end), ur, :)), 1);
            Mr2(i:i+mergeTrials-1, :) = repmat(meanb, [mergeTrials 1]);
        end
    else
        Mr2 = squeeze(Mr(:, ur, :));   % [nTrials x nBins]
    end

    if sum(Mr2, 'all') == 0
        close(fig); ur = ur + 1; continue
    end

    % ==================================================================
    % PANEL 2 (rows 6-16): Combined raster
    % ==================================================================

    ax_raster = subplot(18, 1, [6 14]);

    imagesc(Mr2 .* (1000 / params.bin));   % Display in spk/s
    colormap(flipud(gray(64)));
    hold on;

    % --- Vertical time markers ---
    % x = 0 corresponds to the start of the preBase buffer.
    % x = preBase/params.bin corresponds to the static onset.
    xMax = round(totalStimDur + preBase*2) / params.bin;
    xline(preBase / params.bin,                  'k',   'LineWidth', 1.5); % Stim on
    xline((preBase + staticDur) / params.bin,    '--k', 'LineWidth', 1.2); % Phase transition
    xline((preBase + totalStimDur) / params.bin, 'k',   'LineWidth', 1.5); % Stim off

    if params.MaxVal_1, caxis([0 1]); end

    % --- Angle-boundary horizontal lines ---
    for a = 2:angleN
        yline(angleChangeIdx(a) - 0.5, 'k', 'LineWidth', 2);
    end

    % Phase labels above the raster
    yAbove = nT + nT * 0.07;
    text(preBase/params.bin + (staticDur/params.bin)/2, yAbove, 'Static', ...
        'HorizontalAlignment', 'center', 'FontSize', 7, ...
        'FontName', 'helvetica', 'Clipping', 'off');
    text((preBase+staticDur)/params.bin + (movingDur/params.bin)/2, yAbove, 'Moving', ...
        'HorizontalAlignment', 'center', 'FontSize', 7, ...
        'FontName', 'helvetica', 'Clipping', 'off');

    xlim([0, xMax]);
    xticks([0, preBase/params.bin : 600/params.bin : xMax, ...
        round((totalStimDur+preBase*2)/100)*100 / params.bin]);
    xticklabels([]);

    % Y-ticks: one per angle block, at block midpoint, labeled with trial count
    if numel(uAngle) ~=1
        yticks(angleMidpoints);
        yticklabels(arrayfun(@num2str, nTrialsPerAngle, 'UniformOutput', false));
    end

    ax_raster.YAxis.FontSize = 8;
    ax_raster.YAxis.FontName = 'helvetica';
    ylabel('Trials', 'FontSize', 10, 'FontName', 'helvetica');

    % ==================================================================
    % RIGHT-SIDE ANGLE LABELS
    % ==================================================================

    colGap  = xMax * 0.06;
    tickLen = colGap * 0.40;
    xCol    = xMax + colGap * 0.5;

    text(xCol + tickLen, 0, 'Angle', ...
        'FontSize', 5.5, 'FontWeight', 'bold', 'FontName', 'helvetica', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Clipping', 'off');

    for a = 1:angleN
        rowStart = angleChangeIdx(a);
        rowEnd   = angleChangeIdx(a+1) - 1;
        yMid     = angleMidpoints(a);
        yT       = rowStart - 0.5;
        yB       = rowEnd   + 0.5;

        line([xCol, xCol + tickLen], [yT, yT], ...
            'Color', [0.4 0.4 0.4], 'LineWidth', 0.5, 'Clipping', 'off');
        line([xCol, xCol],           [yT, yB], ...
            'Color', [0.4 0.4 0.4], 'LineWidth', 0.5, 'Clipping', 'off');
        line([xCol, xCol + tickLen], [yB, yB], ...
            'Color', [0.4 0.4 0.4], 'LineWidth', 0.5, 'Clipping', 'off');

        text(xCol + tickLen * 1.5, yMid, sprintf('%.0f°', uAngle(a)), ...
            'FontSize', 7, 'FontWeight', 'bold', 'FontName', 'helvetica', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Clipping', 'off');
    end

    % ==================================================================
    % 7b. Identify best response window
    % ==================================================================

    if params.SelectedWindow

        % Find the angle block with the highest mean firing rate
        meanMr = zeros(1, angleN);
        for a = 1:angleN
            meanMr(a) = mean(Mr2(angleChangeIdx(a):angleChangeIdx(a+1)-1, :), 'all');
        end
        [~, bestAngle_a] = max(meanMr);

        % All trial indices for the best angle block
        trials = angleChangeIdx(bestAngle_a) : angleChangeIdx(bestAngle_a+1) - 1;

        window = 500;   % Sliding-window width (ms)

        % Mr2 rows for the best angle block [nTrialsInBlock x nBins]
        X = Mr2(trials, :);
        X(X > 1) = 1;   % Clip outliers before the search

        [n_rows, n_cols] = size(X);
        nWinPos          = n_cols - round(window / params.bin) + 1;

        % Per-trial mean inside every sliding window [n_rows x nWinPos]
        window_means = zeros(n_rows, nWinPos);
        for col = 1:nWinPos
            window_means(:, col) = mean(X(:, col : col+round(window/params.bin)-1), 2);
        end

        % Best (trial, window-start) pair
        [~, linear_idx]      = max(window_means(:));
        [best_row, best_col] = ind2sub(size(window_means), linear_idx);

        % 'start': window start position in FULL-WINDOW coordinates.
        % This is in the same coordinate system as Mr2 columns:
        %   start = 0            -> very beginning of the preBase buffer
        %   start = preBase      -> static onset
        %   start = preBase+staticDur -> moving onset
        % (i.e., preBase is already included in 'start')
        start = best_col * params.bin;   % (ms, from beginning of recording window)

        bestPhase = 'Static';
        if start >= preBase + staticDur
            bestPhase = 'Moving';
        elseif start >= preBase
            bestPhase = 'Static (stim)';
        end

    else
        [~, bestPhaseIdx] = max([ ...
            max(NeuronResp.Static.NeuronVals(u, :, 4)), ...
            max(NeuronResp.Moving.NeuronVals(u, :, 4))]);
        phaseNames     = ["Static", "Moving"];
        bestPhase      = phaseNames(bestPhaseIdx);
        [~, maxRespIn] = max(NeuronResp.(bestPhase).NeuronVals(u, :, 4));
        start          = NeuronResp.(bestPhase).NeuronVals(u, maxRespIn, 3) * ...
                         NeuronResp.params.binRaster - 20;
        window         = 500;
        bestAngleVal   = NeuronResp.(bestPhase).NeuronVals(u, maxRespIn, 6);
        bestAngle_a    = find(abs(uAngle - bestAngleVal) < 1e-3, 1);
        if isempty(bestAngle_a), bestAngle_a = 1; end
        trials         = angleChangeIdx(bestAngle_a) : angleChangeIdx(bestAngle_a+1) - 1;
        best_row       = 1;   % FIX Bug 1: conservative default
    end

    RasterTrials = trials(best_row);   % Absolute trial index of the best single trial

    % Grey band: full time extent of the best angle block
    if numel(uAngle) ~=1
        yBandTop = angleChangeIdx(bestAngle_a) - 0.5;
        yBandBot = angleChangeIdx(bestAngle_a+1) - 0.5;
        patch([0, xMax, xMax, 0], [yBandTop, yBandTop, yBandBot, yBandBot], ...
            'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end

    % --- Red patch: best trial x best response window ---
    %
    % FIX (Bug 10): 'start' is already in full-window coordinates (column
    % index of Mr2 * params.bin), so the raster x-positions are simply
    % start/params.bin and (start+window)/params.bin.
    %
    % The previous code used (preBase+start)/params.bin, which added an
    % extra preBase offset and shifted the patch to the right by
    % preBase/params.bin bins — misaligning it with the raw trace below.

    patch([start/params.bin,        (start+window)/params.bin, ...
           (start+window)/params.bin, start/params.bin], ...
          [RasterTrials-0.5, RasterTrials-0.5, RasterTrials+0.5, RasterTrials+0.5], ...
          'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % ==================================================================
    % PANEL 3 (rows 17-18): PSTH for the best angle block
    % ==================================================================

    ax_psth = subplot(18, 1, [16 18]);

    MRhist = BuildBurstMatrix(goodU(:, u), round(p.t), ...
        round(directimesSorted - preBase), round(totalStimDur + preBase*2));
    MRhist = squeeze(MRhist(trials, :, :));   % [nTrialsInBestAngle x windowDur_ms]

    [nT2, nB2] = size(MRhist);
    spikeTimes  = repmat(1:nB2, nT2, 1);
    spikeTimes  = spikeTimes(logical(MRhist));

    binWidth = 125;
    if nBins > 300, binWidth = 250; end

    edges      = 1:binWidth:round(totalStimDur + preBase*2);
    psthCounts = histcounts(spikeTimes, edges);
    psthRate   = (psthCounts / (binWidth * nT2)) * 1000;   % [spk/s]

    b = bar(edges(1:end-1), psthRate, 'histc');
    b.FaceColor = 'k';  b.FaceAlpha = 0.3;  b.MarkerEdgeColor = 'none';

    xlim([0, round((totalStimDur + preBase*2) / 100) * 100]);

    try
        ylim([0, max(psthRate) + std(psthRate)]);
    catch
        close(fig); ur = ur + 1; continue
    end

    xticks([0, preBase:600:(totalStimDur+preBase*2), ...
        round((totalStimDur+preBase*2)/100)*100]);
    xline(preBase,                   'LineWidth', 1.5);
    xline(preBase + staticDur, '--', 'LineWidth', 1.2);
    xline(preBase + totalStimDur,    'LineWidth', 1.5);
    xticklabels([-(preBase), 0:600:round((totalStimDur/100))*100, ...
        round((totalStimDur/100))*100 + 2*preBase] ./ 1000);

    ax_psth.XAxis.FontSize = 8;  ax_psth.XAxis.FontName = 'helvetica';
    ax_psth.YAxis.FontSize = 8;  ax_psth.YAxis.FontName = 'helvetica';
    ylabel('[spk/s]', 'FontSize', 10, 'FontName', 'helvetica');
    xlabel('Time [s]',  'FontSize', 10, 'FontName', 'helvetica');

    ylims = ylim;
    yticks([round(ylims(2)/10)*5, ceil(ylims(2)/10)*10]);

    % ==================================================================
    % PANEL 1 (rows 1-3): Raw AP/LFP trace for the best single trial
    % ==================================================================

   if params.plotRaw   
    chan = goodU(1, u);

    subplot(18, 1, [1 3]);

    % --- Raw trace start time ---
    %
    % 'start' is in full-window ms (0 = beginning of preBase buffer, i.e.
    % preBase ms before static onset).
    %
    % Absolute time of static onset for the best trial:
    %   staticOnset_best = directimesSorted(RasterTrials)
    %
    % The window starts at:
    %   staticOnset_best - preBase  (the recording window start)
    %
    % The best response window starts 'start' ms into that recording window:
    %   startTimes = (staticOnset_best - preBase) + start
    %             = staticOnset_best + start - preBase
    %
    % This is the same absolute time the red patch begins at on the raster,
    % so the raw trace and the patch are now guaranteed to be aligned.

    startTimes = directimesSorted(RasterTrials) + start - preBase;   % (ms, absolute)

    spikes = squeeze(BuildBurstMatrix(goodU(:,u), round(p.t), round(startTimes), round(window)));

    [fig, ~, ~] = PlotRawDataNP(obj, fig=fig, chan=chan, ...
        startTimes=startTimes, window=window, spikeTimes=spikes);

    ax_raw = gca;
    ax_raw.YAxis.FontSize = 8;  ax_raw.YAxis.FontName = 'helvetica';
    ax_raw.XAxis.FontSize = 8;  ax_raw.XAxis.FontName = 'helvetica';
    ax_raw.XRuler.TickDirection = 'out';
    ax_raw.XAxisLocation = 'bottom';

    xlims = xlim;
    xticks(0:(xlims(2)/5):xlims(2));
    xticklabels(0:100:window);

    % --- Stimulus-event markers in the raw trace (FIX Bug 11) ---
    %
    % The raw trace runs from absolute time 'startTimes' to 'startTimes+window'.
    % The x-axis is mapped linearly to [0, window] ms by the tick labels above.
    % Scale factor converts ms offsets to x-axis units.
    %
    % Event offsets from the trace start (ms):
    %   static onset     : directimesSorted(RasterTrials) - startTimes
    %                    = directimesSorted(RasterTrials) - (directimesSorted(RasterTrials) + start - preBase)
    %                    = preBase - start
    %   phase transition : preBase - start + staticDur
    %   stim offset      : preBase - start + totalStimDur
    %
    % Only draw events that fall inside [0, window] ms.
    %
    % (The original 'xline(-start/1000)' divided ms by 1000 suggesting a
    % seconds axis but the tick labels are in ms — it reliably produced a
    % marker at ~0 ms regardless of 'start', which was incorrect.)

    scale = xlims(2) / window;   % x-axis units per ms (handles both ms and s axes)

    evOffsets = [preBase - start, ...                  % static onset
                 preBase - start + staticDur, ...      % static->moving transition
                 preBase - start + totalStimDur];      % stim offset
    evStyles  = {'k',  '--k', 'k' };
    evWidths  = [1.5,   1.2,  1.5 ];

    for ev = 1:3
        xEv = evOffsets(ev) * scale;
        if xEv >= xlims(1) && xEv <= xlims(2)
            xline(xEv, evStyles{ev}, 'LineWidth', evWidths(ev));
        end
    end

    xlabel('Time [ms]', 'FontName', 'helvetica', 'FontSize', 10);
    ylabel('[\muV]',    'FontSize', 10,           'FontName', 'helvetica');

   end

    % Title: identity + best angle + per-phase p-values
    bestAngleVal = C(RasterTrials, 2);
    bestTF       = C(RasterTrials, 3);
    bestSF       = C(RasterTrials, 4);

    title({ ...
        sprintf('U.%d  Chan-%d  Phy-%d  |  pS=%.4f  pM=%.4f', ...
            u, channels(ur), phy_IDg(u), pvalsS(u), pvalsM(u)), ...
        sprintf('Best angle: %.0f°  (TF=%.1f Hz, SF=%.3f c/d at best trial)  [%s]', ...
            bestAngleVal, bestTF, bestSF, bestPhase) ...
        });

    % ==================================================================
    % AXES POSITION ADJUSTMENT
    % ==================================================================
    % Shrink raster and PSTH by the same factor to keep their time axes
    % aligned while making room for the angle-label column on the right.

    shrinkFactor = 0.85;
    pos = ax_raster.Position;
    ax_raster.Position = [pos(1), pos(2), pos(3)*shrinkFactor, pos(4)];
    pos = ax_psth.Position;
    ax_psth.Position   = [pos(1), pos(2), pos(3)*shrinkFactor, pos(4)];

    % ==================================================================
    % 7c. Figure layout and export
    % ==================================================================

    set(fig, 'Units', 'centimeters');
    set(fig, 'Position', [20 20 9 12]);

    if params.PaperFig
        obj.printFig(fig, sprintf('%s-Grating-CombinedRaster-Unit%d', ...
            obj.dataObj.recordingName, u), 'PaperFig', params.PaperFig);
    elseif params.overwrite
        obj.printFig(fig, sprintf('%s-Grating-CombinedRaster-Unit%d', ...
            obj.dataObj.recordingName, u));
    end

    %if ur ~= length(eNeuron), close(fig); end

    ur = ur + 1;   % MUST be reached in every code path above

end  % end neuron loop

end  % end plotRaster