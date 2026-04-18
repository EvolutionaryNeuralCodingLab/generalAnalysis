function plotRaster(obj, params)
% plotRaster  Combined static + drifting raster, PSTH, and raw trace.
%
% Plots both grating phases (static and moving) in a single raster for each
% neuron. A solid vertical line marks stimulus onset and offset; a dashed
% vertical line marks the static→moving phase transition.
%
% Trial timing (from ResponseWindow):
%   |-- preBase --|-- staticDur --|-- movingDur --|-- preBase --|
%                 ^               ^               ^
%              stim on       phase change      stim off
%
% ResponseWindow stores:
%   NeuronResp.Onsets(:,1)  = static onset  per trial (ms)
%   NeuronResp.Onsets(:,2)  = moving onset  per trial (ms)
%   NeuronResp.Offsets(:,2) = moving offset per trial (ms)
%   NeuronResp.C  columns   = [stimOnTime, angle, TF, SF]
%
% Usage:
%   obj.plotRaster()
%   obj.plotRaster('AllResponsiveNeurons', true)
%   obj.plotRaster('exNeuronsPhyID', [42 87], 'PaperFig', true)

% --------------------------------------------------------------------------
% BUG FIXES vs original moving-ball plotRaster (carried forward from v1):
%  1. best_row uninitialized when SelectedWindow=false  -> default = 1.
%  2. [nT,nN,nB]=size(Mr2) on a 2-D matrix             -> removed.
%  3. Floating-point filter comparisons use tolerance   -> abs(a-b)<eps.
%  4. yline LineWidth=0.05 (invisible)                  -> 0.5.
%  5. ur/u dual-index confusion                         -> fully documented.
%
% NEW in combined version:
%  6. stimType removed: both phases always shown together.
%  7. Responsive-neuron selection uses min(pStatic, pMoving) to catch
%     neurons that respond to only one phase.
%  8. Three xline markers instead of two (onset, transition, offset).
%  9. Title reports p-values for both phases independently.
% --------------------------------------------------------------------------

arguments (Input)
    obj
    params.overwrite            logical = false           % Overwrite saved figures
    params.analysisTime                 = datetime('now') % Provenance timestamp
    params.inputParams          logical = false           % Print params and exit
    params.preBase                      = 200             % Pre/post-stimulus baseline (ms)
    params.bin                          = 30              % Raster bin size (ms/bin)
    params.exNeurons                    = 1               % Neuron index into good-unit list
    params.exNeuronsPhyID       double  = []              % Override: select by phy cluster ID
    params.AllSomaticNeurons    logical = false           % Plot all good units
    params.AllResponsiveNeurons logical = true           % Plot units with min(pS,pM) < 0.05
    params.SelectedWindow       logical = true            % Auto-detect best response window
    params.MergeNtrials                 = 1               % Trials averaged per raster row
    params.GaussianLength               = 10              % Gaussian smoothing kernel (bins)
    params.Gaussian             logical = false           % Apply Gaussian smoothing
    params.MaxVal_1             logical = true            % Clamp raster colormap to [0 1]
    params.OneAngle             string  = "all"           % Restrict to one angle, e.g. "90"
    params.OneTF                string  = "all"           % Restrict to one TF (Hz)
    params.OneSF                string  = "all"           % Restrict to one SF (c/deg)
    params.PaperFig             logical = false           % High-quality figure export
    params.statType             string  = "maxPermuteTest"% "maxPermuteTest"|"BootstrapPerNeuron"
end

if params.inputParams, disp(params); return; end

% ==========================================================================
% 1. LOAD PRE-COMPUTED RESULTS
% ==========================================================================

% ResponseWindow struct: fields Static, Moving, C, Onsets, Offsets, stimInter, params
NeuronResp = obj.ResponseWindow;

% Load statistics for both phases independently
if params.statType == "BootstrapPerNeuron"
    Stats = obj.BootstrapPerNeuron;
else
    Stats = obj.StatisticsPerNeuron;
end

% Per-neuron p-values for each phase (vectors, length = nGoodUnits)
pvalsS = Stats.Static.pvalsResponse;   % p-value: static phase response
pvalsM = Stats.Moving.pvalsResponse;   % p-value: moving phase response

% For neuron selection: use the more responsive of the two phases
pvalsMin = min(pvalsS, pvalsM);        % Element-wise minimum across phases

% ==========================================================================
% 2. SPIKE SORTING AND STIMULUS TIMING
% ==========================================================================

% Load phy/Kilosort output: struct with fields ic, t, label, phy_ID
p = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);

% Phy cluster IDs restricted to good (somatic) units
phy_IDg = p.phy_ID(string(p.label') == 'good');

% Unit quality labels as string array
label = string(p.label');

% Spike train matrix: rows = channels, cols = time samples; good units only
goodU = p.ic(:, label == 'good');

% --- Derive combined trial timing from ResponseWindow ---

% staticDur: duration of the static grating phase (ms), same for all trials
% Onsets(:,1) = static onset; Onsets(:,2) = moving onset
staticDur = round(mean(NeuronResp.Onsets(:,2) - NeuronResp.Onsets(:,1)));   % (ms)

% totalStimDur: full stimulus duration, static + moving (ms)
% Offsets(:,2) = moving offset (end of stimulus)
totalStimDur = round(mean(NeuronResp.Offsets(:,2) - NeuronResp.Onsets(:,1))); % (ms)

% movingDur: duration of the drifting phase alone (ms)
movingDur = totalStimDur - staticDur;   % (ms)

% Inter-trial interval (ms), used to set the pre-stimulus baseline window
stimInter = NeuronResp.stimInter;

% Pre-stimulus baseline window: 75% of the ITI to leave a guard band
preBase = round(stimInter - stimInter / 4);  % (ms)

% ==========================================================================
% 3. CONDITION MATRIX C
% ==========================================================================

% C is saved pre-sorted in ResponseWindow: columns = [stimOnTime, angle, TF, SF]
% stimOnTime here = static onset time (Onsets(:,1)) for each trial
C = NeuronResp.C;

% Optionally restrict to one grating angle (direction, degrees)
if params.OneAngle ~= "all"
    angleVal = str2double(params.OneAngle);
    C = C(abs(C(:,2) - angleVal) < 1e-3, :);   % Tolerance-based equality
    if isempty(C)
        error('No trials found for OneAngle = "%s"', params.OneAngle);
    end
end

% Optionally restrict to one temporal frequency (Hz)
if params.OneTF ~= "all"
    tfVal = str2double(params.OneTF);
    C = C(abs(C(:,3) - tfVal) < 1e-4, :);
    if isempty(C)
        error('No trials found for OneTF = "%s"', params.OneTF);
    end
end

% Optionally restrict to one spatial frequency (cyc/deg)
if params.OneSF ~= "all"
    sfVal = str2double(params.OneSF);
    C = C(abs(C(:,4) - sfVal) < 1e-4, :);
    if isempty(C)
        error('No trials found for OneSF = "%s"', params.OneSF);
    end
end

% Re-sort after any filtering: angle primary, TF secondary, SF tertiary
[C, ~] = sortrows(C, [2 3 4]);

% Row vector of static-onset times for each (sorted) trial (ms, absolute recording time)
% This is the reference time used for all matrix constructions below
directimesSorted = C(:,1)';

% ==========================================================================
% 4. CONDITION COUNTS
% ==========================================================================

uAngle = unique(C(:,2));   % Unique grating angles (deg)
uTF    = unique(C(:,3));   % Unique temporal frequencies (Hz)
uSF    = unique(C(:,4));   % Unique spatial frequencies (cyc/deg)

angleN = numel(uAngle);    % Number of unique angles
tfN    = numel(uTF);       % Number of unique TFs
sfN    = numel(uSF);       % Number of unique SFs
nT     = size(C, 1);       % Total number of trials

% Number of repeats per unique (angle x TF x SF) combination
trialDivision = nT / (angleN * tfN * sfN);
if mod(trialDivision, 1) ~= 0
    warning('trialDivision is non-integer (%.2f): conditions may be unbalanced.', ...
        trialDivision);
    trialDivision = floor(trialDivision);
end

% ==========================================================================
% 5. RESOLVE NEURON SELECTION
% ==========================================================================

% If phy cluster IDs are provided, convert to indices into goodU.
% This overrides params.exNeurons; phy IDs are stable across re-sorts.
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

% Build eNeuron: set of absolute indices into goodU to iterate over
if params.AllSomaticNeurons
    eNeuron  = 1:size(goodU, 2);            % All good units
    pvalsOut = [eNeuron; pvalsMin(eNeuron)'];
elseif params.AllResponsiveNeurons
    % Include neurons responsive in at least one phase (min p-value < 0.05)
    eNeuron  = find(pvalsMin < 0.05);
    pvalsOut = [eNeuron; pvalsMin(eNeuron)'];
    if isempty(eNeuron)
        fprintf('No responsive neurons found (min p < 0.05 across both phases).\n');
        return
    end
else
    eNeuron  = params.exNeurons;
    pvalsOut = [eNeuron; pvalsMin(eNeuron)'];
end

% ==========================================================================
% 6. BUILD RASTER MATRIX (combined window: preBase + staticDur + movingDur + preBase)
% ==========================================================================

% Mr: [nTrials x nSelectedNeurons x nBins] using params.bin ms bins.
% The window starts preBase ms before the static onset and ends preBase ms
% after the moving offset, so both phases are captured in a single matrix.
Mr = BuildBurstMatrix( ...
    goodU(:, eNeuron), ...                                          % Selected neurons
    round(p.t / params.bin), ...                                    % Spike times in bins
    round((directimesSorted - preBase) / params.bin), ...           % Window start (bins)
    round((totalStimDur + preBase*2) / params.bin));                % Window length (bins)

% Optionally smooth raster across time with a 1-D Gaussian kernel
if params.Gaussian
    Mr = ConvBurstMatrix(Mr, fspecial('gaussian', [1 params.GaussianLength], 3), 'same');
end

% Recording channel for each selected neuron (used in title and raw-trace plot)
channels = goodU(1, eNeuron);

% Total number of time bins in the combined window (for x-axis scaling)
[~, ~, nBins] = size(Mr);

% ==========================================================================
% 7. PER-NEURON FIGURE LOOP
% ==========================================================================
%
% INDEXING:
%   u  = absolute index into goodU   (e.g., goodU(:, u), phy_IDg(u))
%   ur = relative index into eNeuron (1, 2, 3, ...) -> Mr(:, ur, :), channels(ur)
%   Both must be kept in sync; ur is incremented at the end of every branch.

ur = 1;  % Relative index into eNeuron; reset once, incremented every iteration

for u = eNeuron   % u: absolute index into goodU

    fig = figure;

    % ------------------------------------------------------------------
    % 7a. Build 2-D merged raster [nTrials x nBins] for this neuron
    % ------------------------------------------------------------------

    mergeTrials = params.MergeNtrials;
    Mr2 = zeros(nT, nBins);  % Pre-allocate: rows = trials, cols = time bins

    if mergeTrials > 1
        % Average groups of mergeTrials rows; replicate result for display continuity
        for i = 1:mergeTrials:nT
            meanb = mean(squeeze(Mr(i:min(i+mergeTrials-1, end), ur, :)), 1);
            Mr2(i:i+mergeTrials-1, :) = repmat(meanb, [mergeTrials 1]);
        end
    else
        % No merging: squeeze the neuron dimension (ur indexes the selected subset)
        Mr2 = squeeze(Mr(:, ur, :));   % [nTrials x nBins]
    end

    % Skip neuron if it has no spikes in the entire window
    if sum(Mr2, 'all') == 0
        close(fig);
        ur = ur + 1;
        continue
    end

    % ==================================================================
    % PANEL 2 (rows 6-16): Combined raster
    % ==================================================================

    subplot(18, 1, [6 16]);

    % Convert spike counts/bin to approximate spike rate (spk/s) for display
    imagesc(Mr2 .* (1000 / params.bin));
    colormap(flipud(gray(64)));   % Dark pixels = high firing rate

    hold on;

    % --- Three vertical markers (all with identical style for clean figures) ---

    % Stimulus onset (static phase starts)
    xline(preBase / params.bin, 'k', 'LineWidth', 1.5);

    % Phase transition: static -> moving (dashed, thinner)
    xline((preBase + staticDur) / params.bin, '--k', 'LineWidth', 1.2);

    % Stimulus offset (moving phase ends)
    xline((preBase + totalStimDur) / params.bin, 'k', 'LineWidth', 1.5);

    % Clamp colormap to [0 1] for cross-neuron comparability
    if params.MaxVal_1
        caxis([0 1]);
    end

    % --- Horizontal dividing lines between condition blocks ---
    % Thick black = angle boundary, thin black = TF boundary, dashed red = SF boundary
    angleStart = C(1, 2);
    tfStart    = C(1, 3);
    sfStart    = C(1, 4);

    for t = 1:nT
        if angleStart ~= C(t, 2)
            yline(t - 0.5, 'k', 'LineWidth', 2);    % New angle block
            angleStart = C(t, 2);
        end
        if tfStart ~= C(t, 3)
            yline(t - 0.5, 'k', 'LineWidth', 0.5);  % New TF block
            tfStart = C(t, 3);
        end
        if sfStart ~= C(t, 4)
            yline(t - 0.5, '--r', 'LineWidth', 0.5); % New SF block (FIX: was 0.05, invisible)
            sfStart = C(t, 4);
        end
    end

    % Phase labels above the raster (text annotations inside the axes)
    % Position them at the horizontal midpoints of each phase
    ax0 = gca;
    yTop = nT + nT*0.02;  % Just above top row

    text(preBase/params.bin + (staticDur/params.bin)/2, yTop, 'Static', ...
        'HorizontalAlignment', 'center', 'FontSize', 7, 'FontName', 'helvetica', ...
        'Clipping', 'off');
    text((preBase + staticDur)/params.bin + (movingDur/params.bin)/2, yTop, 'Moving', ...
        'HorizontalAlignment', 'center', 'FontSize', 7, 'FontName', 'helvetica', ...
        'Clipping', 'off');

    % X-axis: hide labels (shared time axis is shown on the PSTH below)
    xlim([0, round(totalStimDur + preBase*2) / params.bin]);
    xticks([0, preBase/params.bin : 600/params.bin : (totalStimDur+preBase*2)/params.bin, ...
        round((totalStimDur+preBase*2)/100)*100 / params.bin]);
    xticklabels([]);

    % Y-axis: ticks at the center of each (TF x SF) block within each angle block
    secondaryN = sfN * tfN;   % Number of conditions per angle
    yt = [0];
    for d = 1:angleN
        block_ticks = 1 : trialDivision*2*secondaryN : (nT/angleN)-1+trialDivision*secondaryN;
        yt = [yt, block_ticks + max(yt) + trialDivision - 1]; %#ok<AGROW>
    end
    yt = yt(2:end-1);  % Remove leading sentinel
    yticks(yt);
    yticklabels(repmat( ...
        trialDivision : trialDivision*2*secondaryN : (nT/angleN)-1+trialDivision*secondaryN, ...
        1, angleN));

    ax = gca;
    ax.YAxis.FontSize = 8;
    ax.YAxis.FontName = 'helvetica';
    ylabel('Trials', 'FontSize', 10, 'FontName', 'helvetica');

    % ==================================================================
    % 7b. Identify best response window (searched across the FULL window,
    %     i.e., both phases together — no phase preference is imposed)
    % ==================================================================

    if params.SelectedWindow

        % Step 1: Mean firing rate per condition group -> find best group
        j = 1;
        meanMr = zeros(1, nT / trialDivision);
        for i = 1:trialDivision:nT
            meanMr(j) = mean(Mr2(i:i+trialDivision-1, :), 'all');
            j = j + 1;
        end
        [~, maxRespIn] = max(meanMr);
        maxRespIn = maxRespIn - 1;  % Convert to 0-based offset for trial range arithmetic

        % Step 2: Extract raster of best condition group [trialDivision x nBins]
        X = Mr2(maxRespIn*trialDivision+1 : maxRespIn*trialDivision+trialDivision, :);

        window = 500;  % Response window width for best-bin search (ms)

        X(X > 1) = 1;  % Clip to [0 1] so outliers don't dominate window selection

        [n_rows, n_cols] = size(X);                           % n_rows = trialDivision
        nWinPos = n_cols - round(window / params.bin) + 1;    % Number of window positions

        % Compute mean firing rate inside every sliding window, for every trial
        % window_means: [trialDivision x nWinPos]
        window_means = zeros(n_rows, nWinPos);
        for col = 1:nWinPos
            window_means(:, col) = mean(X(:, col : col + round(window/params.bin) - 1), 2);
        end

        % Find the (trial, window) pair with the global maximum mean rate
        [~, linear_idx]      = max(window_means(:));
        [best_row, best_col] = ind2sub(size(window_means), linear_idx);

        % Convert window start from bins to ms (relative to trial onset, i.e., after preBase)
        start = best_col * params.bin;   % (ms) offset from trial start (0 = static onset)

        % Which phase did the window land in?
        if start < staticDur
            bestPhase = 'Static';
        else
            bestPhase = 'Moving';
        end

    else
        % Use the pre-computed NeuronVals from whichever phase has a higher raw rate
        [~, bestPhaseIdx] = max([ ...
            max(NeuronResp.Static.NeuronVals(u, :, 4)), ...   % Max raw rate in static phase
            max(NeuronResp.Moving.NeuronVals(u, :, 4))]);     % Max raw rate in moving phase

        phaseNames = ["Static", "Moving"];
        bestPhase  = phaseNames(bestPhaseIdx);

        [~, maxRespIn] = max(NeuronResp.(bestPhase).NeuronVals(u, :, 4));
        % Column 3 = MaxWinBin (bin index); convert to ms
        start     = NeuronResp.(bestPhase).NeuronVals(u, maxRespIn, 3) * NeuronResp.params.binRaster - 20;
        window    = 500;
        maxRespIn = maxRespIn - 1;  % 0-based offset
        best_row  = 1;              % FIX: was uninitialized in original; default = 1
    end

    % Absolute trial indices belonging to the best condition group
    trials = maxRespIn*trialDivision+1 : maxRespIn*trialDivision + trialDivision;

    % Highlight the selected condition group (light grey band, full time extent)
    y1 = maxRespIn*trialDivision + trialDivision + 0.5;
    y2 = maxRespIn*trialDivision + 0.5;
    patch([0, (preBase*2+totalStimDur)/params.bin, (preBase*2+totalStimDur)/params.bin, 0], ...
        [y2, y2, y1, y1], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

    % Absolute index of the single best trial (from sliding window search)
    RasterTrials = trials(best_row);

    % Red patch: mark the best trial and its response window in the raster
    % start is in ms relative to static onset; offset by preBase for display bins
    patch([(preBase + start)/params.bin, (preBase + start + window)/params.bin, ...
        (preBase + start + window)/params.bin, (preBase + start)/params.bin], ...
        [RasterTrials-0.5, RasterTrials-0.5, RasterTrials+0.5, RasterTrials+0.5], ...
        'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % ==================================================================
    % PANEL 3 (rows 17-18): PSTH across the full combined window
    % ==================================================================

    subplot(18, 1, [17 18]);

    % Rebuild 1 ms-resolution spike matrix for all trials (needed for PSTH bins)
    MRhist = BuildBurstMatrix(goodU(:, u), ...
        round(p.t), ...                              % 1 ms resolution
        round(directimesSorted - preBase), ...       % Window start per trial
        round(totalStimDur + preBase*2));            % Full combined window length

    % Select only the best condition group for the PSTH
    MRhist = squeeze(MRhist(trials, :, :));   % [nTrialsInGroup x nTimePoints_ms]

    [nT2, nB2] = size(MRhist);  % nT2 = trials in group, nB2 = window duration ms

    % Collect spike times (ms from trial start)
    spikeTimes = repmat(1:nB2, nT2, 1);       % Index grid matching MRhist
    spikeTimes = spikeTimes(logical(MRhist));  % Keep only spike locations

    % PSTH bin width: use wider bins for longer stimuli to reduce noise
    binWidth = 125;   % Default (ms)
    if nBins > 300
        binWidth = 250;
    end

    edges      = 1:binWidth:round(totalStimDur + preBase*2);  % Bin edges (ms)
    psthCounts = histcounts(spikeTimes, edges);                % Spike count per bin

    % Normalize to firing rate [spk/s]: counts / (bin_s * nTrials)
    psthRate = (psthCounts / (binWidth * nT2)) * 1000;

    b = bar(edges(1:end-1), psthRate, 'histc');
    b.FaceColor       = 'k';
    b.FaceAlpha       = 0.3;
    b.MarkerEdgeColor = 'none';

    xlim([0, round((totalStimDur + preBase*2) / 100) * 100]);

    try   % Guard against all-zero PSTH (std=0 makes ylim fail)
        ylim([0, max(psthRate) + std(psthRate)]);
    catch
        close(fig);
        ur = ur + 1;
        continue
    end

    % X-axis ticks at 600 ms intervals; labels converted from ms to seconds
    xticks([0, preBase:600:(totalStimDur+preBase*2), ...
        round((totalStimDur+preBase*2)/100)*100]);

    % Three xline markers matching the raster above
    xline(preBase,                    'LineWidth', 1.5);          % Stim on
    xline(preBase + staticDur,  '--', 'LineWidth', 1.2);          % Phase transition
    xline(preBase + totalStimDur,     'LineWidth', 1.5);          % Stim off

    xticklabels([-(preBase), 0:600:round((totalStimDur/100))*100, ...
        round((totalStimDur/100))*100 + 2*preBase] ./ 1000);

    ax = gca;
    ax.XAxis.FontSize = 8;  ax.XAxis.FontName = 'helvetica';
    ax.YAxis.FontSize = 8;  ax.YAxis.FontName = 'helvetica';
    ylabel('[spk/s]', 'FontSize', 10, 'FontName', 'helvetica');
    xlabel('Time [s]',  'FontSize', 10, 'FontName', 'helvetica');

    ylims = ylim;
    yticks([round(ylims(2)/10)*5, ceil(ylims(2)/10)*10]);  % Two clean ticks

    % ==================================================================
    % PANEL 1 (rows 1-3): Raw AP/LFP trace for the best single trial
    % ==================================================================

    bin3 = 1;  % 1 ms bins for raw trace extraction

    % Build spike matrix around the response window for all group trials
    % start is in ms from static onset; add preBase offset to get absolute ms
    trialM = BuildBurstMatrix(goodU(:, u), ...
        round(p.t / bin3), ...
        round((directimesSorted + start) / bin3), ...   % Align window to response onset
        round(window / bin3));                           % 500 ms window

    TrialM = squeeze(trialM(trials, :, :))';   % [nBins x nTrialsInGroup]

    % best_row already identifies the highest-response trial from the window search
    chan = goodU(1, u);  % Recording channel for this neuron

    subplot(18, 1, [1 3]);

    % Absolute start time of the raw trace (ms in recording time):
    % align to start of the response window (start ms after static onset, minus preBase)
    startTimes = directimesSorted(RasterTrials) + start - preBase;  % (ms)

    % Binary spike vector for the selected trial in the response window
    spikes = squeeze(BuildBurstMatrix(goodU(:, u), round(p.t), round(startTimes), round(window)));

    % Render raw voltage trace with spike overlay
    [fig, ~, ~] = PlotRawDataNP(obj, fig=fig, chan=chan, ...
        startTimes=startTimes, window=window, spikeTimes=spikes);

    ax = gca;
    ax.YAxis.FontSize = 8;  ax.YAxis.FontName = 'helvetica';

    xlims = xlim;
    xticks(0:(xlims(2)/5):xlims(2));   % 5 evenly spaced ticks
    xticklabels(0:100:window);
    ax.XAxis.FontSize = 8;  ax.XAxis.FontName = 'helvetica';
    ax.XRuler.TickDirection = 'out';
    ax.XAxisLocation = 'bottom';

    % Mark the stimulus event nearest to this window
    % (start is the offset from static onset; if start > staticDur, this is the phase change)
    xline(-start / 1000, 'LineWidth', 1.5);

    xlabel('Time [ms]', 'FontName', 'helvetica', 'FontSize', 10);
    ylabel('[\muV]',    'FontSize', 10,           'FontName', 'helvetica');

    % --- Title: neuron ID + best condition + p-values for both phases --------
    bestAngle = C(maxRespIn*trialDivision + 1, 2);  % Angle of best condition group (deg)
    bestTF    = C(maxRespIn*trialDivision + 1, 3);  % TF of best condition (Hz)
    bestSF    = C(maxRespIn*trialDivision + 1, 4);  % SF of best condition (cyc/deg)

    title({ ...
        sprintf('U.%d  Chan-%d  Phy-%d  |  pS=%.4f  pM=%.4f', ...
            u, channels(ur), phy_IDg(u), pvalsS(u), pvalsM(u)), ...
        sprintf('Best: %.0f deg | TF=%.1f Hz | SF=%.3f c/deg | window in [%s]', ...
            bestAngle, bestTF, bestSF, bestPhase) ...
        });

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

    % Keep the last figure open; close all intermediate ones
    if ur ~= length(eNeuron)
        close(fig);
    end

    ur = ur + 1;  % MUST be reached in every code path above

end  % end neuron loop

end  % end plotRaster