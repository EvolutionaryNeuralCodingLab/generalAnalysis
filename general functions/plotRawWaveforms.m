function [fig1, fig2] = plotRawWaveforms(obj, unitIDs, params)
% plotRawWaveforms - Plot raw spike waveforms from KS4 output, Phy-style
%                    Each unit is shown in its own tile at true probe positions.
%                    Optionally plots ACGs for all units in a single tiled figure.
%
% INPUTS:
%   obj     - Visual stimulation object
%   unitIDs - scalar or vector of cluster IDs to plot  e.g. 42  or  [3 7 42]
%
% OPTIONAL NAME-VALUE PARAMS:
%   nWaveforms    - number of random waveforms to plot      (default: 100)
%   nChanAround   - nearest channels around max amp channel (default: 10)
%   nPre          - samples before spike peak               (default: 20)
%   nPost         - samples after spike peak                (default: 61)
%   showCorr      - plot auto-correlogram figure            (default: false)
%   corrWin       - correlogram half-window in ms           (default: 100)
%   corrBin       - correlogram bin size in ms              (default: 1)
%
% EXAMPLES:
%   plotRawWaveforms(obj, 42)
%   plotRawWaveforms(obj, [3 7 42], nWaveforms=200, nChanAround=6)
%   plotRawWaveforms(obj, [3 7 42], showCorr=true, corrWin=50, corrBin=0.5)

arguments (Input)
    obj
    unitIDs               (1,:) double
    params.nWaveforms     (1,1) double  = 100
    params.nChanAround    (1,1) double  = 8
    params.nPre           (1,1) double  = 20
    params.nPost          (1,1) double  = 61
    params.showCorr       (1,1) logical = false
    params.corrWin        (1,1) double  = 100
    params.corrBin        (1,1) double  = 1
end

nUnits = numel(unitIDs);

%% Paths
ksDir        = obj.spikeSortingFolder;
recordingDir = obj.dataObj.recordingDir;

%% Settings from obj
n_channels  = str2double(obj.dataObj.nSavedChansImec);
sample_rate = obj.dataObj.samplingFrequency;
uV_per_bit  = unique(obj.dataObj.MicrovoltsPerAD);
chPos       = obj.dataObj.chLayoutPositions;   % [2 x nAllCh]: row1=x, row2=y

fprintf('Settings — nCh: %d | Fs: %d Hz | uV/bit: %.4f\n', ...
    n_channels, sample_rate, uV_per_bit);

%% Find binary file
binFiles = dir(fullfile(recordingDir, '*.bin'));
if isempty(binFiles), binFiles = dir(fullfile(recordingDir, '*.dat')); end
if isempty(binFiles), error('No .bin or .dat file found in: %s', recordingDir); end
binPath = fullfile(recordingDir, binFiles(1).name);
fprintf('Using binary file: %s\n', binPath);

%% Load KS4 output (once, shared across all units)
spike_times    = readNPY(fullfile(ksDir, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
templates      = readNPY(fullfile(ksDir, 'templates.npy'));         % [nUnits x T x nCh]
chan_map        = readNPY(fullfile(ksDir, 'channel_map.npy'));       % [nCh x 1], 0-indexed
chan_pos        = readNPY(fullfile(ksDir, 'channel_positions.npy')); % [nCh x 2]

unit_ids_ks = (0 : size(templates, 1) - 1)';

%% Probe pitch (shared across all units)
x_unique  = unique(chPos(1,:));
y_unique  = unique(chPos(2,:));
x_spacing = min(diff(sort(x_unique)));
y_spacing = min(diff(sort(y_unique)));
if isempty(x_spacing) || numel(x_unique) == 1, x_spacing = 32; end
if isempty(y_spacing) || numel(y_unique) == 1, y_spacing = 20; end

t_ms = (-params.nPre : params.nPost) / sample_rate * 1000;

%% Colours
col_default = [0.25 0.45 0.75];   % blue
col_best    = [0.85 0.20 0.15];   % red

%% ---- Extract data for each unit ----
finfo        = dir(binPath);
n_samp_total = finfo.bytes / (n_channels * 2);
fid          = fopen(binPath, 'rb');

unitData = struct();   % will hold per-unit results

for ui = 1:nUnits
    unitID = unitIDs(ui);

    % Template index
    tmpl_idx = find(unit_ids_ks == unitID);
    if isempty(tmpl_idx)
        warning('Unit %d not found in templates.npy, skipping.', unitID);
        unitData(ui).valid = false;
        continue
    end

    % Best channel by p2p on template
    unit_template = squeeze(templates(tmpl_idx, :, :));   % [T x nCh]
    p2p           = max(unit_template) - min(unit_template);
    [~, best_tmpl_chan] = max(p2p);

    % nChanAround nearest channels by Euclidean distance on probe
    best_xy  = chan_pos(best_tmpl_chan, :);
    dists    = sqrt(sum((chan_pos - best_xy).^2, 2));
    [~, sorted_idx]  = sort(dists, 'ascend');
    chan_indices      = sorted_idx(1 : min(params.nChanAround + 1, numel(dists)))';
    n_chans_plot      = numel(chan_indices);
    best_local_idx    = find(chan_indices == best_tmpl_chan);

    bin_chans    = chan_map(chan_indices) + 1;   % 1-indexed
    best_bin_chan = bin_chans(best_local_idx);

    % Spike times for this unit
    st = double(spike_times(spike_clusters == unitID));
    if numel(st) < 2
        warning('Unit %d has fewer than 2 spikes, skipping.', unitID);
        unitData(ui).valid = false;
        continue
    end

    % Random subsample
    idx    = randperm(numel(st), min(params.nWaveforms, numel(st)));
    st_sub = st(idx);
    fprintf('Unit %d: %d total spikes, extracting %d waveforms\n', ...
        unitID, numel(st), numel(st_sub));

    % Extract waveforms
    waveform_len = params.nPre + params.nPost + 1;
    waveforms    = NaN(n_chans_plot, waveform_len, numel(st_sub));

    for si = 1:numel(st_sub)
        s0 = st_sub(si) - params.nPre;
        s1 = st_sub(si) + params.nPost;
        if s0 < 1 || s1 > n_samp_total, continue; end
        fseek(fid, (s0 - 1) * n_channels * 2, 'bof');
        raw = fread(fid, [n_channels, waveform_len], '*int16');
        if size(raw, 2) < waveform_len, continue; end
        waveforms(:, :, si) = double(raw(bin_chans, :)) * uV_per_bit;
    end

    % Baseline subtract
    baseline  = mean(waveforms(:, 1:params.nPre, :), 2, 'omitnan');
    waveforms = waveforms - baseline;

    % Store
    unitData(ui).valid         = true;
    unitData(ui).unitID        = unitID;
    unitData(ui).waveforms     = waveforms;
    % Exclude outlier waveforms based on peak-to-peak MAD
    % Compute p2p amplitude for each waveform (max across channels and time)
    wf_p2p   = squeeze(max(max(waveforms,[],1),[],2) - ...
        min(min(waveforms,[],1),[],2));  % [1 x nWaveforms]
    wf_median = median(wf_p2p, 'omitnan');
    wf_mad    = median(abs(wf_p2p - wf_median), 'omitnan');
    inlier_mask = abs(wf_p2p - wf_median) < 5 * wf_mad;  % 5-MAD threshold
    fprintf('Unit %d: %d/%d waveforms kept for envelope (outlier rejection)\n', ...
        unitID, sum(inlier_mask), numel(inlier_mask));

    unitData(ui).mean_wf  = mean(waveforms(:,:,inlier_mask), 3, 'omitnan');
    unitData(ui).std_wf   = std(waveforms(:,:,inlier_mask),  0, 3, 'omitnan');
    unitData(ui).bin_chans     = bin_chans;
    unitData(ui).best_bin_chan = best_bin_chan;
    unitData(ui).best_local_idx= best_local_idx;
    unitData(ui).n_chans_plot  = n_chans_plot;
    unitData(ui).ch_x          = chPos(1, bin_chans);
    unitData(ui).ch_y          = chPos(2, bin_chans);
    unitData(ui).st            = st;
    unitData(ui).n_wf          = numel(st_sub);

    % ACG
    if params.showCorr
        [unitData(ui).ccg_counts, unitData(ui).ccg_bins] = ...
            computeACG(st, sample_rate, params.corrWin, params.corrBin);
    end
end
fclose(fid);

%% ---- Waveform figure: one tile per unit ----
% Determine tiled layout dimensions
nCols = min(nUnits, 4);
nRows = ceil(nUnits / nCols);

fig1 = figure('Color', 'w', 'Name', 'Waveforms');
wf_tl = tiledlayout(fig1, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
title(wf_tl, 'Raw Waveforms', 'FontSize', 13, 'FontWeight', 'bold');

for ui = 1:nUnits
    if ~unitData(ui).valid, continue; end

    d            = unitData(ui);
    mean_wf      = d.mean_wf;
    std_wf       = d.std_wf;
    ch_x         = d.ch_x;
    ch_y         = d.ch_y;
    bin_chans    = d.bin_chans;
    best_local_idx = d.best_local_idx;
    n_chans_plot = d.n_chans_plot;

    % Per-unit amplitude scale: use mean±std envelope to prevent overlap
    % on noisy units (large std compresses the scale automatically)
    upper_env = max(mean_wf + std_wf, [], 2);   % [nCh x 1]
    lower_env = min(mean_wf - std_wf, [], 2);
    max_extent = max(upper_env - lower_env);
    if max_extent == 0, max_extent = 1; end
    amp_scale = 0.8 * y_spacing / max_extent;
    t_scale   = 0.8 * x_spacing / (t_ms(end) - t_ms(1));

    % Scale bar µV: round max amplitude to nearest 50 µV
    sb_uv = max(50, round(max_extent / 50) * 50);

    ax = nexttile(wf_tl);
    hold(ax, 'on');

    for ci = 1:n_chans_plot
        cx = ch_x(ci);
        cy = ch_y(ci);
        col = col_default;
        if ci == best_local_idx, col = col_best; end

        x_wf = cx + t_ms * t_scale;

        % Individual waveforms
        wf_ci = squeeze(d.waveforms(ci, :, :));
        plot(ax, x_wf, cy + wf_ci * amp_scale, ...
            'Color', [col, 0.12], 'LineWidth', 0.5);

        % Std shading
        upper = cy + (mean_wf(ci,:) + std_wf(ci,:)) * amp_scale;
        lower = cy + (mean_wf(ci,:) - std_wf(ci,:)) * amp_scale;
        fill(ax, [x_wf, fliplr(x_wf)], [upper, fliplr(lower)], ...
            col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

       % Mean waveform (black), with coloured std shading
        plot(ax, x_wf, cy + mean_wf(ci,:) * amp_scale, ...
            'Color', 'k', 'LineWidth', 2);

        % Channel label (two rows, left of waveform start)
        text(ax, x_wf(1) - 2, cy, ...
            sprintf('ch%d\n(%g,%g)', bin_chans(ci), cx, cy), ...
            'FontSize', 6, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', 'Color', col);
    end

    % L-scale bar: bottom-right channel of this unit
    sb_ms   = 1;    % sb_uv already set above
    sb_xlen = sb_ms * t_scale;
    sb_ylen = sb_uv * amp_scale;

    [~, br_ci] = min(ch_y - ch_x * 1e-6);
    sb_ox = ch_x(br_ci) + t_ms(end) * t_scale + 0.2 * x_spacing;
    sb_oy = ch_y(br_ci);

    plot(ax, [sb_ox, sb_ox],           [sb_oy, sb_oy - sb_ylen], 'k', 'LineWidth', 2);
    plot(ax, [sb_ox, sb_ox + sb_xlen], [sb_oy, sb_oy],           'k', 'LineWidth', 2);
    text(ax, sb_ox - 2, sb_oy - sb_ylen/2, sprintf('%d µV', sb_uv), ...
        'FontSize', 7, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Rotation', 90);
    text(ax, sb_ox + sb_xlen/2, sb_oy + 2, sprintf('%d ms', sb_ms), ...
        'FontSize', 7, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

    title(ax, sprintf('Unit %d  |  ch%d  |  n=%d', ...
        d.unitID, d.best_bin_chan, d.n_wf), 'FontSize', 9);
    axis(ax, 'tight');
    axis(ax, 'off');
end

%% ---- ACG figure: one tile per unit ----
if params.showCorr
    fig2 = figure('Color', 'w', 'Name', 'ACGs');
    acg_tl = tiledlayout(fig2, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(acg_tl, sprintf('ACG | RP 2 ms | bin %.1f ms | win ±%d ms', ...
    params.corrBin, params.corrWin), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel(acg_tl, 'Lag (ms)');
    ylabel(acg_tl, 'Spike count');

    for ui = 1:nUnits
        if ~unitData(ui).valid, continue; end
        d = unitData(ui);

        ax_c = nexttile(acg_tl);
        bar(ax_c, d.ccg_bins, d.ccg_counts, 1, ...
            'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none');
        hold(ax_c, 'on');
        xline(ax_c, 0, '--k', 'Alpha', 0.4);

        ylims = ylim(ax_c);
        patch(ax_c, [-2 2 2 -2], [0 0 ylims(2) ylims(2)], ...
            'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

        xlim(ax_c, [-params.corrWin params.corrWin]);
        title(ax_c, sprintf('Unit %d', d.unitID), 'FontSize', 9);
        box(ax_c, 'off');
    end
else
    fig2 = [];
end

end % main function


%% =========================================================================
function [counts, bin_centers] = computeACG(spike_times_samples, fs, win_ms, bin_ms)
st_ms       = spike_times_samples / fs * 1000;
edges       = -win_ms : bin_ms : win_ms;
bin_centers = edges(1:end-1) + bin_ms / 2;
counts      = zeros(1, numel(bin_centers));
for i = 1:numel(st_ms)
    diffs    = st_ms - st_ms(i);
    diffs(i) = NaN;
    diffs    = diffs(diffs > -win_ms & diffs < win_ms);
    counts   = counts + histcounts(diffs, edges);
end
end