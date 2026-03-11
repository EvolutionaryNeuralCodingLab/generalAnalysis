function plotRawWaveforms(obj, unitID, params)
% plotRawWaveforms - Plot raw spike waveforms from KS4 output, Phy-style
%                    Channels are drawn at their true probe x/y positions.
%                    Optionally plots an auto-correlogram in a separate figure.
%
% INPUTS:
%   obj    - Visual stimulation object 
%   unitID - cluster ID to plot (single unit)
%
% OPTIONAL NAME-VALUE PARAMS:
%   nWaveforms    - number of random waveforms to plot      (default: 100)
%   nChanAround   - channels above/below max amp channel    (default: 4)
%   nPre          - samples before spike peak               (default: 20)
%   nPost         - samples after spike peak                (default: 61)
%   showCorr      - plot auto-correlogram                   (default: false)
%   corrWin       - correlogram half-window in ms           (default: 100)
%   corrBin       - correlogram bin size in ms              (default: 1)
%
% EXAMPLES:
%   plotRawWaveforms(obj, 42)
%   plotRawWaveforms(obj, 42, nWaveforms=200, nChanAround=6)
%   plotRawWaveforms(obj, 42, showCorr=true, corrWin=50, corrBin=0.5)

arguments (Input)
    obj
    unitID                (1,1) double
    params.nWaveforms     (1,1) double  = 100
    params.nChanAround    (1,1) double  = 10
    params.nPre           (1,1) double  = 20
    params.nPost          (1,1) double  = 61
    params.showCorr       (1,1) logical = false
    params.corrWin        (1,1) double  = 100
    params.corrBin        (1,1) double  = 1
end

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

%% Load KS4 output
spike_times    = readNPY(fullfile(ksDir, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
templates      = readNPY(fullfile(ksDir, 'templates.npy'));         % [nUnits x T x nCh]
chan_map        = readNPY(fullfile(ksDir, 'channel_map.npy'));       % [nCh x 1], 0-indexed
chan_pos        = readNPY(fullfile(ksDir, 'channel_positions.npy')); % [nCh x 2]

%% Find template index for this unit
unit_ids = (0 : size(templates, 1) - 1)';
tmpl_idx = find(unit_ids == unitID);
if isempty(tmpl_idx)
    error('Unit %d not found in templates.npy', unitID);
end

%% Find best channel (max peak-to-peak across template channels)
unit_template = squeeze(templates(tmpl_idx, :, :));  % [T x nCh]
p2p           = max(unit_template) - min(unit_template);
[~, best_tmpl_chan] = max(p2p);

% Get probe positions for all template channels via chan_map
% chan_pos is [nTemplateCh x 2]: col1=x, col2=y (from KS4, in µm)
% Find nChanAround closest channels to best channel by Euclidean distance
best_xy  = chan_pos(best_tmpl_chan, :);                          % [1 x 2]
dists    = sqrt(sum((chan_pos - best_xy).^2, 2));                % [nTemplateCh x 1]
[~, sorted_idx] = sort(dists, 'ascend');
chan_indices = sorted_idx(1 : min(params.nChanAround + 1, numel(dists)))';
n_chans_plot = numel(chan_indices);

% Index of best channel within the plotted subset
best_local_idx = find(chan_indices == best_tmpl_chan);

% Map to binary file channels (1-indexed for MATLAB)
bin_chans    = chan_map(chan_indices) + 1;   % [n_chans_plot x 1], 1-indexed
best_bin_chan = bin_chans(best_local_idx);

%% Get spike times for this unit
st = double(spike_times(spike_clusters == unitID));
if numel(st) < 2, error('Unit %d has fewer than 2 spikes.', unitID); end
fprintf('Unit %d: %d total spikes, extracting %d waveforms\n', ...
    unitID, numel(st), min(params.nWaveforms, numel(st)));

idx    = randperm(numel(st), min(params.nWaveforms, numel(st)));
st_sub = st(idx);

%% Extract waveforms from binary
waveform_len = params.nPre + params.nPost + 1;
finfo        = dir(binPath);
n_samp_total = finfo.bytes / (n_channels * 2);
fid          = fopen(binPath, 'rb');

waveforms = NaN(n_chans_plot, waveform_len, numel(st_sub));

for si = 1:numel(st_sub)
    s0 = st_sub(si) - params.nPre;
    s1 = st_sub(si) + params.nPost;
    if s0 < 1 || s1 > n_samp_total, continue; end

    fseek(fid, (s0 - 1) * n_channels * 2, 'bof');
    raw = fread(fid, [n_channels, waveform_len], '*int16');
    if size(raw, 2) < waveform_len, continue; end

    waveforms(:, :, si) = double(raw(bin_chans, :)) * uV_per_bit;
end
fclose(fid);

% Baseline subtract (mean of pre-spike window)
baseline  = mean(waveforms(:, 1:params.nPre, :), 2, 'omitnan');
waveforms = waveforms - baseline;

%% Compute correlogram if requested
if params.showCorr
    [ccg_counts, ccg_bins] = computeACG(st, sample_rate, params.corrWin, params.corrBin);
end

%% ---- Spatial positions for plotted channels ----
% chPos is [2 x nAllCh]: row 1 = x (shank col), row 2 = y (depth)
ch_x = chPos(1, bin_chans);   % [1 x n_chans_plot]
ch_y = chPos(2, bin_chans);   % [1 x n_chans_plot]

% Detect inter-channel pitch from all channels on the probe
x_unique  = unique(chPos(1,:));
y_unique  = unique(chPos(2,:));
x_spacing = min(diff(sort(x_unique)));
y_spacing = min(diff(sort(y_unique)));

if isempty(x_spacing) || numel(x_unique) == 1
    x_spacing = 32;   % fallback NP1 column pitch
end
if isempty(y_spacing) || numel(y_unique) == 1
    y_spacing = 20;   % fallback NP1 row pitch
end

% Time axis scaled to fit in x_spacing (80% fill)
t_ms    = (-params.nPre : params.nPost) / sample_rate * 1000;
t_scale = 0.8 * x_spacing / (t_ms(end) - t_ms(1));   % µm per ms

% Amplitude scale: normalise so max p2p fits in y_spacing (80% fill)
mean_wf = mean(waveforms, 3, 'omitnan');
std_wf  = std(waveforms,  0, 3, 'omitnan');
max_p2p = max(max(mean_wf, [], 2) - min(mean_wf, [], 2));
if max_p2p == 0, max_p2p = 1; end
amp_scale = 0.8 * y_spacing / max_p2p;   % µm per µV

%% ---- Colours: best channel = red, all others = blue ----
col_default = [0.25 0.45 0.75];   % blue
col_best    = [0.85 0.20 0.15];   % red

%% ---- Waveform figure ----
figure('Color', 'w', 'Name', sprintf('Unit %d — Waveforms', unitID));
ax = axes('Color', 'w');
hold(ax, 'on');

for ci = 1:n_chans_plot
    cx = ch_x(ci);
    cy = ch_y(ci);

    if ci == best_local_idx
        col = col_best;
    else
        col = col_default;
    end

    x_wf = cx + t_ms * t_scale;

    % Individual waveforms (translucent)
    wf_ci = squeeze(waveforms(ci, :, :));   % [nSamples x nWaveforms]
    y_wf  = cy + wf_ci * amp_scale;
    plot(ax, x_wf, y_wf, 'Color', [col, 0.12], 'LineWidth', 0.5);

    % Std shading
    upper = cy + (mean_wf(ci,:) + std_wf(ci,:)) * amp_scale;
    lower = cy + (mean_wf(ci,:) - std_wf(ci,:)) * amp_scale;
    fill(ax, [x_wf, fliplr(x_wf)], [upper, fliplr(lower)], ...
        col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Mean waveform
    y_mean = cy + mean_wf(ci,:) * amp_scale;
    plot(ax, x_wf, y_mean, 'k', 'LineWidth', 2);

    % Channel label: two rows, just left of waveform start
    text(ax, x_wf(1) - 2, cy + amp_scale * 0, ...
        sprintf('ch%d\n(%g, %g)', bin_chans(ci), cx, cy), ...
        'FontSize', 7, 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', 'Color', col);
end

%% ---- L-shaped scale bar ----
% Fixed scale: 2 ms horizontal, 200 µV vertical
sb_ms  = 1;      % ms
sb_uv  = 200;    % µV
sb_xlen = sb_ms  * t_scale;    % µm
sb_ylen = sb_uv  * amp_scale;  % µm

% Position: to the right of the bottom-right waveform, at the same y level
[~, bottom_right_ci] = min(ch_y - ch_x * 1e-6);  % lowest y, rightmost x as tiebreak
br_cx = ch_x(bottom_right_ci);
br_cy = ch_y(bottom_right_ci);

sb_gap = 0.2 * x_spacing;                        % horizontal gap from last waveform
sb_ox  = br_cx + t_ms(end) * t_scale + sb_gap;   % L corner x: just right of waveform end
sb_oy  = br_cy;                                   % L corner y: same level as that channel

% Draw L: vertical arm then horizontal arm, meeting at bottom-left corner
plot(ax, [sb_ox, sb_ox],          [sb_oy, sb_oy - sb_ylen], 'k', 'LineWidth', 2);
plot(ax, [sb_ox, sb_ox + sb_xlen],[sb_oy, sb_oy],           'k', 'LineWidth', 2);

% Labels
text(ax, sb_ox - 2, sb_oy - sb_ylen/2, sprintf('%d µV', sb_uv), ...
    'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation',90);
text(ax, sb_ox + sb_xlen/2, sb_oy + 2, sprintf('%d ms', sb_ms), ...
    'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

%% Axis cosmetics — no tick marks, no box
title(ax, sprintf('Unit %d  |  %d waveforms  |  best ch: %d', ...
    unitID, numel(st_sub), best_bin_chan), 'FontSize', 12);
set(ax, 'XTick', [], 'YTick', [], 'YDir', 'normal');
axis(ax, 'tight');
box(ax, 'off');
axis(ax, 'off');   % hide axes entirely — scale bar carries all metric info

%% ---- ACG figure (separate) ----
if params.showCorr
    figure('Color', 'w', 'Name', sprintf('Unit %d — ACG', unitID));
    ax_corr = axes;

    bar(ax_corr, ccg_bins, ccg_counts, 1, ...
        'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none');
    hold(ax_corr, 'on');
    xline(ax_corr, 0, '--k', 'Alpha', 0.4);

    % Shade refractory period (±2 ms)
    ylims = ylim(ax_corr);
    patch(ax_corr, [-2 2 2 -2], [0 0 ylims(2) ylims(2)], ...
        'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

    xlabel(ax_corr, 'Lag (ms)');
    ylabel(ax_corr, 'Spike count');
    title(ax_corr, sprintf('Unit %d ACG | RP 2 ms | bin %.1f ms | win ±%d ms', ...
        unitID, params.corrBin, params.corrWin), 'FontSize', 12);
    xlim(ax_corr, [-params.corrWin params.corrWin]);
    box(ax_corr, 'off');
end

end % main function


%% =========================================================================
function [counts, bin_centers] = computeACG(spike_times_samples, fs, win_ms, bin_ms)
% Compute auto-correlogram for a single unit
%   spike_times_samples - spike times in samples
%   fs                  - sampling rate (Hz)
%   win_ms              - half-window in ms
%   bin_ms              - bin size in ms

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