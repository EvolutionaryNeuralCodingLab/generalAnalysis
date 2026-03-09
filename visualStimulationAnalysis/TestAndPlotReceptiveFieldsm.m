%% =========================================================================
%  EXAMPLE: Receptive Field Calculation + Spatial Tuning Index
%  Assumes 'obj' is already loaded (your recording object)
% =========================================================================



%% -----------------------------------------------------------------------
% 0. Setup
% -----------------------------------------------------------------------
clear; clc;

% Load your recording object however you normally do it
%
ex = [85]; %97
NP = loadNPclassFromTable(ex); %73 81
vs = linearlyMovingBallAnalysis(NP,Session=1);

%% -----------------------------------------------------------------------
% 1. Calculate Receptive Fields
% -----------------------------------------------------------------------
fprintf('=== STEP 1: Calculating Receptive Fields ===\n')

RFresults = vs.CalculateReceptiveFieldsV2( ...
    'overwrite',            true,      ...  % use saved file if it exists
    'speed',                1,          ...  % select speed index
    'delay',                250,        ...  % ms delay for RF snapshot
    'reduceFactor',         20,         ...  % screen downsampling factor
    'nShuffle',             100,        ...  % number of shuffle iterations
    'shuffleMode',          'within',   ...  % shuffle within stimulus category
    'AllResponsiveNeurons', true,       ...  % only use responsive neurons
    'testConvolution',      false);          % set true to debug with one neuron

% --- Unpack outputs ---
% Real RF:     [nDir x nSize x nLum x Y x Y x nNeurons]  smoothed
RFreal    = RFresults.RFuDirSizeLumFilt;

% Shuffled RF: [Y x Y x nNeurons x nShuffle]              smoothed
RFshuff   = RFresults.RFuShuffSTFilt;

% Unsmoothed condition RF (for inspection):
% [nDir x nSize x nLum x Y x Y x nNeurons]
RFraw     = RFresults.RFuSTDirSizeLum;

fprintf('RF matrix size (real):    %s\n', mat2str(size(RFreal)))
fprintf('RF matrix size (shuffled):%s\n', mat2str(size(RFshuff)))

%% -----------------------------------------------------------------------
% 2. Quick sanity check: plot RF of one neuron, one condition
% -----------------------------------------------------------------------
exNeuron = 1;    % example neuron index
exDir    = 1;    % direction index
exSize   = 1;    % size index
exLum    = 2;    % luminosity index (e.g. bright ball)

figure('Name', 'Sanity Check: Single Neuron RF');
tiledlayout(1, 3, 'TileSpacing', 'compact');

nexttile
imagesc(squeeze(RFreal(exDir, exSize, exLum, :, :, exNeuron)));
colorbar; axis image; colormap hot;
title(sprintf('Real RF — Neuron %d', exNeuron));
xlabel('X (pixels)'); ylabel('Y (pixels)');

nexttile
imagesc(mean(squeeze(RFshuff(:,:,exNeuron,:)), 3));  % mean across shuffles
colorbar; axis image; colormap hot;
title('Mean Shuffled RF');
xlabel('X (pixels)'); ylabel('Y (pixels)');

nexttile
% Difference: real minus mean shuffle
imagesc(squeeze(RFreal(exDir, exSize, exLum, :, :, exNeuron)) - ...
        mean(squeeze(RFshuff(:,:,exNeuron,:)), 3));
colorbar; axis image; colormap(bluewhitered);  % or any diverging colormap
title('Real minus Mean Shuffle');
xlabel('X (pixels)'); ylabel('Y (pixels)');

%% -----------------------------------------------------------------------
% 3. Calculate Spatial Tuning Index
% -----------------------------------------------------------------------
fprintf('\n=== STEP 2: Calculating Spatial Tuning Index ===\n')

STIresults = vs.CalculateSpatialTuningIndex( ...
    'OneDirection',  'all',   ...   % compute for all directions
    'OneLuminosity', 'all',   ...   % compute for all luminosities
    'OneSize',       'all',   ...   % compute for all sizes
    'plotResult',    false);        % we'll do custom plots below

% STIresults: [nDir x nSize x nLum x nNeurons]  Z-score
fprintf('STI matrix size: %s\n', mat2str(size(STIresults)))

%% -----------------------------------------------------------------------
% 4. Inspect STI distributions
% -----------------------------------------------------------------------

% Load condition labels (to label axes)
responses = vs.ResponseWindow;
fieldName = 'Speed1';
uDir  = unique(responses.(fieldName).C(:,2));
uSize = unique(responses.(fieldName).C(:,4));
uLum  = unique(responses.(fieldName).C(:,6));

nDir  = numel(uDir);
nSize = numel(uSize);
nLum  = numel(uLum);

%% --- 4a. Population histogram of STI (collapsed across all conditions) ---
figure('Name', 'STI Population Distribution');
STI_flat = STIresults(:);
STI_flat = STI_flat(~isnan(STI_flat));

histogram(STI_flat, 40, 'Normalization', 'probability', ...
    'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none');
hold on
xline(0,    'k--', 'LineWidth', 1.5, 'Label', 'Chance');
xline(median(STI_flat), 'r-', 'LineWidth', 2, ...
      'Label', sprintf('Median = %.2f', median(STI_flat)));
xline(1.96, 'g--', 'LineWidth', 1.5, 'Label', 'p<0.05');
xline(3.09, 'm--', 'LineWidth', 1.5, 'Label', 'p<0.001');
xlabel('Spatial Tuning Index (Z-score)');
ylabel('Proportion of neurons');
title('Population STI — all conditions');
box off

fprintf('\nPopulation STI summary:\n')
fprintf('  Median Z-score:         %.3f\n', median(STI_flat))
fprintf('  %% neurons Z > 1.96:    %.1f%%\n', 100*mean(STI_flat > 1.96))
fprintf('  %% neurons Z > 3.09:    %.1f%%\n', 100*mean(STI_flat > 3.09))

%% --- 4b. STI as a function of ball SIZE (collapsed across dir and lum) ---
figure('Name', 'STI vs Ball Size');
STI_bySize = squeeze(mean(mean(STIresults, 1, 'omitnan'), 3, 'omitnan'));
% STI_bySize: [nSize x nNeurons]

% Boxplot across neurons for each size
boxplot(STI_bySize', uSize, 'Labels', arrayfun(@num2str, uSize, 'UniformOutput', false));
hold on
yline(0,    'k--');
yline(1.96, 'g--', 'p<0.05');
xlabel('Ball radius (pixels)');
ylabel('STI (Z-score)');
title('STI vs Ball Size');
box off

%% --- 4c. STI as a function of DIRECTION (collapsed across size and lum) ---
figure('Name', 'STI vs Direction');
dirLabels = {'Right (0)', 'Up (π/2)', 'Left (π)', 'Down (3π/2)'};  % adjust to your convention
STI_byDir = squeeze(mean(mean(STIresults, 2, 'omitnan'), 3, 'omitnan'));
% STI_byDir: [nDir x nNeurons]

boxplot(STI_byDir', 'Labels', dirLabels(1:nDir));
hold on
yline(0,    'k--');
yline(1.96, 'g--', 'p<0.05');
xlabel('Direction of motion');
ylabel('STI (Z-score)');
title('STI vs Direction');
box off

%% --- 4d. Heatmap: mean STI across [Direction x Size] (collapsed over lum) ---
figure('Name', 'STI Heatmap: Direction x Size');
STI_DirSize = squeeze(mean(mean(STIresults, 3, 'omitnan'), 4, 'omitnan'));
% STI_DirSize: [nDir x nSize]

imagesc(uSize, 1:nDir, STI_DirSize);
colorbar;
clim([0 max(STI_DirSize(:))]);
colormap hot;
yticks(1:nDir);
yticklabels(dirLabels(1:nDir));
xlabel('Ball radius (pixels)');
ylabel('Direction');
title('Mean STI — Direction x Size');

%% -----------------------------------------------------------------------
% 5. Identify significantly tuned neurons
% -----------------------------------------------------------------------
% Threshold: Z > 1.96 (p < 0.05, one-tailed) in at least one condition

zThresh = 1.96;

% Best condition per neuron (max STI across all conditions)
STI_maxPerNeuron = squeeze(max(STIresults, [], [1 2 3], 'omitnan'));
% STI_maxPerNeuron: [nNeurons x 1]

isTuned = STI_maxPerNeuron > zThresh;

fprintf('\nNeurons with STI > %.2f in at least one condition: %d / %d (%.1f%%)\n', ...
    zThresh, sum(isTuned), numel(isTuned), 100*mean(isTuned))

%% --- 5a. Compare tuned vs untuned neurons ---
figure('Name', 'Tuned vs Untuned');
hold on
histogram(STI_maxPerNeuron(~isTuned), 30, 'Normalization', 'probability', ...
    'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'DisplayName', 'Not tuned');
histogram(STI_maxPerNeuron( isTuned), 30, 'Normalization', 'probability', ...
    'FaceColor', [0.2 0.6 0.3], 'EdgeColor', 'none', 'DisplayName', 'Tuned');
xline(zThresh, 'r--', 'LineWidth', 2, 'Label', 'Threshold');
legend;
xlabel('Max STI across conditions (Z-score)');
ylabel('Proportion');
title('Tuned vs Untuned Neurons');
box off

%% -----------------------------------------------------------------------
% 6. Deep dive: best condition per tuned neuron
% -----------------------------------------------------------------------
% Find which direction/size/lum gives the highest STI for each neuron

[~, bestCondLinIdx] = max(reshape(STIresults, [], size(STIresults,4)), [], 1);
[bestDir, bestSize, bestLum] = ind2sub([nDir nSize nLum], bestCondLinIdx);

figure('Name', 'Best Condition Distribution');
tiledlayout(1, 3, 'TileSpacing', 'compact');

nexttile
histogram(bestDir(isTuned), 1:nDir+1, 'Normalization', 'probability', ...
    'FaceColor', [0.2 0.5 0.8]);
xticks(1.5:nDir+0.5);
xticklabels(dirLabels(1:nDir));
xlabel('Best direction');
ylabel('Proportion of tuned neurons');
title('Best Direction');

nexttile
histogram(uSize(bestSize(isTuned)), numel(uSize), 'Normalization', 'probability', ...
    'FaceColor', [0.8 0.4 0.2]);
xlabel('Best size (px)');
title('Best Size');

nexttile
histogram(uLum(bestLum(isTuned)), numel(uLum), 'Normalization', 'probability', ...
    'FaceColor', [0.5 0.2 0.8]);
xlabel('Best luminosity');
title('Best Luminosity');

%% -----------------------------------------------------------------------
% 7. Plot RF gallery for tuned neurons
% -----------------------------------------------------------------------
% Show real RF at best condition for each tuned neuron

tunedIdx = find(isTuned);
nToPlot  = min(12, numel(tunedIdx));   % show up to 12 neurons

figure('Name', 'RF Gallery — Tuned Neurons');
tiledlayout(3, 4, 'TileSpacing', 'compact');

for k = 1:nToPlot
    ui  = tunedIdx(k);
    bd  = bestDir(ui);
    bs  = bestSize(ui);
    bl  = bestLum(ui);

    nexttile
    imagesc(squeeze(RFreal(bd, bs, bl, :, :, ui)));
    axis image; colormap hot; colorbar;
    title(sprintf('N%d | Z=%.1f\nDir=%d Sz=%dpx L=%d', ...
        ui, STIresults(bd,bs,bl,ui), ...
        uDir(bd), uSize(bs), uLum(bl)), 'FontSize', 7);
    axis off
end
sgtitle('Receptive Fields at Best Condition (Tuned Neurons)')

%% -----------------------------------------------------------------------
% 8. Compare shuffle modes (optional diagnostic)
% -----------------------------------------------------------------------
% Run STI with each shuffle mode and compare distributions.
% Useful to verify that 'within' gives a more conservative baseline
% than 'global' for condition-selective neurons.

fprintf('\n=== STEP 3 (optional): Shuffle mode comparison ===\n')

runComparison = false;   % set to true to run

if runComparison

    modes     = {'within', 'global', 'within_bins_only'};
    colors    = {[0.2 0.5 0.8], [0.8 0.3 0.2], [0.2 0.7 0.4]};
    STI_modes = cell(1,3);

    for m = 1:3
        fprintf('Running shuffle mode: %s\n', modes{m})
        RF_m = obj.CalculateReceptiveFields( ...
            'overwrite',   true,      ...
            'nShuffle',    100,       ...
            'shuffleMode', modes{m}, ...
            'AllResponsiveNeurons', true);

        STI_modes{m} = obj.CalculateSpatialTuningIndex( ...
            'OneDirection',  'all', ...
            'OneLuminosity', 'all', ...
            'OneSize',       'all');
    end

    figure('Name', 'Shuffle Mode Comparison');
    hold on
    for m = 1:3
        flat = STI_modes{m}(:);
        flat = flat(~isnan(flat));
        histogram(flat, 40, 'Normalization', 'probability', ...
            'FaceColor', colors{m}, 'EdgeColor', 'none', ...
            'FaceAlpha', 0.5, 'DisplayName', modes{m});
    end
    xline(1.96, 'k--', 'p<0.05');
    legend;
    xlabel('STI (Z-score)');
    ylabel('Proportion');
    title('Effect of Shuffle Mode on STI Distribution');
    box off

    % Print % significant per mode
    for m = 1:3
        flat = STI_modes{m}(:);
        fprintf('  %s: %.1f%% neurons Z>1.96\n', modes{m}, 100*mean(flat > 1.96, 'omitnan'))
    end
end

fprintf('\n=== All done ===\n')