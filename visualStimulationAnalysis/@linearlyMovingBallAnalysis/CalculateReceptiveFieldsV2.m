function results = CalculateReceptiveFieldsV2(obj, params)

arguments (Input)
    obj
    params.overwrite            logical  = false
    params.analysisTime                  = datetime('now')
    params.inputParams                   = false
    params.preBase                       = 200
    params.bin                           = 15
    params.exNeurons                     = 0
    params.AllResponsiveNeurons logical  = true
    params.fixedWindow          logical  = false
    params.speed                         = 1
    params.noEyeMoves           logical  = false
    params.delay                         = 250
    params.nShuffle                      = 100
    params.testConvolution      logical  = false
    params.reduceFactor                  = 20
    params.useGPU               logical  = false
    params.shuffleMode          string   = "within"  % 'within', 'global', 'within_bins_only'
end

if params.inputParams, disp(params); return; end

%% -----------------------------------------------------------------------
% 0. File check
% -----------------------------------------------------------------------
filename = obj.getAnalysisFileName;
[filepath, name, ~] = fileparts(filename);
filename = fullfile(filepath, name);
savePath = sprintf('%s-Speed-%d-V2.mat', filename, params.speed);

if isfile(savePath) && ~params.overwrite
    if nargout == 1
        fprintf('Loading saved results from file.\n');
        results = load(savePath);
    else
        fprintf('Analysis already exists (use overwrite to recalculate).\n');
    end
    return
end

%% -----------------------------------------------------------------------
% 1. Load neurons and spike sorting
% -----------------------------------------------------------------------
NeuronResp = obj.ResponseWindow;
Stats      = obj.ShufflingAnalysis;
goodU      = NeuronResp.goodU;
p_spk      = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);

fieldName  = sprintf('Speed%d', params.speed);
pvals      = Stats.(fieldName).pvalsResponse;
C          = NeuronResp.(fieldName).C;
stimDur    = NeuronResp.(fieldName).stimDur;

if params.AllResponsiveNeurons
    respU = find(pvals < 0.05);
    if isempty(respU)
        fprintf('No responsive neurons.\n'); return
    end
end
if params.exNeurons > 0,   respU = params.exNeurons; end
if params.testConvolution, respU = 1;                end

nNeurons = numel(respU);
fprintf('Processing %d responsive neurons.\n', nNeurons);

%% -----------------------------------------------------------------------
% 2. Speed and trajectory setup
% -----------------------------------------------------------------------
selecFrames = sort(unique(obj.VST.nFrames), 'descend');
try
    selecFrames = selecFrames(params.speed);
    selecSpeed  = obj.VST.speed(params.speed);
    speedIndx   = params.speed;
catch
    selecFrames = selecFrames(end);
    selecSpeed  = obj.VST.speed(end);
    speedIndx   = numel(obj.VST.speed);
end
fprintf('Selecting speed of %d frames/sec.\n', selecSpeed)

Xpos = obj.VST.ballTrajectoriesX;
Ypos = obj.VST.ballTrajectoriesY;
if size(Xpos, 1) > 1
    Xpos = Xpos(speedIndx, :, :, 1:selecFrames);
    Ypos = Ypos(speedIndx, :, :, 1:selecFrames);
end

sizeN = numel(unique(obj.VST.ballSizes));
sizeX = size(Xpos);   % [speeds x offsets x directions x frames]

%% -----------------------------------------------------------------------
% 3. Sort trial matrix C
% -----------------------------------------------------------------------
% Original C columns: [stimOn directions offsets sizes speeds orientations]
% Reorganise to:      [stimOn directions offsets speeds  sizes  orientations]
trialDivisionVid = size(C,1) / numel(unique(C(:,2))) / numel(unique(C(:,3))) ...
                             / numel(unique(C(:,4))) / numel(unique(C(:,5)));

C = [C(:,1) C(:,2) C(:,3) C(:,5) C(:,4) C(:,6)];
[C, ~]           = sortrows(C, [2 3 4 5 6]);
directimesSorted = C(:,1)';
sizeV            = C(:,5);

trialDivision = size(C,1) / numel(unique(C(:,2))) / numel(unique(C(:,3))) ...
                           / numel(unique(C(:,4))) / numel(unique(C(:,5))) ...
                           / numel(unique(C(:,6)));

offsetN = numel(unique(C(:,3)));

%% -----------------------------------------------------------------------
% 4. Screen grid
% -----------------------------------------------------------------------
coorRect     = obj.VST.rect;
reduceFactor = min([params.reduceFactor, min(sizeV)]);
redCoorX     = round(coorRect(3) / reduceFactor);
redCoorY     = round(coorRect(4) / reduceFactor);
[x, y]       = meshgrid(1:redCoorX, 1:redCoorY);
x            = fliplr(x);
xCol         = (redCoorX - redCoorY)/2 + 1 : (redCoorX - redCoorY)/2 + redCoorY;

fprintf('Reduced screen size: %d x %d pixels.\n', redCoorX, redCoorY)

%% -----------------------------------------------------------------------
% 5. Build spike raster
% -----------------------------------------------------------------------
msPerFrame = stimDur / sizeX(4);

Raster = BuildBurstMatrix(goodU(:,respU), ...
    round(p_spk.t   / msPerFrame), ...
    round(directimesSorted / msPerFrame), ...
    round(stimDur   / msPerFrame));          % [nT x nNeurons x nBins]

[nT, ~, nB] = size(Raster);
spikeSum    = single(Raster) ./ msPerFrame .* 1000;   % [nT x nNeurons x nBins] sp/s

fprintf('Raster size: %d trials x %d neurons x %d bins.\n', nT, nNeurons, nB)

%% -----------------------------------------------------------------------
% 6. Vectorised video pre-computation
%    videoTrials: [nUniqueTrials x redCoorY x redCoorY x nFrames]
% -----------------------------------------------------------------------
% Build ChangePosX/Y: [totalTrials x nFrames]
ChangePosX = zeros(sizeX(1)*sizeX(2)*sizeX(3)*sizeN*trialDivisionVid, sizeX(4));
ChangePosY = zeros(size(ChangePosX));
j = 1;
for d  = 1:sizeX(3)    % directions
    for of = 1:sizeX(2) % offsets
        for sp = 1:sizeX(1) % speeds
            rowX = squeeze(Xpos(sp,of,d,:))';
            rowY = squeeze(Ypos(sp,of,d,:))';
            ChangePosX(j:j+sizeN*trialDivisionVid-1, :) = repmat(rowX, sizeN*trialDivisionVid, 1);
            ChangePosY(j:j+sizeN*trialDivisionVid-1, :) = repmat(rowY, sizeN*trialDivisionVid, 1);
            j = j + sizeN*trialDivisionVid;
        end
    end
end

% Vectorised mask generation across all unique trials and frames
nUniqueTrials = nT / trialDivisionVid;
videoTrials   = zeros(nUniqueTrials, redCoorY, redCoorY, sizeX(4), 'single');

fprintf('Pre-computing video frames...\n')
for f = 1:sizeX(4)
    cx    = ChangePosX(1:trialDivisionVid:nT, f) / reduceFactor;   % [nUnique x 1]
    cy    = ChangePosY(1:trialDivisionVid:nT, f) / reduceFactor;
    radii = sizeV(1:trialDivisionVid:nT) / 2 / reduceFactor + 0.5; % [nUnique x 1]

    % Broadcasting: [Y x X x 1] vs [1 x 1 x nUnique]
    dx   = reshape(x,     redCoorY, redCoorX, 1) - reshape(cx,    1, 1, []);
    dy   = reshape(y,     redCoorY, redCoorX, 1) - reshape(cy,    1, 1, []);
    dist = sqrt(dx.^2 + dy.^2);                                     % [Y x X x nUnique]
    mask = dist <= reshape(radii, 1, 1, []);                        % logical [Y x X x nUnique]
    mask = permute(mask(:, xCol, :), [3 1 2]);                      % [nUnique x Y x Y]
    videoTrials(:,:,:,f) = single(mask);
end
fprintf('Video pre-computation complete.\n')

%% -----------------------------------------------------------------------
% 7. Unique condition labels
% -----------------------------------------------------------------------
Udir  = unique(C(:,2));
Usize = unique(C(:,5));
Ulum  = unique(C(:,6));
nDir  = numel(Udir);
nSize = numel(Usize);
nLum  = numel(Ulum);

fprintf('Conditions — Directions: %d | Sizes: %d | Luminosities: %d\n', nDir, nSize, nLum)

%% -----------------------------------------------------------------------
% 8. Pre-allocate accumulators
% -----------------------------------------------------------------------
RFuDirSizeLum = zeros(nDir, nSize, nLum, redCoorY, redCoorY, sizeX(4), nNeurons, 'single');
RFu           = zeros(redCoorY, redCoorY, sizeX(4), nNeurons, 'single');
RFuShuffSum   = zeros(redCoorY, redCoorY, sizeX(4), nNeurons, params.nShuffle, 'single');

normFactor    = trialDivision / nT;

%% -----------------------------------------------------------------------
% 9. Validate shuffle mode
% -----------------------------------------------------------------------
validModes = ["within", "global", "within_bins_only"];
if ~any(params.shuffleMode == validModes)
    error('Invalid shuffleMode "%s". Choose: within, global, within_bins_only.', params.shuffleMode)
end
fprintf('Shuffle mode: %s | nShuffle: %d\n', params.shuffleMode, params.nShuffle)

% For global mode only: pre-generate shuffle indices once (reused across blocks)
if params.shuffleMode == "global"
    globalIdx1 = zeros(nT, params.nShuffle, 'uint32');
    globalIdx3 = zeros(nB, params.nShuffle, 'uint32');
    for s = 1:params.nShuffle
        globalIdx1(:,s) = randperm(nT);
        globalIdx3(:,s) = randperm(nB);
    end
else
    globalIdx1 = [];   % not used
    globalIdx3 = [];
end

%% -----------------------------------------------------------------------
% 10. Main convolution loop
% -----------------------------------------------------------------------
fprintf('Starting convolution...\n')
tic

nBlocks = nT / trialDivision;
p_iter  = 1;

for i = 1:trialDivision:nT

    % --- Condition indices for this block ---
    dirIdx  = find(Udir  == C(i,2));
    sizeIdx = find(Usize == C(i,5));
    lumIdx  = find(Ulum  == C(i,6));

    % --- Video for this block ---
    vidIdx = ceil(p_iter / 2);
    vid    = squeeze(videoTrials(vidIdx,:,:,:));    % [Y x Y x F]

    % --- Real spike mean: [nNeurons x nBins] ---
    spkBlock = spikeSum(i:i+trialDivision-1, :, :);         % [tDiv x nN x F]
    spkMean  = squeeze(mean(spkBlock, 1, 'omitnan'));        % [nN x F]

    needsFlip = (C(i,2) == Udir(1)) || (C(i,2) == Udir(3));

    % --- Convolve real RF: all neurons vectorised ---
    Co = convVideoSpikes(vid, spkMean, needsFlip);           % [Y x Y x F x nN]

    % --- Accumulate real RF ---
    RFuDirSizeLum(dirIdx, sizeIdx, lumIdx, :,:,:,:) = ...
        squeeze(RFuDirSizeLum(dirIdx, sizeIdx, lumIdx, :,:,:,:)) + Co * normFactor;
    RFu = RFu + Co * normFactor;

    % --- Shuffled spike means: [nN x F x nShuffle] ---
    spkMeanShuff = computeShuffleMeans(Raster, i, trialDivision, ...
                                       params.nShuffle, msPerFrame, ...
                                       params.shuffleMode, ...
                                       globalIdx1, globalIdx3);

    % --- Convolve shuffled RF: all neurons and shuffles vectorised ---
    CoShuff = convVideoSpikesShuff(vid, spkMeanShuff, needsFlip); % [Y x Y x F x nN x nShuffle]

    % --- Accumulate shuffled RF ---
    RFuShuffSum = RFuShuffSum + CoShuff * normFactor;

    if mod(p_iter, 5) == 0
        fprintf('  Block %d/%d (%.1f%%)\n', p_iter, nBlocks, 100*p_iter/nBlocks);
    end
    p_iter = p_iter + 1;
end
toc

%% -----------------------------------------------------------------------
% 11. Select delay time point
% -----------------------------------------------------------------------
L             = size(spikeSum, 3);
time_zero_idx = ceil(L / 2);
delay_idx     = time_zero_idx + round(params.delay / msPerFrame);

% [nDir x nSize x nLum x Y x Y x nNeurons]
RFuSTDirSizeLum = squeeze(RFuDirSizeLum(:,:,:,:,:,delay_idx,:));

% [Y x Y x nNeurons]
RFuST = squeeze(RFu(:,:,delay_idx,:));

% [Y x Y x nNeurons x nShuffle]
RFuShuffST = squeeze(RFuShuffSum(:,:,delay_idx,:,:));

%% -----------------------------------------------------------------------
% 12. Gaussian smoothing — parfor over neurons (outermost dim)
% -----------------------------------------------------------------------
gaussSigma = redCoorY / offsetN;
gaussSize  = floor(redCoorY / (offsetN / 2));
TwoDGauss  = fspecial('gaussian', gaussSize, gaussSigma);

[d1,d2,d3,d4,d5,~] = size(RFuSTDirSizeLum);
RFuDirSizeLumFilt   = zeros(size(RFuSTDirSizeLum), 'single');

% Also smooth shuffled RFs: [Y x Y x nNeurons x nShuffle]
[sY1, sY2, ~, nS] = size(RFuShuffST);
RFuShuffSTFilt     = zeros(sY1, sY2, nNeurons, nS, 'single');

fprintf('Applying Gaussian smoothing...\n')

parfor ui = 1:nNeurons

    % --- Smooth per-condition RF ---
    tmp = zeros(d1, d2, d3, d4, d5, 'single');
    for d = 1:d1
        for s = 1:d2
            for l = 1:d3
                slice        = squeeze(RFuSTDirSizeLum(d,s,l,:,:,ui));
                tmp(d,s,l,:,:) = conv2(single(slice), TwoDGauss, 'same');
            end
        end
    end
    RFuDirSizeLumFilt(:,:,:,:,:,ui) = tmp;

    % --- Smooth shuffled RF ---
    tmpShuff = zeros(sY1, sY2, nS, 'single');
    for s = 1:nS
        tmpShuff(:,:,s) = conv2(squeeze(RFuShuffST(:,:,ui,s)), TwoDGauss, 'same');
    end
    RFuShuffSTFilt(:,:,ui,:) = tmpShuff;

end

fprintf('Smoothing complete.\n')

%% -----------------------------------------------------------------------
% 13. Save and return
% -----------------------------------------------------------------------
fprintf('Saving results to file.\n');
S.params              = params;
S.RFuDirSizeLumFilt   = RFuDirSizeLumFilt;   % [nDir x nSize x nLum x Y x Y x nNeurons]  — smoothed
S.RFuSTDirSizeLum     = RFuSTDirSizeLum;     % [nDir x nSize x nLum x Y x Y x nNeurons]  — unsmoothed
S.RFuST               = RFuST;               % [Y x Y x nNeurons]                         — all conditions collapsed
S.RFuShuffST          = RFuShuffST;          % [Y x Y x nNeurons x nShuffle]              — unsmoothed shuffled
S.RFuShuffSTFilt      = RFuShuffSTFilt;      % [Y x Y x nNeurons x nShuffle]              — smoothed shuffled
save(savePath, '-struct', 'S');
results = S;

fprintf('Done.\n')

end   % <<<< end of main function


%% =========================================================================
%  LOCAL HELPER FUNCTIONS
%% =========================================================================

function spkMeanShuff = computeShuffleMeans(Raster, iStart, tDiv, nShuffle, ...
                                             msPerFrame, shuffleMode, ...
                                             globalIdx1, globalIdx3)
% Compute shuffled spike means for one trial block.
%
% Inputs:
%   Raster       [nT x nNeurons x nBins]
%   iStart       first row index of current block
%   tDiv         number of trials per block
%   nShuffle     number of shuffle iterations
%   msPerFrame   ms per bin (for sp/s conversion)
%   shuffleMode  'within' | 'global' | 'within_bins_only'
%   globalIdx1   [nT x nShuffle] uint32 — only used in 'global' mode
%   globalIdx3   [nBins x nShuffle] uint32 — only used in 'global' mode
%
% Output:
%   spkMeanShuff [nNeurons x nBins x nShuffle] single, in sp/s

nN          = size(Raster, 2);
F           = size(Raster, 3);
nT          = size(Raster, 1);
block_rows  = iStart : iStart + tDiv - 1;

spkMeanShuff = zeros(nN, F, nShuffle, 'single');

switch shuffleMode

    case 'within'
        % ------------------------------------------------------------------
        % Shuffle trial order AND bin order within the current block only.
        % Preserves: firing rate distribution per condition.
        % Destroys:  spatial relationship between spikes and ball position.
        % Recommended for spatial Z-score baseline.
        % Note: if tDiv is small (e.g. <10), increase nShuffle to 100-200
        % to get a stable null distribution from the limited pool.
        % ------------------------------------------------------------------
        for s = 1:nShuffle
            idx_trials           = block_rows(randperm(tDiv));
            idx_bins             = randperm(F);
            shuf_block           = Raster(idx_trials, :, idx_bins);   % [tDiv x nN x F]
            spkMeanShuff(:,:,s)  = squeeze(mean(single(shuf_block), 1)) ./ msPerFrame .* 1000;
        end

    case 'global'
        % ------------------------------------------------------------------
        % Shuffle trial order AND bin order across ALL trials.
        % Preserves: overall population firing rate.
        % Destroys:  both spatial structure AND condition selectivity.
        % Faster but will inflate Z-scores for condition-selective neurons.
        % Use only if neurons are not condition-selective, or for comparison.
        % ------------------------------------------------------------------
        for s = 1:nShuffle
            shuffled_rows        = globalIdx1(block_rows, s);         % rows drawn from full nT pool
            idx_bins             = globalIdx3(:, s)';
            shuf_block           = Raster(shuffled_rows, :, idx_bins);
            spkMeanShuff(:,:,s)  = squeeze(mean(single(shuf_block), 1)) ./ msPerFrame .* 1000;
        end

    case 'within_bins_only'
        % ------------------------------------------------------------------
        % Shuffle ONLY time bins; keep trial identity fixed within block.
        % Preserves: trial-to-trial firing rate variability + condition identity.
        % Destroys:  only the temporal/spatial relationship within each trial.
        % Most conservative null — use when tDiv is very small (e.g. 3-5 trials)
        % and you cannot afford to shuffle trial identity.
        % ------------------------------------------------------------------
        for s = 1:nShuffle
            idx_bins             = randperm(F);
            shuf_block           = Raster(block_rows, :, idx_bins);
            spkMeanShuff(:,:,s)  = squeeze(mean(single(shuf_block), 1)) ./ msPerFrame .* 1000;
        end

end
end


function Co = convVideoSpikes(vid, spkMean, needsFlip)
% FFT-based convolution of a video with spike kernels across all neurons.
%
% Inputs:
%   vid       [Y x Y x F]     single — stimulus video for one trial block
%   spkMean   [nNeurons x F]  single — mean spike response
%   needsFlip logical         — whether to flip spatial dims (replaces rot90(...,2))
%
% Output:
%   Co        [Y x Y x F x nNeurons] single

[Y1, Y2, F] = size(vid);
nN          = size(spkMean, 1);

% Zero-pad to next power of 2 for FFT efficiency
Flen = 2^nextpow2(2*F - 1);

% Transform along frame dimension
Vf = fft(vid,      Flen, 3);                   % [Y x Y x Flen]
Sf = fft(spkMean', Flen, 1);                   % [Flen x nN]

% Broadcast multiply across neurons
% Vf: [Y x Y x Flen x 1]  .* Sf: [1 x 1 x Flen x nN]
Vf4    = repmat(Vf,                  [1, 1,  1,    nN]);
Sf4    = repmat(reshape(Sf,[1,1,Flen,nN]), [Y1,Y2, 1,     1]);
result = real(ifft(Vf4 .* Sf4, [], 3));        % [Y x Y x Flen x nN]

% Trim to 'same' size
offset = floor(F / 2);
Co     = result(:, :, offset+1:offset+F, :);   % [Y x Y x F x nN]

if needsFlip
    Co = Co(end:-1:1, end:-1:1, :, :);         % flip both spatial dims (= rot90 x2)
end
end


function CoShuff = convVideoSpikesShuff(vid, spkMeanShuff, needsFlip)
% FFT-based convolution of a video with shuffled spike kernels.
% Loops over shuffles; each shuffle is vectorised over neurons.
%
% Inputs:
%   vid           [Y x Y x F]            single
%   spkMeanShuff  [nNeurons x F x nShuffle] single
%   needsFlip     logical
%
% Output:
%   CoShuff       [Y x Y x F x nNeurons x nShuffle] single

[Y1, Y2, F] = size(vid);
[nN, ~, nS] = size(spkMeanShuff);
CoShuff     = zeros(Y1, Y2, F, nN, nS, 'single');

Flen = 2^nextpow2(2*F - 1);
Vf   = fft(vid, Flen, 3);                      % [Y x Y x Flen] — compute once

for s = 1:nS
    Sf     = fft(squeeze(spkMeanShuff(:,:,s))', Flen, 1);     % [Flen x nN]
    Vf4    = repmat(Vf,                  [1, 1,  1,    nN]);
    Sf4    = repmat(reshape(Sf,[1,1,Flen,nN]), [Y1,Y2, 1,     1]);
    result = real(ifft(Vf4 .* Sf4, [], 3));
    offset = floor(F / 2);
    Co_s   = result(:, :, offset+1:offset+F, :);

    if needsFlip
        Co_s = Co_s(end:-1:1, end:-1:1, :, :);
    end
    CoShuff(:,:,:,:,s) = Co_s;
end
end