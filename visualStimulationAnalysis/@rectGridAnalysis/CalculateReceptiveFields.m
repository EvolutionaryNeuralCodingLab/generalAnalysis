function results = CalculateReceptiveFields(obj,params)

arguments (Input)
    obj
    params.overwrite logical = false
    params.analysisTime = datetime('now')
    params.inputParams = false
    params.preBase = 200
    params.bin = 10
    params.exNeurons = 0
    params.AllResponsiveNeurons = true
    params.fixedWindow = false
    params.noEyeMoves = false
    params.delay = 250
    params.nShuffle = 20                % Number of shuffles to generate shuffled receptive fields
    params.testConvolution = false
    params.reduceFactor = 20;
    params.duration = 300;              % Response window (ms)
    params.durationOff = 3000;          % Off-response window (ms)
    params.offsetR = 50;                % Response after onset of stim (ms)
    params.TakeAllStimDur = true        % Use whole stim window for RF calculation
    params.statType string = "maxPermutationTest"
    params.nGrid = 9
end

if params.inputParams, disp(params), return, end

filename = obj.getAnalysisFileName;

% -------------------------------------------------------------------------
% Load from file if it exists and overwrite is false
% -------------------------------------------------------------------------
if isfile(obj.getAnalysisFileName) && ~params.overwrite
    if nargout == 1
        fprintf('Loading saved results from file.\n');
        results = load(filename);
    else
        fprintf('Analysis already exists (use overwrite option to recalculate).\n');
    end
    return
end

NeuronResp = obj.ResponseWindow;

% Select statistics struct for p-values based on statType parameter
if params.statType == "BootstrapPerNeuron"
    Stats = obj.BootstrapPerNeuron;
else
    Stats = obj.StatisticsPerNeuron;
end

% Extract spike-sorted unit data: phy IDs, labels, and spike train matrix
p       = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
phy_IDg = p.phy_ID(string(p.label') == 'good');  % phy IDs of good units
label   = string(p.label');
goodU   = p.ic(:, label == 'good');               % spike train matrix for good units

% Get p-values and trial metadata from response window
pvals    = Stats.pvalsResponse;
C        = NeuronResp.C;        % trial condition matrix: [stimOn, pos, size, lum, ...]
stimDur  = NeuronResp.stimDur;  % stimulus duration (ms)

% Select all statistically responsive neurons (p < 0.05)
if params.AllResponsiveNeurons
    respU = find(pvals < 0.05);
    if isempty(respU)
        fprintf('No responsive neurons.\n')
        return
    end
else
    respU = 1:size(goodU,2);
end

% Override with manually specified neuron indices if provided
if params.exNeurons > 0
    respU = params.exNeurons;
end

% Extract stimulus layout: positions, sizes, luminosities
seqMatrix = obj.VST.pos;
sizes     = obj.VST.tilingRatios;
uSize     = unique(sizes);
nSize     = length(uSize);
uLums     = unique(obj.VST.rectLuminosity(obj.VST.luminosities));
nLums     = length(uLums);

% Number of trial repetitions per unique condition (pos x size x lum)
trialDiv = length(seqMatrix) / length(unique(seqMatrix)) / nSize / nLums;

% Sorted stimulus onset times
directimesSorted = C(:, 1)';

% Use full stimulus duration as response window if TakeAllStimDur is set
if params.TakeAllStimDur
    params.offsetR   = 0;
    params.duration  = stimDur;
    params.durationOff = NeuronResp.stimInter;
end

% Use the shorter of on/off windows to keep matrix sizes consistent
durationMin = min([params.duration params.durationOff]);

% Build spike count matrices for on-response (Mr) and off-response (Mro)
[Mr]  = BuildBurstMatrix(goodU, round(p.t / params.bin), ...
    round((directimesSorted + params.offsetR) / params.bin), ...
    round(durationMin / params.bin));
[Mro] = BuildBurstMatrix(goodU, round(p.t / params.bin), ...
    round((directimesSorted + stimDur) / params.bin), ...
    round(durationMin / params.bin));

[nT, nN, NB] = size(Mr);  % nT=trials, nN=neurons, NB=time bins

% -------------------------------------------------------------------------
% Build shuffle distributions for both on and off responses
% Each shuffle randomly permutes trial order (dim 1) and time bin order
% (dim 3) independently, breaking stimulus-response associations
% -------------------------------------------------------------------------
nShuffle = params.nShuffle;

% On-response raster in spike/s, used as base for on-shuffles
RasterOn  = Mr  .* (1000 / params.bin);  % [nT, nN, NB]
% Off-response raster in spike/s, used as base for off-shuffles
RasterOff = Mro .* (1000 / params.bin);  % [nT, nN, NB]

% Pre-allocate shuffle arrays: [nT, nN, NB, nShuffle]
shuffledDataOn  = zeros(nT, nN, NB, nShuffle);
shuffledDataOff = zeros(nT, nN, NB, nShuffle);

for i = 1:nShuffle
    idx1 = randperm(nT);   % shuffle trial order
    idx3 = randperm(NB);   % shuffle time bin order within trial

    shuffledDataOn( :,:,:,i) = RasterOn( idx1, :, idx3);
    shuffledDataOff(:,:,:,i) = RasterOff(idx1, :, idx3);
end

% -------------------------------------------------------------------------
% NOTE: Original code used shuffledData (on only). We now use
% shuffledDataOn / shuffledDataOff to keep on and off shuffles separate
% and symmetric. gridSpikeRateShuff still uses shuffledDataOn for
% backward compatibility (on response only).
% -------------------------------------------------------------------------

if params.noEyeMoves
    % Eye movement exclusion path — not implemented (see commented code above)
else

    % Build per-condition mean spike rate: [2, nLums, nSize, nPos, nN, NB]
    % dim 1: on(1) / off(2) response
    % NOTE: dim order here is [onOff, nLums, nSize, ...] — not [onOff, nSize, nLums]
    % (matches MrC indexing: MrC(o, uLums==..., uSize==..., j, u, :))
    MrC = zeros(2, nLums, nSize, round(nT / trialDiv), nN, NB);

    % Stack on and off into a single 4D array for unified indexing
    MRtotal        = zeros(2, size(Mr,1), size(Mr,2), size(Mr,3));
    MRtotal(1,:,:,:) = Mr;   % on-response
    MRtotal(2,:,:,:) = Mro;  % off-response

    % Average over trialDiv repetitions for each unique condition
    for u = 1:size(goodU, 2)
        for o = 1:2
            j = 1;
            for i = 1:trialDiv:nT
                % Mean over trialDiv reps, convert to spikes/s
                meanR = mean(squeeze(MRtotal(o, i:i+trialDiv-1, u, :))) .* (1000 / params.bin);
                MrC(o, uLums == C(i,4), uSize == C(i,3), j, u, :) = meanR;
                j = j + 1;
            end
        end
    end

    % Average over time bins to get mean rate per condition: [2, nLums, nSize, nPos, nN]
    MrMean = mean(MrC, 6);

end

% -------------------------------------------------------------------------
% Build position-indexed video frames: circle mask for each stimulus position
% -------------------------------------------------------------------------
screenSide = obj.VST.rect;
screenRed  = screenSide(4) / params.reduceFactor;  % reduced screen resolution
[x, y]     = meshgrid(1:screenRed, 1:screenRed);

pxyScreen  = zeros(screenRed, screenRed);                          % cumulative position coverage
VideoScreen = zeros(screenRed, screenRed, size(C,1) / trialDiv);  % per-position stimulus mask

rectData = obj.VST.rectData;

% Store reduced-coordinate centres for each unique position (for grid binning)
XcStore = zeros(1, size(C,1) / trialDiv);
YcStore = zeros(1, size(C,1) / trialDiv);

j = 1;
for i = 1:trialDiv:length(C)
    xyScreen = zeros(screenRed, screenRed)';

    % Compute circle centre in reduced pixel coordinates
    Xc = round((rectData.X2{1,C(i,3)}(C(i,2)) - rectData.X1{1,C(i,3)}(C(i,2))) / 2) + ...
         rectData.X1{1,C(i,3)}(C(i,2));
    Xc = Xc / params.reduceFactor;

    Yc = round((rectData.Y4{1,C(i,3)}(C(i,2)) - rectData.Y1{1,C(i,3)}(C(i,2))) / 2) + ...
         rectData.Y1{1,C(i,3)}(C(i,2));
    Yc = Yc / params.reduceFactor;

    XcStore(j) = Xc;
    YcStore(j) = Yc;

    % Circle radius in reduced pixels
    r = round((rectData.X2{1,C(i,3)}(C(i,2)) - rectData.X1{1,C(i,3)}(C(i,2))) / 2) / params.reduceFactor;

    % Binary circle mask: 1 inside stimulus, 0 outside
    distances               = sqrt((x - Xc).^2 + (y - Yc).^2);
    xyScreen(distances <= r) = 1;

    VideoScreen(:,:,j) = xyScreen';
    pxyScreen          = pxyScreen + xyScreen;
    j = j + 1;
end

% -------------------------------------------------------------------------
% Spike rate grid map: bin trials into nGrid x nGrid spatial grid
% -------------------------------------------------------------------------
nGrid  = params.nGrid;

% Grid edges in reduced pixel coordinates
xEdges = linspace(0, screenSide(3) / params.reduceFactor, nGrid + 1);
yEdges = linspace(0, screenSide(4) / params.reduceFactor, nGrid + 1);

% [nGrid, nGrid, nN, 2(on/off), nSize, nLums]
gridSpikeRate      = zeros(nGrid, nGrid, nN, 2, nSize, nLums);
% [nGrid, nGrid, nN, nShuffle, 2(on/off), nSize, nLums]
gridSpikeRateShuff = zeros(nGrid, nGrid, nN, nShuffle, 2, nSize, nLums);
trialCount         = zeros(nGrid, nGrid, nSize, nLums);

jj = 1;
for i = 1:trialDiv:nT

    % Bin stimulus centre into grid cell
    xBin = discretize(XcStore(jj), xEdges);
    yBin = discretize(YcStore(jj), yEdges);

    if isnan(xBin) || isnan(yBin)
        jj = jj + 1;
        continue
    end

    sizeIdx = find(uSize == C(i,3));
    lumIdx  = find(uLums == C(i,4));

    trialCount(yBin, xBin, sizeIdx, lumIdx) = trialCount(yBin, xBin, sizeIdx, lumIdx) + 1;

    % Mean on/off rate over trialDiv reps, convert to spikes/s: [1,1,nN]
    onRate  = reshape(mean(mean(Mr( i:i+trialDiv-1,:,:), 1), 3) .* (1000/params.bin), [1 1 nN]);
    offRate = reshape(mean(mean(Mro(i:i+trialDiv-1,:,:), 1), 3) .* (1000/params.bin), [1 1 nN]);

    gridSpikeRate(yBin, xBin, :, 1, sizeIdx, lumIdx) = gridSpikeRate(yBin, xBin, :, 1, sizeIdx, lumIdx) + onRate;
    gridSpikeRate(yBin, xBin, :, 2, sizeIdx, lumIdx) = gridSpikeRate(yBin, xBin, :, 2, sizeIdx, lumIdx) + offRate;

    % Accumulate shuffle spike rates (on response only, for grid)
    for s = 1:nShuffle
        shuffRate = reshape(mean(mean(shuffledDataOn(i:i+trialDiv-1,:,:,s), 1), 3), [1 1 nN]);
        gridSpikeRateShuff(yBin, xBin, :, s, 1, sizeIdx, lumIdx) = ...
            gridSpikeRateShuff(yBin, xBin, :, s, 1, sizeIdx, lumIdx) + shuffRate;
    end

    jj = jj + 1;
end

% Normalize grid maps by trial count per cell
for si = 1:nSize
    for li = 1:nLums
        tc = max(trialCount(:,:,si,li), 1);  % [nGrid, nGrid] — avoid divide-by-zero
        for s = 1:nShuffle
            for oi = 1:2
                gridSpikeRateShuff(:,:,:,s,oi,si,li) = gridSpikeRateShuff(:,:,:,s,oi,si,li) ./ tc;
            end
        end
        for oi = 1:2
            gridSpikeRate(:,:,:,oi,si,li) = gridSpikeRate(:,:,:,oi,si,li) ./ tc;
        end
    end
end

% -------------------------------------------------------------------------
% RF via point multiplication: weight VideoScreen frames by mean spike rate
% VD:  [2, nLums, nSize, screenRed, screenRed, nPos, nN]  (after repmat)
% Res: [2, nLums, nSize, 1,         1,         nPos, nN]
% -------------------------------------------------------------------------
nPos = size(VideoScreen, 3);  % number of unique stimulus positions

% Reshape VideoScreen to broadcast across onOff/lum/size/neuron dims
VD = reshape(VideoScreen, [1, 1, 1, screenRed, screenRed, nPos]);
VD = repmat(VD, [1, 1, 1, 1, 1, 1, nN]);  % [1, 1, 1, screenRed, screenRed, nPos, nN]

NanPos        = isnan(MrMean);
MrMean(NanPos) = 0;

% Reshape MrMean to align position and neuron dims with VD
Res = reshape(MrMean, [2, nLums, nSize, 1, 1, nPos, nN]);

% Weighted average across positions: [2, nLums, nSize, screenRed, screenRed, nN]
RFu = reshape(mean(VD .* Res, 6), [2, nLums, nSize, screenRed, screenRed, nN]);

% -------------------------------------------------------------------------
% Shuffle RF: same computation as RFu but using shuffled spike rates
% Result: RFuShuffMean [2, nLums, nSize, screenRed, screenRed, nN]
%         averaged across nShuffle, directly comparable to RFu
% -------------------------------------------------------------------------

% Pre-compute mean shuffle rate per trial by averaging over time bins
% [nT, nN, nShuffle] — avoids recomputing inside the shuffle loop

shuffMeanRateOn  = reshape(mean(shuffledDataOn,  3), [nT, nN, nShuffle]);  % [nT, nN, nShuffle]
shuffMeanRateOff = reshape(mean(shuffledDataOff, 3), [nT, nN, nShuffle]);  % [nT, nN, nShuffle]

% Accumulate shuffle RFs across all nShuffle — [2, nLums, nSize, screenRed, screenRed, nN, nShuffle]
RFuShuffAll = zeros(2, nLums, nSize, screenRed, screenRed, nN, nShuffle);

for sh = 1:nShuffle

    % Build per-condition mean shuffle rate: [2, nLums, nSize, nPos, nN]
    MrMeanShuff_sh = zeros(2, nLums, nSize, nPos, nN);
    j = 1;
    for i = 1:trialDiv:nT
        li = find(uLums == C(i,4));
        si = find(uSize == C(i,3));

        % Average over trialDiv reps for this position
        MrMeanShuff_sh(1, li, si, j, :) = mean(shuffMeanRateOn( i:i+trialDiv-1, :, sh), 1);  % on
        MrMeanShuff_sh(2, li, si, j, :) = mean(shuffMeanRateOff(i:i+trialDiv-1, :, sh), 1);  % off
        j = j + 1;
    end

    % Set NaNs to zero before multiplication (same as real RF path)
    MrMeanShuff_sh(isnan(MrMeanShuff_sh)) = 0;

    % Reshape to align with VD: [2, nLums, nSize, 1, 1, nPos, nN]
    ResSh = reshape(MrMeanShuff_sh, [2, nLums, nSize, 1, 1, nPos, nN]);

    % Weighted average across positions: [2, nLums, nSize, screenRed, screenRed, nN]
    RFuShuffAll(:,:,:,:,:,:,sh) = reshape(mean(VD .* ResSh, 6), ...
        [2, nLums, nSize, screenRed, screenRed, nN]);

end

% Average shuffle RFs across shuffles: [2, nLums, nSize, screenRed, screenRed, nN]
% This is the shuffle baseline, directly comparable to RFu
RFuShuffMean = mean(RFuShuffAll, 7);

% -------------------------------------------------------------------------
% Apply 2D Gaussian smoothing to both RFu and RFuShuffMean
% Sigma and kernel size scale with screen resolution and position layout
% -------------------------------------------------------------------------
offsetN = sqrt(max(seqMatrix));  % number of positions along one screen axis

TwoDGaussian = fspecial('gaussian', ...
    floor(size(RFu, 4) / (offsetN / 2)), ...
    screenRed / offsetN);

RFuFilt          = zeros(size(RFu));       % smoothed RF
RFuShuffMeanFilt = zeros(size(RFu));       % smoothed shuffle RF baseline

for d = 1:size(RFu, 1)     % on/off
    for s = 1:size(RFu, 2)  % lums
        for l = 1:size(RFu, 3)  % sizes
            for ui = 1:size(RFu, 6)  % neurons
                slice                      = squeeze(RFu(d,s,l,:,:,ui));
                RFuFilt(d,s,l,:,:,ui)      = conv2(slice, TwoDGaussian, 'same');

                sliceSh                         = squeeze(RFuShuffMean(d,s,l,:,:,ui));
                RFuShuffMeanFilt(d,s,l,:,:,ui)  = conv2(sliceSh, TwoDGaussian, 'same');
            end
        end
    end
end

% -------------------------------------------------------------------------
% Save results
% -------------------------------------------------------------------------
S.RFu              = RFu;               % [2, nLums, nSize, screenRed, screenRed, nN]  — NOTE: dim2=lums, dim3=size
S.RFuFilt          = RFuFilt;           % Gaussian-smoothed version of RFu
S.RFuShuffMean     = RFuShuffMean;      % Shuffle baseline RF [same dims as RFu]
S.RFuShuffMeanFilt = RFuShuffMeanFilt;  % Gaussian-smoothed shuffle baseline
S.shuffledData     = shuffledDataOn;    % Kept for backward compatibility (on-response shuffles)
S.shuffledDataOff    = shuffledDataOff;
S.gridSpikeRate      = gridSpikeRate;      % [nGrid, nGrid, nN, 2, nSize, nLums]
S.gridSpikeRateShuff = gridSpikeRateShuff; % [nGrid, nGrid, nN, nShuffle, 2, nSize, nLums]
S.params           = params;

save(filename, '-struct', 'S');
fprintf('Saved results to %s\n', filename);

results = S;

end