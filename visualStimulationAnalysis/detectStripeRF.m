function result = detectStripeRF(rfImage, options)
% detectStripeRF  Detect diagonal excitatory stripe patterns in a 2D RF.
%
%   Filters that must all pass for isStripe = true:
%     (i)   anisotropy score significant vs pixel-permutation null
%           AND above minStripeScore
%     (ii)  stripe orientation within [minStripeAngle, maxStripeAngle]
%           in the requested diagonal direction
%     (iii) stripe is excitatory (positive projection peak)
%     (iv)  [opt-in, requireContiguity] the excitatory pixels form a
%           spatially contiguous, elongated cluster — not just a single
%           strong pixel. The Radon-variance ratio in (i) can be driven
%           by one outlier pixel interacting with background noise on a
%           coarse grid (e.g. the 9x9 rectGrid), independent of whether
%           the response is actually spatially elongated. See below.
%
%   Stripeness metric:  S = max(var) / min(var)   (anisotropy)
%
%   Angle convention (stripeAngle ∈ [0, 180)):
%       0°   = horizontal
%       45°  = bottom-left → top-right  (BL→TR)
%       90°  = vertical
%       135° = bottom-right → top-left  (BR→TL)
%
%   NAME-VALUE OPTIONS
%       angleStep           projection angle step (default 5°)
%       nSurrogates         pixel-permutation surrogates (default 1000)
%       alpha               significance threshold (default 0.05)
%       minStripeScore      absolute floor on S (default 2)
%       rngSeed             RNG seed (default 42; NaN to skip)
%       diagonalOnly        apply orientation filter (default true)
%       diagonalDirection   "BLtoTR" | "BRtoTL" | "both" (default "BLtoTR")
%       minStripeAngle      lower bound of accepted stripe angle (default 5°)
%       maxStripeAngle      upper bound of accepted stripe angle (default 60°)
%                           For BLtoTR these apply directly.
%                           For BRtoTL they are mirrored: (180-max, 180-min).
%       requirePositive     require excitatory peak (default true)
%       requireContiguity   apply the spatial-contiguity filter (iv) (default false)
%       contiguityThreshSD  active-pixel threshold, in SDs above the RF's
%                            own mean (default 1.5). Recomputed per RF, not
%                            a fixed value, so it scales with each image's
%                            own noise level.
%       minStripeCells       minimum size (pixels) of the largest connected
%                            excitatory cluster (default 3 — the geometric
%                            minimum needed to define a direction at all;
%                            2 pixels is just an edge that could point
%                            anywhere by chance).
%       minElongation        minimum major/minor axis ratio of that cluster
%                            (default 1.5 — excludes roughly circular/blob
%                            clusters, keeps elongated ones).
%       contiguityImage      optional separate image to run the contiguity
%                            check (iv) on instead of rfImage (default: use
%                            rfImage). MUST be supplied when rfImage has been
%                            smoothed/interpolated (e.g. a convolution-based
%                            RF map) — blurring manufactures spatial
%                            correlation regardless of whether the true,
%                            non-interpolated data actually forms a
%                            contiguous cluster, so checking contiguity on
%                            a smoothed image is close to a no-op filter.

arguments
    rfImage       (:,:) double
    options.angleStep          (1,1) double  = 5
    options.nSurrogates        (1,1) double  = 1000
    options.alpha              (1,1) double  = 0.05
    options.minStripeScore     (1,1) double  = 2
    options.rngSeed            (1,1) double  = 42
    options.diagonalOnly       (1,1) logical = true
    options.diagonalDirection  (1,1) string  = "BLtoTR"
    options.minStripeAngle     (1,1) double  = 20      % degrees — lower bound
    options.maxStripeAngle     (1,1) double  = 60     % degrees — upper bound
    options.requirePositive    (1,1) logical = true
    options.requireContiguity  (1,1) logical = false
    options.contiguityThreshSD (1,1) double  = 1.5
    options.minStripeCells     (1,1) double  = 3
    options.minElongation      (1,1) double  = 1.5
    options.contiguityImage    (:,:) double  = []   % optional: check contiguity on a
                                                      % different (typically raw, non-
                                                      % interpolated) image than rfImage.
                                                      % Required whenever rfImage has been
                                                      % smoothed/interpolated, since blurring
                                                      % manufactures spatial correlation
                                                      % regardless of the true underlying
                                                      % signal. Defaults to rfImage.
end

% -------------------------------------------------------------------------
% 0.  RNG
% -------------------------------------------------------------------------
if ~isnan(options.rngSeed)
    rng(options.rngSeed, 'twister');
end

% -------------------------------------------------------------------------
% 1.  Projection angles
% -------------------------------------------------------------------------
angles = 0 : options.angleStep : (180 - options.angleStep);

% -------------------------------------------------------------------------
% 2.  Stripeness on the real RF
% -------------------------------------------------------------------------
[stripeScore, bestProjAngle, varProfile, bestProj] = ...
    computeStripeness(rfImage, angles);

% -------------------------------------------------------------------------
% 3.  Pixel-permutation null
% -------------------------------------------------------------------------
nPix       = numel(rfImage);
imgSize    = size(rfImage);
nSurr      = options.nSurrogates;
nullScores = zeros(nSurr, 1);

for si = 1:nSurr
    surrogate      = rfImage(randperm(nPix));
    surrogate      = reshape(surrogate, imgSize);
    nullScores(si) = computeStripeness(surrogate, angles);
end

% -------------------------------------------------------------------------
% 4.  Filters
% -------------------------------------------------------------------------

% (i) Significance + absolute floor
pval        = mean(nullScores >= stripeScore);
passedScore = (pval < options.alpha) && (stripeScore >= options.minStripeScore);

% Stripe orientation
stripeAngle = mod(bestProjAngle + 90, 180);

% (ii) Angle within accepted range
if options.diagonalOnly
    lo = options.minStripeAngle;   % e.g. 5°
    hi = options.maxStripeAngle;   % e.g. 60°

    % BL→TR band:  stripeAngle ∈ [lo, hi]
    inBLtoTR = (stripeAngle >= lo) && (stripeAngle <= hi);
    % BR→TL band:  mirror → stripeAngle ∈ [180−hi, 180−lo]
    inBRtoTL = (stripeAngle >= (180 - hi)) && (stripeAngle <= (180 - lo));

    switch options.diagonalDirection
        case "BLtoTR"
            passedDiagonal = inBLtoTR;
        case "BRtoTL"
            passedDiagonal = inBRtoTL;
        case "both"
            passedDiagonal = inBLtoTR || inBRtoTL;
        otherwise
            error('Unknown diagonalDirection: %s', options.diagonalDirection);
    end
else
    passedDiagonal = true;
end

% (iii) Excitatory peak
[~, peakIdx] = max(abs(bestProj));
peakValue    = bestProj(peakIdx);
peakSign     = sign(peakValue);
if options.requirePositive
    passedPositive = peakValue > 0;
else
    passedPositive = true;
end

% (iv) Spatial contiguity — opt-in, independent of the Radon anisotropy
% metric above. See checkContiguity for the mechanism. Uses
% contiguityImage if the caller supplied one (e.g. the raw, non-
% interpolated grid), otherwise falls back to rfImage itself.
if options.requireContiguity
    if isempty(options.contiguityImage)
        imgForContiguity = rfImage;
    else
        imgForContiguity = options.contiguityImage;
    end
    [passedContiguity, clusterSize, clusterElongation] = checkContiguity( ...
        imgForContiguity, options.contiguityThreshSD, options.minStripeCells, options.minElongation);
else
    passedContiguity  = true;   % filter disabled: does not affect isStripe
    clusterSize       = NaN;
    clusterElongation = NaN;
end

% Final
isStripe = passedScore && passedDiagonal && passedPositive && passedContiguity;

% -------------------------------------------------------------------------
% 5.  Package
% -------------------------------------------------------------------------
result.isStripe          = isStripe;
result.passedScore       = passedScore;
result.passedDiagonal    = passedDiagonal;
result.passedPositive    = passedPositive;
result.passedContiguity  = passedContiguity;
result.clusterSize       = clusterSize;         % NaN if requireContiguity=false
result.clusterElongation = clusterElongation;   % NaN if requireContiguity=false
result.stripeScore     = stripeScore;
result.stripePval      = pval;
result.stripeAngle     = stripeAngle;
result.projAngle       = bestProjAngle;
result.peakSign        = peakSign;
result.varProfile      = varProfile;
result.angles          = angles;
result.nullScores      = nullScores;

end


% =========================================================================
function [score, bestAngle, varProfile, bestProj] = computeStripeness(img, angles)

R          = radon(img, angles);
varProfile = var(R, 0, 1);

minVar = min(varProfile);
maxVar = max(varProfile);

if minVar <= 0 || maxVar <= 0
    score     = 1;
    bestAngle = NaN;
    bestProj  = R(:, 1);
    return
end

score = maxVar / minVar;

[~, maxIdx] = max(varProfile);
bestAngle   = angles(maxIdx);
bestProj    = R(:, maxIdx);

end


% =========================================================================
function [passed, clusterSize, elongation] = checkContiguity(rfImage, threshSD, minCells, minElong)
% checkContiguity  Require the excitatory response to form a real,
% spatially contiguous, elongated cluster of pixels — not a single
% strong pixel. This is independent of (and complements) the Radon
% anisotropy ratio: a lone bright pixel on a coarse grid can still
% produce a high maxVar/minVar ratio by chance interaction with
% background noise, without the response actually being stripe-shaped.
%
% Mechanism:
%   1. Threshold rfImage at mean + threshSD*SD (recomputed per RF, so it
%      scales with each image's own noise level rather than a fixed value).
%   2. Label 8-connected components of that binary mask (diagonal
%      neighbors count as connected, since we're looking for diagonal
%      stripes).
%   3. Take the largest component and require it has >= minCells pixels
%      (below that, a "direction" isn't geometrically well-defined) and
%      an elongation (major/minor axis ratio) >= minElong (excludes
%      round/blob-shaped clusters).

thresh = mean(rfImage(:)) + threshSD * std(rfImage(:));
mask   = rfImage > thresh;

if ~any(mask(:))
    passed = false; clusterSize = 0; elongation = 0;
    return
end

cc = bwconncomp(mask, 8);                 % 8-connectivity: diagonal neighbors count
sizes = cellfun(@numel, cc.PixelIdxList);
[clusterSize, biggest] = max(sizes);

if clusterSize < minCells
    passed = false; elongation = 0;
    return
end

props = regionprops(cc, 'MajorAxisLength', 'MinorAxisLength');
majorAx = props(biggest).MajorAxisLength;
minorAx = props(biggest).MinorAxisLength;
if minorAx <= 0
    elongation = Inf;              % perfectly thin line — maximally elongated
else
    elongation = majorAx / minorAx;
end

passed = elongation >= minElong;

end