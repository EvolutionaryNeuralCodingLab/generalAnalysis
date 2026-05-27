function result = detectStripeRF(rfImage, options)
% detectStripeRF  Detect diagonal excitatory stripe patterns in a 2D RF.
%
%   Uses the Radon transform to detect whether the RF contains a stripe
%   (ridge/line) pattern.  Four filters must all pass for isStripe=true:
%       (i)   the anisotropy score is significant vs. a pixel-permutation
%             null AND exceeds an optional absolute floor (minStripeScore)
%       (ii)  the stripe orientation lies in the requested diagonal
%             direction (BLtoTR by default)
%       (iii) the stripe orientation is at least minAngleFromHorizontal
%             degrees away from horizontal AND vertical
%       (iv)  the stripe is excitatory (peak of best projection is positive)
%
% ------------------------------------------------------------------------
%   Stripeness metric (anisotropy):
%       S = max(var_across_angles) / min(var_across_angles)
%   For a stripe, perpendicular projections are peaked (high variance)
%   while parallel projections are flat (low variance) — so S >> 1.
%   For blobs and noise, all projections look similar — so S ≈ 1.
% ------------------------------------------------------------------------
%   Angle conventions:
%       stripeAngle ∈ [0, 180), measured counter-clockwise from horizontal
%         0°   → horizontal
%         45°  → bottom-left → top-right (BL→TR)
%         90°  → vertical
%         135° → bottom-right → top-left (BR→TL)
%
%       projAngle (in the Radon variance profile) is perpendicular:
%         stripeAngle = mod(projAngle + 90, 180)
% ------------------------------------------------------------------------
%
%   INPUTS
%       rfImage    2D matrix (e.g. 54×54 for MB, 9×9 for RG).
%                  Should be shuffle-subtracted so values are signed.
%
%   NAME-VALUE OPTIONS
%       angleStep               degrees between tested projection angles (5)
%       nSurrogates             pixel-permutation surrogates (1000)
%       alpha                   significance threshold (0.05)
%       minStripeScore          absolute floor on S (2)
%       rngSeed                 RNG seed for reproducibility (42; NaN = skip)
%       diagonalOnly            apply orientation filter (true)
%       diagonalDirection       "BLtoTR" | "BRtoTL" | "both"   default "BLtoTR"
%       minAngleFromHorizontal  degrees of margin from cardinal axes (3)
%                               accepted band is (minAng, 90 − minAng)
%                               for each chosen diagonal
%       requirePositive         require excitatory peak (true)
%
%   OUTPUT
%       result struct with:
%           isStripe         logical — all filters pass
%           passedScore      logical
%           passedDiagonal   logical
%           passedPositive   logical
%           stripeScore      double  — max(var)/min(var)
%           stripePval       double
%           stripeAngle      double  — degrees, 0–180
%           projAngle        double  — Radon angle at peak variance
%           peakSign         double  — +1 / −1
%           varProfile       [1×nAngles]
%           angles           [1×nAngles]
%           nullScores       [nSurrogates×1]

arguments
    rfImage       (:,:) double
    options.angleStep              (1,1) double  = 5
    options.nSurrogates            (1,1) double  = 1000
    options.alpha                  (1,1) double  = 0.05
    options.minStripeScore         (1,1) double  = 2
    options.rngSeed                (1,1) double  = 42
    options.diagonalOnly           (1,1) logical = true
    options.diagonalDirection      (1,1) string  = "BLtoTR"   % "BLtoTR", "BRtoTL", or "both"
    options.minAngleFromHorizontal (1,1) double  = 3
    options.requirePositive        (1,1) logical = true
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

% Stripe orientation (perpendicular to peak projection angle)
stripeAngle = mod(bestProjAngle + 90, 180);

% (ii)+(iii) Diagonal direction with margin from horizontal/vertical
if options.diagonalOnly
    minA = options.minAngleFromHorizontal;

    % BL→TR band: stripeAngle ∈ (minA, 90 − minA)
    inBLtoTR = (stripeAngle > minA)        && (stripeAngle < 90 - minA);
    % BR→TL band: stripeAngle ∈ (90 + minA, 180 − minA)
    inBRtoTL = (stripeAngle > 90 + minA)   && (stripeAngle < 180 - minA);

    switch options.diagonalDirection
        case "BLtoTR"
            passedDiagonal = inBLtoTR;
        case "BRtoTL"
            passedDiagonal = inBRtoTL;
        case "both"
            passedDiagonal = inBLtoTR || inBRtoTL;
        otherwise
            error('Unknown diagonalDirection: %s. Use "BLtoTR", "BRtoTL", or "both".', ...
                  options.diagonalDirection);
    end
else
    passedDiagonal = true;
end

% (iv) Excitatory peak
[~, peakIdx] = max(abs(bestProj));
peakValue    = bestProj(peakIdx);
peakSign     = sign(peakValue);
if options.requirePositive
    passedPositive = peakValue > 0;
else
    passedPositive = true;
end

% Final classification
isStripe = passedScore && passedDiagonal && passedPositive;

% -------------------------------------------------------------------------
% 5.  Package results
% -------------------------------------------------------------------------
result.isStripe        = isStripe;
result.passedScore     = passedScore;
result.passedDiagonal  = passedDiagonal;
result.passedPositive  = passedPositive;
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
%  LOCAL FUNCTION
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