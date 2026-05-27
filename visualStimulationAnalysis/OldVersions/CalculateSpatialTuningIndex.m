function [spatialTuningIndex] = CalculateSpatialTuningIndex(obj, params)
arguments (Input)
    obj
    params.overwrite         logical  = false
    params.analysisTime               = datetime('now')
    params.inputParams                = false
    params.OneDirection  string       = "all"
    params.OneLuminosity string       = "all"
    params.OneSize       string       = "all"
    params.plotResult    logical      = false
    params.fieldName     string       = "default"   % fieldname for moving ball responses struct
end

%% -----------------------------------------------------------------------
% 1. Load receptive fields (real and shuffled)
% -----------------------------------------------------------------------
% Assumption:
%   RFs          -> real RF matrix (stimulus dimensions x X x Y x N)
%   RFs_shuffled -> shuffled RF matrix (stimulus dimensions x X x Y x N x nShuffles)
fprintf('Loading receptive fields...\n');
[RFs, RFs_shuffled] = obj.CalculateReceptiveFields();

%% -----------------------------------------------------------------------
% 2. Identify stimulus type and load parameters
% -----------------------------------------------------------------------
isStatic = contains(obj.stimName, 'rectGrid',           'IgnoreCase', true);
isMoving = contains(obj.stimName, 'linearlyMovingBall', 'IgnoreCase', true);

if ~isStatic && ~isMoving
    error('CalculateSpatialTuningIndex: unrecognized stimName "%s".', obj.stimName);
end

responses = obj.ResponseWindow;

if isStatic
    % Static ball layout: R x L x S x X x Y x N
    % Shuffled layout:    R x L x S x X x Y x N x nShuffles
    uSize = unique(responses.C(:,3));
    uLum  = unique(responses.C(:,4));

elseif isMoving
    % Moving ball layout: D x S x L x X x Y x N
    % Shuffled layout:    D x S x L x X x Y x N x nShuffles
    fn = params.fieldName;
    uDir  = unique(responses.(fn).C(:,2));
    uSize = unique(responses.(fn).C(:,4));
    uLum  = unique(responses.(fn).C(:,6));
end

%% -----------------------------------------------------------------------
% 3. Filter dimensions based on params
% -----------------------------------------------------------------------
[sizeIdx] = filterDimension(uSize, params.OneSize,  'Size');
[lumIdx ] = filterDimension(uLum,  params.OneLuminosity, 'Luminosity');

if isMoving
    [dirIdx] = filterDimension(uDir, params.OneDirection, 'Direction');
end

%% -----------------------------------------------------------------------
% 4. Get matrix dimensions
% -----------------------------------------------------------------------
szRF     = size(RFs);
nNeurons = szRF(end);       % N is always last dim of real RF

szShuff    = size(RFs_shuffled);
nShuffles  = szShuff(end);  % nShuffles is always last dim of shuffled RF

fprintf('Neurons: %d | Shuffle iterations: %d\n', nNeurons, nShuffles);

%% -----------------------------------------------------------------------
% 5. Pre-allocate output
% -----------------------------------------------------------------------
if isStatic
    nR     = szRF(1);
    nLums  = numel(lumIdx);
    nSizes = numel(sizeIdx);
    % Output: R x L x S x N  (Z-score)
    spatialTuningIndex = nan(nR, nLums, nSizes, nNeurons);

elseif isMoving
    nDirs  = numel(dirIdx);
    nSizes = numel(sizeIdx);
    nLums  = numel(lumIdx);
    % Output: D x S x L x N  (Z-score)
    spatialTuningIndex = nan(nDirs, nSizes, nLums, nNeurons);
end

%% -----------------------------------------------------------------------
% 6. Main computation loop
% -----------------------------------------------------------------------
fprintf('Computing spatial tuning index for stimulus: %s\n', obj.stimName);

for n = 1:nNeurons
    if isStatic
        for ri = 1:nR
            for li = 1:numel(lumIdx)
                L = lumIdx(li);
                for si = 1:numel(sizeIdx)
                    S = sizeIdx(si);

                    % Real RF map: [X x Y]
                    map_real = squeeze(RFs(ri, L, S, :, :, n));

                    % Shuffled RF maps: [X x Y x nShuffles]
                    map_shuff = squeeze(RFs_shuffled(ri, L, S, :, :, n, :));

                    spatialTuningIndex(ri, li, si, n) = computeZScore(map_real, map_shuff);
                end
            end
        end

    elseif isMoving
        for di = 1:numel(dirIdx)
            D = dirIdx(di);
            for si = 1:numel(sizeIdx)
                S = sizeIdx(si);
                for li = 1:numel(lumIdx)
                    L = lumIdx(li);

                    % Real RF map: [X x Y]
                    map_real = squeeze(RFs(D, S, L, :, :, n));

                    % Shuffled RF maps: [X x Y x nShuffles]
                    map_shuff = squeeze(RFs_shuffled(D, S, L, :, :, n, :));

                    spatialTuningIndex(di, si, li, n) = computeZScore(map_real, map_shuff);
                end
            end
        end
    end
end

%% -----------------------------------------------------------------------
% 7. Optional summary plot
% -----------------------------------------------------------------------
if params.plotResult
    figure;
    tiledlayout(1,2);

    nexttile;
    histogram(spatialTuningIndex(:), 40, 'Normalization', 'probability', ...
              'FaceColor', [0.2 0.5 0.8]);
    xlabel('Spatial Tuning Index (Z-score)');
    ylabel('Proportion');
    title(sprintf('Population — %s', obj.stimName));
    xline(0, 'k--', 'Chance');
    xline(median(spatialTuningIndex(:), 'omitnan'), 'r--', 'Median');

    nexttile;
    imagesc(squeeze(mean(spatialTuningIndex, 4, 'omitnan')));  % mean across neurons
    colorbar;
    title('Mean Z-score across neurons');
    xlabel('Size index'); ylabel('Lum / Response index');
end

fprintf('Done. Spatial tuning index computed for %d neurons.\n', nNeurons);
end

%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function zScore = computeZScore(map_real, map_shuff)
% Computes Z-score of spatial entropy of real RF map against
% the null distribution built from shuffled RF maps.
%
% map_real  : [X x Y]           - real receptive field map
% map_shuff : [X x Y x nShuff]  - shuffled receptive field maps

    % Entropy of the real map
    H_real = spatialEntropy(map_real);

    % Entropy of each shuffled map -> null distribution
    nShuff = size(map_shuff, 3);
    H_null = nan(1, nShuff);
    for k = 1:nShuff
        H_null(k) = spatialEntropy(squeeze(map_shuff(:,:,k)));
    end

    mu_null  = mean(H_null,  'omitnan');
    std_null = std(H_null,   'omitnan');

    if std_null == 0
        zScore = 0;  % no variability in null -> undefined, return 0
    else
        zScore = (H_real - mu_null) / std_null;
    end

    % Note on sign convention:
    % High positive Z: real map is MORE spatially dispersed than chance
    % High negative Z: real map is MORE spatially concentrated than chance
    %   -> a well-tuned (focal) RF will give a NEGATIVE Z-score.
    % You may want to flip the sign: zScore = -zScore;
    % so that higher values = more spatially tuned. See note below.
    zScore = -zScore;  % flipped: higher Z = more focal/tuned RF
end


function H = spatialEntropy(map)
% Computes normalized spatial entropy of a 2D pixel map.
% The map is treated as a probability distribution over pixel locations.
%
% Normalization: H is divided by log2(nPixels) so H in [0,1]:
%   H = 0  -> all activity concentrated in one pixel (perfectly focal RF)
%   H = 1  -> activity uniformly spread across all pixels (no tuning)

    vals = double(map(:));

    % Shift to non-negative (in case of negative values after smoothing)
    vals = vals - min(vals);

    total = sum(vals);
    if total == 0
        H = 1;  % undefined map -> treat as maximally diffuse
        return;
    end

    % Normalize to probability distribution
    p = vals / total;

    % Remove zeros (0*log(0) = 0 by convention but log(0) = -Inf)
    p = p(p > 0);

    % Shannon entropy
    H_raw = -sum(p .* log2(p));

    % Normalize by maximum possible entropy
    nPixels = numel(map(:));
    H_max   = log2(nPixels);

    H = H_raw / H_max;
end


function idx = filterDimension(uVals, paramVal, dimName)
% Returns indices into uVals based on a string filter parameter.
% paramVal = "all" returns all indices; otherwise parses as number.

    if paramVal == "all"
        idx = 1:numel(uVals);
    else
        numVal = str2double(paramVal);
        idx    = find(uVals == numVal);
        if isempty(idx)
            error('Requested %s "%s" not found in stimulus set.', dimName, paramVal);
        end
    end
end