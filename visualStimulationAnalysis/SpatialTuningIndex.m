function results = SpatialTuningIndex(exList, params)

arguments
    exList   double
    params.stimTypes  (1,:) string  = ["linearlyMovingBall", "rectGrid"]
    params.topPercent double        = 10
    params.overwrite  logical       = false
    params.statType   string        = "maxPermuteTest"
    params.speed      double        = 1
    params.plot       logical       = true
    params.indexType  string        = "L_amplitude_diff"  % L_amplitude_diff, L_amplitude_ratio, L_geometric, L_combined
    params.onOff      double        = 1             % 1=on, 2=off (rectGrid only; ignored for linearlyMovingBall)
    params.sizeIdx    double        = 1
    params.lumIdx     double        = 1
    params.nBoot      double        = 10000
    params.yLegend    char          = 'Spatial Tuning Index'
    params.yMaxVis    double        = 1
    params.Alpha      double        = 0.4
    params.PaperFig   logical       = false
    params.useRF      logical       = true   % If true: use RF (convolution-based) maps for both stim types.
                                             % If false: use gridSpikeRate (trial-binned) maps.
                                             % Recommended: true for linearlyMovingBall (avoids Y-offset bug),
                                             % and true for rectGrid for cross-stimulus comparability.
    params.prefDir    logical       = true  % If true (requires useRF=true): use each neuron's preferred
                                             % direction RF for linearlyMovingBall instead of averaging
                                             % across all directions. Preferred direction is defined as the
                                             % direction with the highest spike rate in NeuronVals.
                                             % Avoids deflating the spatial tuning index by averaging over
                                             % non-preferred directions.
end

% Guard: prefDir requires useRF — it operates on RFuSTDirSizeLum which is
% only available in the RF path
if params.prefDir && ~params.useRF
    error('prefDir=true requires useRF=true. The preferred direction RF is only available in the RF path.');
end

% -------------------------------------------------------------------------
% Build save path
% -------------------------------------------------------------------------
NP_first = loadNPclassFromTable(exList(1));  % Load first experiment to extract file path

% Build path to combined analysis directory based on first stim type
switch params.stimTypes(1)
    case "rectGrid"
        vs_first = rectGridAnalysis(NP_first);
    case "linearlyMovingBall"
        vs_first = linearlyMovingBallAnalysis(NP_first);
end

% Extract base path up to 'lizards' folder
p = extractBefore(vs_first.getAnalysisFileName, 'lizards');
p = [p 'lizards'];

% Create combined analysis directory if it does not exist
if ~exist([p '\Combined_lizard_analysis'], 'dir')
    cd(p)
    mkdir Combined_lizard_analysis
end
saveDir = [p '\Combined_lizard_analysis'];

% Build filename encoding experiment range, stim types, RF mode, and prefDir mode
% so that different parameter combinations never share a cache file
stimLabel  = strjoin(params.stimTypes, '-');
rfLabel    = '';
if params.useRF,   rfLabel   = '_RF';      end
prefLabel  = '';
if params.prefDir, prefLabel = '_prefDir'; end
nameOfFile = sprintf('\\Ex_%d-%d_SpatialTuningIndex_%s%s%s.mat', ...
    exList(1), exList(end), stimLabel, rfLabel, prefLabel);

% -------------------------------------------------------------------------
% Decide whether to compute or load from cache
% -------------------------------------------------------------------------
if exist([saveDir nameOfFile], 'file') == 2 && ~params.overwrite
    S = load([saveDir nameOfFile]);
    if isequal(S.expList, exList)
        fprintf('Loading saved SpatialTuningIndex from:\n  %s\n', [saveDir nameOfFile]);
        tbl       = S.tbl;
        goto_plot = true;
    else
        fprintf('Experiment list mismatch — recomputing.\n');
        goto_plot = false;
    end
else
    goto_plot = false;
end

% =========================================================================
% COMPUTE
% =========================================================================
if ~goto_plot

    nExp  = numel(exList);
    nStim = numel(params.stimTypes);

    % Guard: useRF must apply to all stim types — mixed inputs are not allowed
    % because RF (convolution-based) and gridSpikeRate (trial-binned) measures
    % are not directly comparable and must not be mixed across stim types
    if params.useRF
        supportedRF = ["rectGrid", "linearlyMovingBall"];
        unsupported = params.stimTypes(~ismember(params.stimTypes, supportedRF));
        if ~isempty(unsupported)
            error(['useRF=true is not supported for stim type(s): %s.\n' ...
                   'Either add RF support for these types or set useRF=false.'], ...
                   strjoin(unsupported, ', '));
        end
        if params.prefDir
            fprintf('RF mode (preferred direction): using preferred-direction RF for linearlyMovingBall.\n');
        else
            fprintf('RF mode: using convolution-based RF maps for all stim types.\n');
        end
    else
        fprintf('Grid mode: using trial-binned gridSpikeRate for all stim types.\n');
    end

    tbl = table();  % Accumulates one row per neuron x condition x stim type x experiment

    for ei = 1:nExp

        ex = exList(ei);
        fprintf('\n=== Experiment %d ===\n', ex);

        % Load NeuropixelsClass object for this experiment
        try
            NP = loadNPclassFromTable(ex);
        catch ME
            warning('Could not load experiment %d: %s', ex, ME.message);
            continue
        end

        % Use linearlyMovingBall to extract spike sorting info (shared across stim types)
        obj_s = linearlyMovingBallAnalysis(NP);

        % Extract animal name from recording name (first underscore-delimited token)
        nameParts  = split(NP.recordingName, '_');
        animalName = nameParts{1};

        % ----------------------------------------------------------
        % Get phy IDs for all good units (same spike sorting for all stim types)
        % ----------------------------------------------------------
        p_s     = NP.convertPhySorting2tIc(obj_s.spikeSortingFolder);
        phy_IDg = p_s.phy_ID(string(p_s.label') == 'good');  % phy IDs of good units

        % Stores responsive unit indices and phy IDs per stim type
        respPhyIDs_all = cell(1, nStim);
        respU_all      = cell(1, nStim);

        % ----------------------------------------------------------
        % Find responsive neurons for each stim type
        % ----------------------------------------------------------
        for s = 1:nStim
            stimType = params.stimTypes(s);
            try
                switch stimType
                    case "rectGrid"
                        obj_s = rectGridAnalysis(NP);
                    case "linearlyMovingBall"
                        obj_s = linearlyMovingBallAnalysis(NP);
                end

                % Select statistical test output
                if params.statType == "BootstrapPerNeuron"
                    Stats = obj_s.BootstrapPerNeuron;
                elseif params.statType == "maxPermuteTest"
                     Stats = obj_s.StatisticsPerNeuron;
                else
                    Stats = obj_s.ShufflingAnalysis;
                end

                % Extract p-values (linearlyMovingBall has per-speed fields)
                try
                    switch stimType
                        case "linearlyMovingBall"
                            fieldName = sprintf('Speed%d', params.speed);
                            pvals = Stats.(fieldName).pvalsResponse;
                        otherwise
                            pvals = Stats.pvalsResponse;
                    end
                catch
                    pvals = Stats.pvalsResponse;
                end

                % Indices of significantly responsive neurons (into phy_IDg)
                respU             = find(pvals < 0.05);
                respU_all{s}      = respU;
                respPhyIDs_all{s} = phy_IDg(respU);
                fprintf('  [%s] %d responsive neuron(s).\n', stimType, numel(respU));

            catch ME
                warning('Could not get pvals for %s exp %d: %s', stimType, ex, ME.message);
                respU_all{s}      = [];
                respPhyIDs_all{s} = [];
            end
        end

        % ----------------------------------------------------------
        % Keep only neurons responsive to ALL stim types (intersection)
        % ----------------------------------------------------------
        sharedPhyIDs = respPhyIDs_all{1};
        for s = 2:nStim
            sharedPhyIDs = intersect(sharedPhyIDs, respPhyIDs_all{s});
        end

        if isempty(sharedPhyIDs)
            fprintf('  No neurons responsive to all stim types in exp %d — skipping.\n', ex);
            continue
        end

        fprintf('  %d neuron(s) responsive to all stim types in exp %d.\n', numel(sharedPhyIDs), ex);

        % ----------------------------------------------------------
        % Loop over stim types and compute spatial tuning index
        % ----------------------------------------------------------
        for s = 1:nStim

            stimType = params.stimTypes(s);

            % Flag to track whether neuronIdx was already applied inside
            % the RF block (prefDir path) to avoid double-indexing
            alreadyIndexed = false;

            % Build stimulus-specific analysis object
            try
                switch stimType
                    case "rectGrid"
                        obj = rectGridAnalysis(NP);
                    case "linearlyMovingBall"
                        obj = linearlyMovingBallAnalysis(NP);
                    otherwise
                        error('Unknown stimType: %s', stimType);
                end
            catch ME
                warning('Could not build %s for exp %d: %s', stimType, ex, ME.message);
                continue
            end

            % Load pre-computed receptive field results from file
            S_rf = obj.CalculateReceptiveFields;

            % ----------------------------------------------------------
            % Build gridSpikeRateSelected and gridShuffMean
            % Two paths: RF-based (convolution) or grid-binned (gridSpikeRate)
            % ----------------------------------------------------------

            if params.useRF

                switch stimType

                    case "linearlyMovingBall"
                        % -------------------------------------------------
                        % Moving ball RF path
                        % RFuSTDirSizeLum: [nDir, nSize, nLum, rfY, rfX, nN]
                        % RFuShuffST:      [rfY, rfX, nN] — already shuffle-averaged
                        % -------------------------------------------------

                        % Read RF dimensions upfront — used by both prefDir and default paths
                        nDir_rf  = size(S_rf.RFuSTDirSizeLum, 1);
                        nSize_rf = size(S_rf.RFuSTDirSizeLum, 2);
                        nLum_rf  = size(S_rf.RFuSTDirSizeLum, 3);
                        rfY      = size(S_rf.RFuSTDirSizeLum, 4);  % typically 54
                        nN_rf    = size(S_rf.RFuSTDirSizeLum, 6);
                        nGrid     = 9;
                        blockSize = rfY / nGrid;  % typically 6

                        if params.prefDir
                            % -------------------------------------------------
                            % Preferred direction path:
                            % Select per-neuron RF slice at preferred direction,
                            % defined as the direction with the highest spike rate
                            % in NeuronVals — chosen independently of RFuSTDirSizeLum
                            % to avoid circularity.
                            %
                            % NeuronVals: [nGoodUnits, nConditions, nFeatures]
                            %   dim3=1: spike rate
                            %   dim3=6: direction in radians
                            % dim 1 of NeuronVals indexes ALL good units,
                            % so we need to map via respU_all to get responsive ones.
                            %
                            % Speed used during RF computation may differ from
                            % params.speed if only one speed was recorded —
                            % use S_rf.params.speed as ground truth.
                            % -------------------------------------------------
                            rfSpeed    = S_rf.params.speed;
                            rfField    = sprintf('Speed%d', rfSpeed);
                            NeuronResp = obj.ResponseWindow;
                            NeuronVals = NeuronResp.(rfField).NeuronVals;
                            % NeuronVals: [nGoodUnits, nConditions, nFeatures]

                            dirLabels  = NeuronVals(:,:,6);  % direction in radians per condition
                            spikeRates = NeuronVals(:,:,1);  % spike rate per condition

                            % Unique directions — same across all neurons, read from row 1
                            uDirs = unique(dirLabels(1,:));  % [1, nDir]

                            % Max spike rate per direction per good unit.
                            % Max (not mean) across conditions sharing a direction avoids
                            % dilution from other conditions (size, lum) at that direction
                            nGoodUnits     = size(NeuronVals, 1);
                            maxRespPerDir  = zeros(nGoodUnits, numel(uDirs));
                            for d = 1:numel(uDirs)
                                % Logical mask: conditions with this direction
                                dirMask = dirLabels(1,:) == uDirs(d);
                                maxRespPerDir(:,d) = max(spikeRates(:, dirMask), [], 2);
                            end

                            % Preferred direction index (1:nDir) for each good unit
                            [~, prefDirIdxAllGood] = max(maxRespPerDir, [], 2);  % [nGoodUnits, 1]

                            % Map onto responsive neurons of this stim type
                            prefDirIdxResp = prefDirIdxAllGood(respU_all{s});  % [nRespUnits, 1]

                            % Map onto shared neurons
                            [~, neuronIdx]   = ismember(sharedPhyIDs, phy_IDg(respU_all{s}));
                            prefDirIdxShared = prefDirIdxResp(neuronIdx);  % [nShared, 1]
                            nShared          = numel(neuronIdx);

                            fprintf('  [prefDir] Preferred directions for shared neurons: %s\n', ...
                                num2str(prefDirIdxShared'));

                            % Build rfRaw by selecting each neuron's preferred direction slice
                            % Result: [nSize, nLum, rfY, rfX, nShared] — no direction dim
                            rfRaw = zeros(nSize_rf, nLum_rf, rfY, rfY, nShared);
                            for u = 1:nShared
                                % Slice preferred direction for this neuron and squeeze
                                % direction singleton: [1,nSize,nLum,rfY,rfX,1] -> [nSize,nLum,rfY,rfX]
                                rfRaw(:,:,:,:,u) = squeeze( ...
                                    S_rf.RFuSTDirSizeLum(prefDirIdxShared(u),:,:,:,:,neuronIdx(u)));
                            end
                            % rfRaw: [nSize, nLum, rfY, rfX, nShared]

                            % Shuffle: select same responsive+shared neuron subset
                            rfShuff = S_rf.RFuShuffST;          % [rfY, rfX, nRespUnits] — already responsive-only
                            rfShuff = rfShuff(:,:, neuronIdx);  % [rfY, rfX, nShared] — select shared neurons

                            % neuronIdx already applied above — skip post-switch indexing
                            alreadyIndexed = true;

                        else
                            % -------------------------------------------------
                            % Default path: average RF across all directions
                            % rfRaw: [nSize, nLum, rfY, rfX, nN_rf]
                            % -------------------------------------------------
                            rfRaw = mean(S_rf.RFuSTDirSizeLum, 1);  % [1, nSize, nLum, rfY, rfX, nN]
                            rfRaw = reshape(rfRaw, [nSize_rf, nLum_rf, rfY, rfY, nN_rf]);

                            rfShuff = S_rf.RFuShuffST;  % [rfY, rfX, nN_rf]
                            nShared = nN_rf;             % neuronIdx applied after switch

                        end

                        % -------------------------------------------------
                        % Downsample rfY x rfY -> nGrid x nGrid
                        % rfRaw is now always [nSize, nLum, rfY, rfX, nShared]
                        % rfShuff is          [rfY, rfX, nShared]
                        % -------------------------------------------------
                        rfDown      = zeros(nGrid, nGrid, nShared, nSize_rf, nLum_rf);
                        rfShuffDown = zeros(nGrid, nGrid, nShared);

                        for bi = 1:nGrid
                            for bj = 1:nGrid
                                % Row and column pixel indices for this 6x6 block
                                rr = (bi-1)*blockSize + (1:blockSize);
                                cc = (bj-1)*blockSize + (1:blockSize);

                                % Average over spatial block dims 3 and 4 of rfRaw
                                % [nSize, nLum, blockSize, blockSize, nShared] -> [nSize, nLum, 1, 1, nShared]
                                block = mean(mean(rfRaw(:,:,rr,cc,:), 3), 4);
                                block = reshape(block, [nSize_rf, nLum_rf, nShared]);  % [nSize, nLum, nShared]
                                block = permute(block, [3, 1, 2]);                     % [nShared, nSize, nLum]
                                rfDown(bi,bj,:,:,:) = reshape(block, [1, 1, nShared, nSize_rf, nLum_rf]);

                                % Shuffle has no size/lum dim — average block directly
                                rfShuffDown(bi,bj,:) = mean(mean(rfShuff(rr,cc,:), 1), 2);
                            end
                        end

                        % Tile shuffle baseline across nSize and nLum to match rfDown
                        rfShuffDown = repmat( ...
                            reshape(rfShuffDown, [nGrid, nGrid, nShared, 1, 1]), ...
                            [1, 1, 1, nSize_rf, nLum_rf]);  % [nGrid, nGrid, nShared, nSize, nLum]

                        gridSpikeRateSelected = rfDown;       % [nGrid, nGrid, nShared, nSize, nLum]
                        gridShuffMean         = rfShuffDown;  % [nGrid, nGrid, nShared, nSize, nLum]

                    case "rectGrid"
                        % -------------------------------------------------
                        % rectGrid RF path
                        % RFu:          [2(onOff), nLums, nSize, screenRed, screenRed, nN]
                        % RFuShuffMean: same dimensions
                        % IMPORTANT: dim2=nLums, dim3=nSize (not the other way around)
                        % -------------------------------------------------

                        % Select on or off response (dim 1); keep remaining dims explicit
                        rfFull      = S_rf.RFu(         params.onOff, :, :, :, :, :);  % [1, nLums, nSize, screenRed, screenRed, nN]
                        rfShuffFull = S_rf.RFuShuffMean(params.onOff, :, :, :, :, :);  % same

                        nLums_rf  = size(rfFull, 2);
                        nSize_rf  = size(rfFull, 3);
                        screenRed = size(rfFull, 4);  % reduced screen resolution
                        nN_rf     = size(rfFull, 6);

                        % Collapse the leading singleton onOff dim
                        rfRaw      = reshape(rfFull,      [nLums_rf, nSize_rf, screenRed, screenRed, nN_rf]);
                        rfShuffRaw = reshape(rfShuffFull, [nLums_rf, nSize_rf, screenRed, screenRed, nN_rf]);
                        % rfRaw: [nLums, nSize, screenRed, screenRed, nN]

                        nGrid     = 9;
                        blockSize = screenRed / nGrid;

                        rfDown      = zeros(nGrid, nGrid, nN_rf, nSize_rf, nLums_rf);
                        rfShuffDown = zeros(nGrid, nGrid, nN_rf, nSize_rf, nLums_rf);

                        for bi = 1:nGrid
                            for bj = 1:nGrid
                                % Row and column pixel indices for this spatial block
                                rr = (bi-1)*blockSize + (1:blockSize);
                                cc = (bj-1)*blockSize + (1:blockSize);

                                % Average over spatial dims 3 and 4
                                % [nLums, nSize, blockSize, blockSize, nN] -> [nLums, nSize, 1, 1, nN]
                                block      = mean(mean(rfRaw(     :,:,rr,cc,:), 3), 4);
                                blockShuff = mean(mean(rfShuffRaw(:,:,rr,cc,:), 3), 4);

                                % Reshape and permute to [nN, nSize, nLums]
                                % Note: after reshape input order is [nLums, nSize, nN]
                                block      = reshape(block,      [nLums_rf, nSize_rf, nN_rf]);
                                blockShuff = reshape(blockShuff, [nLums_rf, nSize_rf, nN_rf]);
                                block      = permute(block,      [3, 2, 1]);  % [nN, nSize, nLums]
                                blockShuff = permute(blockShuff, [3, 2, 1]);

                                rfDown(bi,bj,:,:,:)      = reshape(block,      [1, 1, nN_rf, nSize_rf, nLums_rf]);
                                rfShuffDown(bi,bj,:,:,:) = reshape(blockShuff, [1, 1, nN_rf, nSize_rf, nLums_rf]);
                            end
                        end

                        gridSpikeRateSelected = rfDown;       % [nGrid, nGrid, nN, nSize, nLum]
                        gridShuffMean         = rfShuffDown;  % [nGrid, nGrid, nN, nSize, nLum]

                end % switch stimType (RF path)

                % Apply shared neuron indexing after the switch block.
                % Skipped for linearlyMovingBall+prefDir since neuronIdx
                % was already applied per-neuron inside that branch.
                if ~alreadyIndexed
                    [~, neuronIdx] = ismember(sharedPhyIDs, phy_IDg(respU_all{s}));
                    gridSpikeRateSelected = gridSpikeRateSelected(:,:,neuronIdx,:,:);
                    gridShuffMean         = gridShuffMean(:,:,neuronIdx,:,:);
                end

            else
                % ----------------------------------------------------------
                % Standard path: use gridSpikeRate / gridSpikeRateShuff
                % Note: linearlyMovingBall may have structural zeros due to
                % Y-offset bug in CalculateReceptiveFields — use useRF=true
                % ----------------------------------------------------------
                if stimType == "linearlyMovingBall"
                    warning(['gridSpikeRate for linearlyMovingBall may contain structural zeros ' ...
                             'due to Y-offset bug in CalculateReceptiveFields. ' ...
                             'Consider using useRF=true.']);
                end

                gridSpikeRate      = S_rf.gridSpikeRate;
                gridSpikeRateShuff = S_rf.gridSpikeRateShuff;

                switch stimType
                    case "rectGrid"
                        % gridSpikeRate:      [nGrid, nGrid, nN, 2, nSize, nLum]
                        % gridSpikeRateShuff: [nGrid, nGrid, nN, nShuffle, 2, nSize, nLum]
                        gridSpikeRateSelected = gridSpikeRate(:,:,:,params.onOff,:,:);
                        gridShuffSelected     = gridSpikeRateShuff(:,:,:,:,params.onOff,:,:);

                        % Remove onOff singleton: [9,9,nN,1,nSize,nLum] -> [9,9,nN,nSize,nLum]
                        gridSpikeRateSelected = reshape(gridSpikeRateSelected, ...
                            [size(gridSpikeRateSelected,1), size(gridSpikeRateSelected,2), ...
                             size(gridSpikeRateSelected,3), size(gridSpikeRateSelected,5), ...
                             size(gridSpikeRateSelected,6)]);

                        % Remove onOff singleton: [9,9,nN,nShuffle,1,nSize,nLum] -> [9,9,nN,nShuffle,nSize,nLum]
                        gridShuffSelected = reshape(gridShuffSelected, ...
                            [size(gridShuffSelected,1), size(gridShuffSelected,2), ...
                             size(gridShuffSelected,3), size(gridShuffSelected,4), ...
                             size(gridShuffSelected,6), size(gridShuffSelected,7)]);

                    case "linearlyMovingBall"
                        gridSpikeRateSelected = gridSpikeRate;       % [nGrid, nGrid, nN, nSize, nLum]
                        gridShuffSelected     = gridSpikeRateShuff;  % [nGrid, nGrid, nN, nShuffle, nSize, nLum]
                end

                % Map sharedPhyIDs onto indices of this stim's responsive neurons
                [~, neuronIdx] = ismember(sharedPhyIDs, phy_IDg(respU_all{s}));

                gridSpikeRateSelected = gridSpikeRateSelected(:,:,neuronIdx,:,:);
                gridShuffSelected     = gridShuffSelected(:,:,neuronIdx,:,:,:);

                % Average shuffle dimension (dim 4) to get baseline map
                gridShuffMean = mean(gridShuffSelected, 4);  % [nGrid, nGrid, nN, 1, nSize, nLum]

            end % useRF / standard path

            % ----------------------------------------------------------
            % Get dimensions and reshape to canonical [nGrid,nGrid,nN,nSize,nLum]
            % ----------------------------------------------------------
            nN    = size(gridSpikeRateSelected, 3);
            nSize = size(gridSpikeRateSelected, 4);
            nLum  = size(gridSpikeRateSelected, 5);
            nGrid = size(gridSpikeRateSelected, 1);

            fprintf('gridSpikeRateSelected size: %s\n', num2str(size(gridSpikeRateSelected)));

            gridSpikeRateSelected = reshape(gridSpikeRateSelected, [nGrid, nGrid, nN, nSize, nLum]);
            gridShuffMean         = reshape(gridShuffMean,         [nGrid, nGrid, nN, nSize, nLum]);

            nCells  = nGrid * nGrid;
            maxDist = sqrt(2) * (nGrid - 1);  % maximum possible distance between two grid cells

            % ----------------------------------------------------------
            % Compute spatial tuning indices per neuron, size, and lum
            % ----------------------------------------------------------
            for si = 1:nSize
                for li = 1:nLum

                    % Flatten spatial dims: [nCells, nN]
                    rateFlat      = reshape(gridSpikeRateSelected(:,:,:,si,li), [nCells, nN]);
                    rateFlatShuff = reshape(gridShuffMean(:,:,:,si,li),          [nCells, nN]);

                    L_amplitude_diff  = zeros(nN, 1);
                    L_amplitude_ratio = zeros(nN, 1);
                    L_geometric       = zeros(nN, 1);
                    L_combined        = zeros(nN, 1);

                    for u = 1:nN

                        rateVec      = rateFlat(:, u);       % spike rate at each grid cell
                        rateVecShuff = rateFlatShuff(:, u);  % shuffle baseline at each grid cell

                        % Threshold for top-percent most active grid cells
                        threshold      = prctile(rateVec,      100 - params.topPercent);
                        thresholdShuff = prctile(rateVecShuff, 100 - params.topPercent);

                        % Indices of top and rest cells for real and shuffle maps
                        topIdx       = find(rateVec      >= threshold);
                        topIdxShuff  = find(rateVecShuff >= thresholdShuff);
                        restIdx      = setdiff(1:nCells, topIdx);
                        restIdxShuff = setdiff(1:nCells, topIdxShuff);

                        % Mean rates in top and rest regions for real and shuffle.
                        % Guard against empty restIdx: occurs when topPercent is large
                        % enough that all cells exceed the threshold (all tied at zero).
                        % In that case there is no spatial contrast — set meanRest = 0.
                        meanTop  = mean(rateVec(topIdx));
                        meanAll  = mean(rateVec);
                        if isempty(restIdx)
                            meanRest = 0;
                        else
                            meanRest = mean(rateVec(restIdx));
                        end

                        meanTopShuff = mean(rateVecShuff(topIdxShuff));
                        meanAllShuff = mean(rateVecShuff);
                        if isempty(restIdxShuff)
                            meanRestShuff = 0;
                        else
                            meanRestShuff = mean(rateVecShuff(restIdxShuff));
                        end

                        % Guard against division by zero in normalisation
                        if meanAll      == 0, meanAll      = eps; end
                        if meanAllShuff == 0, meanAllShuff = eps; end

                        % L_amplitude_diff: normalised contrast, shuffle-subtracted
                        L_amplitude_diff(u) = ...
                            (meanTop - meanRest) / meanAll - ...
                            (meanTopShuff - meanRestShuff) / meanAllShuff;

                        % L_amplitude_ratio: real contrast divided by shuffle contrast
                        shuffleNorm = (meanTopShuff - meanRestShuff) / meanAllShuff;
                        if shuffleNorm == 0, shuffleNorm = eps; end
                        L_amplitude_ratio(u) = ((meanTop - meanRest) / meanAll) / shuffleNorm;

                        % L_geometric: clustering of top cells (low spread = high tuning)
                        % Convert linear indices to [row, col] grid coordinates
                        [rowIdx,      colIdx]      = ind2sub([nGrid nGrid], topIdx);
                        [rowIdxShuff, colIdxShuff] = ind2sub([nGrid nGrid], topIdxShuff);

                        % Mean pairwise distance among top cells, normalised by max possible distance
                        if size(rowIdx, 1) > 1
                            D = mean(pdist([rowIdx, colIdx], 'euclidean')) / maxDist;
                        else
                            D = 0;  % single top cell: perfectly localised by definition
                        end
                        if size(rowIdxShuff, 1) > 1
                            DShuff = mean(pdist([rowIdxShuff, colIdxShuff], 'euclidean')) / maxDist;
                        else
                            DShuff = 0;
                        end

                        % Shuffle-corrected geometric index: positive = more clustered than chance
                        L_geometric(u) = (1 - D) - (1 - DShuff);

                        % L_combined: product of amplitude and geometric indices
                        L_combined(u)  = L_amplitude_diff(u) * L_geometric(u);

                    end % neuron loop

                    % Check for NaN indices and report their source for debugging
                    nanMask = isnan(L_amplitude_diff) | isnan(L_amplitude_ratio) | ...
                              isnan(L_geometric)      | isnan(L_combined);
                    if any(nanMask)
                        fprintf('  WARNING: %d/%d neurons have NaN index values [stim=%s, si=%d, li=%d]\n', ...
                            sum(nanMask), nN, char(stimType), si, li);
                    end

                    % Build one table row per neuron for this condition
                    rows                   = table();
                    rows.L_amplitude_diff  = L_amplitude_diff;
                    rows.L_amplitude_ratio = L_amplitude_ratio;
                    rows.L_geometric       = L_geometric;
                    rows.L_combined        = L_combined;
                    rows.stimulus          = categorical(repmat({char(stimType)}, nN, 1));
                    rows.experimentNum     = categorical(repmat(ex,               nN, 1));
                    rows.animal            = categorical(repmat({animalName},     nN, 1));
                    rows.NeurID            = (1:nN)';
                    % Store actual phy cluster ID for each neuron.
                    % After neuronIdx selection, neuron u in dim 3 corresponds to sharedPhyIDs(u).
                    rows.phyID             = sharedPhyIDs(:);
                    rows.onOff             = repmat(params.onOff, nN, 1);  % meaningful for rectGrid; stored for consistency
                    rows.sizeIdx           = repmat(si, nN, 1);
                    rows.lumIdx            = repmat(li, nN, 1);

                    tbl = [tbl; rows];  %#ok<AGROW>

                end % lum loop
            end % size loop

            fprintf('  [%s] Indices computed. %d neurons.\n', stimType, nN);

        end % stim loop
    end % exp loop

    % Remove unused categorical levels introduced by partial data
    tbl.stimulus      = removecats(tbl.stimulus);
    tbl.animal        = removecats(tbl.animal);
    tbl.experimentNum = removecats(tbl.experimentNum);

    % Cache results to disk
    S.expList = exList;
    S.tbl     = tbl;
    S.params  = params;
    save([saveDir nameOfFile], '-struct', 'S');
    fprintf('\nSaved to:\n  %s\n', [saveDir nameOfFile]);

end % compute block

results.tbl = tbl;

% =========================================================================
% TOP UNIT TABLES
% Top 20% of neurons globally by params.indexType, for each stim type
% separately. Uses the same condition filter as the plot (onOff/sizeIdx/lumIdx).
% =========================================================================

% Filter to the requested condition — same as plot filter
idxCond = tbl.onOff   == params.onOff   & ...
          tbl.sizeIdx == params.sizeIdx & ...
          tbl.lumIdx  == params.lumIdx;

tblCond       = tbl(idxCond, :);
tblCond.value = tblCond.(params.indexType);  % column to rank on

% Identify the stim label for SB (rectGrid) and MB (linearlyMovingBall)
sbLabel = 'rectGrid';
mbLabel = 'linearlyMovingBall';

% Build one top-unit table per stim type
for tt = 1:2

    if tt == 1
        stimLabel = sbLabel;
        outField  = 'topUnitsSB';
    else
        stimLabel = mbLabel;
        outField  = 'topUnitsMB';
    end

    % Check this stim type was actually computed
    if ~any(tblCond.stimulus == stimLabel)
        fprintf('  No data for %s — skipping top unit table.\n', stimLabel);
        results.(outField) = table();
        continue
    end

    % Subset to this stim type
    tblStim = tblCond(tblCond.stimulus == stimLabel, :);

    % Global threshold: top 20% across all animals and insertions
    globalThreshold = prctile(tblStim.value, 80);

    % Select top units and sort descending by index value
    topMask = tblStim.value >= globalThreshold;
    tblTop  = sortrows(tblStim(topMask, :), 'value', 'descend');

    % Build clean output table with one row per top unit
    outTbl              = table();
    outTbl.animal       = tblTop.animal;
    outTbl.experimentNum = tblTop.experimentNum;
    outTbl.phyID        = tblTop.phyID;      % phy cluster ID (Kilosort/phy)
    outTbl.indexValue   = tblTop.value;      % spatial tuning index value

    fprintf('  [%s] %d top units (top 20%%, threshold=%.4f).\n', ...
        stimLabel, height(outTbl), globalThreshold);

    results.(outField) = outTbl;

end

% =========================================================================
% PLOT
% =========================================================================
if params.plot

    % Filter table to the requested on/off, size, and lum condition
    idx = tbl.onOff   == params.onOff   & ...
          tbl.sizeIdx == params.sizeIdx & ...
          tbl.lumIdx  == params.lumIdx;

    tblPlot       = tbl(idx, :);
    tblPlot.value = tblPlot.(params.indexType);  % select which index to plot

    % ----------------------------------------------------------
    % Compute hierarchical bootstrap p-value for the comparison pair
    % ----------------------------------------------------------
    pairs = {char(params.stimTypes(1)), char(params.stimTypes(2))};  % 1x2 cell

    ps = zeros(size(pairs, 1), 1);
    j  = 1;

    for i = 1:size(pairs, 1)
        diffs   = [];
        insers  = [];
        animals = [];

        % Compute per-neuron differences within each insertion
        for ins = unique(tblPlot.experimentNum)'
            idx1 = tblPlot.experimentNum == categorical(ins) & tblPlot.stimulus == pairs{i,1};
            idx2 = tblPlot.experimentNum == categorical(ins) & tblPlot.stimulus == pairs{i,2};

            V1 = tblPlot.value(idx1);
            V2 = tblPlot.value(idx2);

            if isempty(V1) || isempty(V2)
                continue
            end

            animal  = unique(tblPlot.animal(idx1));
            diffs   = [diffs;   V1 - V2];                                    %#ok<AGROW>
            insers  = [insers;  double(repmat(ins,    size(V1,1), 1))];      %#ok<AGROW>
            animals = [animals; double(repmat(animal, size(V1,1), 1))];      %#ok<AGROW>
        end

        if isempty(diffs)
            ps(j) = NaN;
        else
            % Hierarchical bootstrap: respects nested structure
            % (neurons within insertions within animals)
            bootDiff = hierBoot(diffs, params.nBoot, insers, animals);
            ps(j)    = mean(bootDiff <= 0);  % one-tailed p: P(stim1 <= stim2)
        end
        j = j + 1;
    end

    % ----------------------------------------------------------
    % Plot swarm with bootstrap confidence intervals
    % ----------------------------------------------------------
    V1max = max(tblPlot.value, [], 'omitnan');  % data max for y-axis scaling

    fprintf('Length of ps: %d\n', numel(ps));
    fprintf('Size of pairs: %s\n', num2str(size(pairs)));

    tblPlot.insertion = tblPlot.experimentNum;  % rename for plotting compatibility

    [fig, ~] = plotSwarmBootstrapWithComparisons(tblPlot, pairs, ps, {'value'}, ...
        yLegend     = params.yLegend, ...
        yMaxVis     = max(params.yMaxVis, V1max), ...
        diff        = false, ...
        Alpha       = params.Alpha, ...
        plotMeanSem = true);

    title(sprintf('%s — %s  (onOff=%d, size=%d, lum=%d, RF=%d, prefDir=%d)', ...
        params.indexType, strjoin(params.stimTypes, '/'), ...
        params.onOff, params.sizeIdx, params.lumIdx, params.useRF, params.prefDir), ...
        'FontSize', 9);

    % Save publication-quality figure if requested
    if params.PaperFig
        vs_first.printFig(fig, sprintf('SpatialTuningIndex-%s-%s-RF%d-prefDir%d', ...
            params.indexType, strjoin(params.stimTypes, '-'), params.useRF, params.prefDir), ...
            PaperFig = params.PaperFig);
    end

    results.fig = fig;
    results.ps  = ps;

end

end