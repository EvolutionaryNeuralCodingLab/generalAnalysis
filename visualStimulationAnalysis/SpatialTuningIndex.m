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
    params.prefDir    logical       = true   % If true (requires useRF=true): use each neuron's preferred
                                             % direction RF for linearlyMovingBall instead of averaging
                                             % across all directions. Preferred direction is defined as the
                                             % direction with the highest spike rate in NeuronVals.
                                             % Avoids deflating the spatial tuning index by averaging over
                                             % non-preferred directions.
    params.allResponsive logical    = false  % If true: compute index for ALL neurons responsive to each
                                             % stim type independently — no intersection across stim types.
                                             % Each swarm column will have a different number of neurons.
                                             % Difference plot is not available in this mode.
                                             % P-value is computed via two-sample hierarchical bootstrap.
    params.unionResponsive logical = false  % If true: compute index for all neurons responsive
                                            % to EITHER stim type (union). Same neuron set used
                                            % for both stim types, so paired diff is still valid.
                                            % P-value uses paired hierBoot on differences.
    params.plotRFs         logical = false  % If true: generate multi-page PDF
                                              % showing each neuron's full-resolution
                                              % receptive field, sorted by tuning index
                                              % (descending). Requires prefDir=true for
                                              % linearlyMovingBall.
    params.subtractShuffle logical = true   % If true: subtract the shuffle-mean baseline
                                           % from each RF before plotting, so the colour
                                           % scale reflects signal above noise.
                                           % If false: plot the raw (unsubtracted) RF.
    params.plotRFunion logical = false     % If true: RF pages show all neurons responsive
                                           % to EITHER stim type (union). Two PDFs are
                                           % generated with the same neuron set, one per
                                           % stim type. Each PDF is sorted by that stim
                                           % type's tuning index; neurons not responsive
                                           % to the plotted stim type are annotated "n/r"
                                           % and sink to the bottom.
end

% -------------------------------------------------------------------------
% Parameter guards
% -------------------------------------------------------------------------
% prefDir requires the RF path
if params.prefDir && ~params.useRF
    error('prefDir=true requires useRF=true. The preferred direction RF is only available in the RF path.');
end

% allResponsive is incompatible with diff plotting — neurons are unpaired
if params.allResponsive
    fprintf('allResponsive=true: each stim type will include all its own responsive neurons.\n');
    fprintf('  Difference plot is not available in this mode.\n');
    fprintf('  P-value computed via two-sample hierarchical bootstrap.\n');
end

% unionResponsive and allResponsive are mutually exclusive
if params.unionResponsive && params.allResponsive
    error('unionResponsive and allResponsive cannot both be true — choose one neuron selection mode.');
end

if params.unionResponsive
    fprintf('unionResponsive=true: using all neurons responsive to either stim type (union).\n');
    fprintf('  Paired difference plot is available since the same neuron set is used for both stim types.\n');
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

% Build filename encoding all parameter modes that affect computation
% so different parameter combinations never collide on disk
stimLabel    = strjoin(params.stimTypes, '-');
rfLabel      = '';   if params.useRF,        rfLabel      = '_RF';          end
prefLabel    = '';   if params.prefDir,       prefLabel    = '_prefDir';     end
allRespLabel = '';   if params.allResponsive, allRespLabel = '_allResp';     end
unionLabel   = '';   if params.unionResponsive, unionLabel = '_union'; end
nameOfFile   = sprintf('\\Ex_%d-%d_SpatialTuningIndex_%s%s%s%s%s.mat', ...
    exList(1), exList(end), stimLabel, rfLabel, prefLabel, allRespLabel, unionLabel);

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
        % Determine which neurons to include per stim type:
        %   allResponsive=false: intersection across all stim types (paired)
        %   allResponsive=true:  each stim type uses its own responsive set
        % ----------------------------------------------------------
        if ~params.allResponsive && ~params.unionResponsive

            % Paired mode: intersection of responsive neurons across stim types
            sharedPhyIDs = respPhyIDs_all{1};
            for s = 2:nStim
                sharedPhyIDs = intersect(sharedPhyIDs, respPhyIDs_all{s});
            end

            if isempty(sharedPhyIDs)
                fprintf('  No neurons responsive to all stim types in exp %d — skipping.\n', ex);
                continue
            end
            fprintf('  %d neuron(s) responsive to all stim types in exp %d.\n', numel(sharedPhyIDs), ex);

            % Same set for all stim types
            sharedPhyIDs_perStim = repmat({sharedPhyIDs}, 1, nStim);

        elseif params.unionResponsive

            % Union mode: neurons responsive to ANY stim type
            % Same union set applied to all stim types — neurons not
            % responsive to a given stim will still have their RF/grid
            % map computed, which may be weak but is not excluded.
            % Paired diff is still valid since every neuron has an index
            % for both stim types.
            unionPhyIDs = respPhyIDs_all{1};
            for s = 2:nStim
                unionPhyIDs = union(unionPhyIDs, respPhyIDs_all{s});
            end

            if isempty(unionPhyIDs)
                fprintf('  No responsive neurons for any stim type in exp %d — skipping.\n', ex);
                continue
            end
            fprintf('  %d neuron(s) in union (responsive to at least one stim type) in exp %d.\n', ...
                numel(unionPhyIDs), ex);

            % Same union set for all stim types
            sharedPhyIDs_perStim = repmat({unionPhyIDs}, 1, nStim);

        else

            % Unpaired mode: each stim type uses its own full responsive set
            anyStimsHaveNeurons = any(cellfun(@(x) ~isempty(x), respPhyIDs_all));
            if ~anyStimsHaveNeurons
                fprintf('  No responsive neurons for any stim type in exp %d — skipping.\n', ex);
                continue
            end
            sharedPhyIDs_perStim = respPhyIDs_all;
            for s = 1:nStim
                fprintf('  [%s] %d neuron(s) (all responsive, unpaired).\n', ...
                    params.stimTypes(s), numel(sharedPhyIDs_perStim{s}));
            end

        end

        % ----------------------------------------------------------
        % Loop over stim types and compute spatial tuning index
        % ----------------------------------------------------------
        for s = 1:nStim

            stimType     = params.stimTypes(s);
            sharedPhyIDs = sharedPhyIDs_perStim{s};  % neurons for THIS stim type

            if isempty(sharedPhyIDs)
                fprintf('  [%s] No neurons — skipping.\n', char(stimType));
                continue
            end

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
                        nDir_rf  = size(S_rf.RFuSTDirSizeLum, 1); %#ok<NASGU>
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
                            % dim 1 of NeuronVals indexes ALL good units.
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
                            nGoodUnits    = size(NeuronVals, 1);
                            maxRespPerDir = zeros(nGoodUnits, numel(uDirs));
                            for d = 1:numel(uDirs)
                                dirMask = dirLabels(1,:) == uDirs(d);
                                maxRespPerDir(:,d) = max(spikeRates(:, dirMask), [], 2);
                            end

                            % Preferred direction index (1:nDir) for each good unit
                            [~, prefDirIdxAllGood] = max(maxRespPerDir, [], 2);  % [nGoodUnits, 1]

    

                            % Map sharedPhyIDs (which are a subset of respU_all{s}'s phy IDs)
                            % onto indices within the responsive set
                            [~, neuronIdx] = ismember(sharedPhyIDs, phy_IDg);
                            prefDirIdxShared = prefDirIdxAllGood(neuronIdx);  % [nShared, 1]
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

                            % Shuffle: RFuShuffST is already responsive-only
                            % just select shared neurons
                            rfShuff = S_rf.RFuShuffST;          % [rfY, rfX, nRespUnits]
                            rfShuff = rfShuff(:,:, neuronIdx);  % [rfY, rfX, nShared]

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
                        rfFull      = S_rf.RFu(         params.onOff, :, :, :, :, :);
                        rfShuffFull = S_rf.RFuShuffMean(params.onOff, :, :, :, :, :);

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
                                rr = (bi-1)*blockSize + (1:blockSize);
                                cc = (bj-1)*blockSize + (1:blockSize);

                                % [nLums, nSize, blockSize, blockSize, nN] -> [nLums, nSize, 1, 1, nN]
                                block      = mean(mean(rfRaw(     :,:,rr,cc,:), 3), 4);
                                blockShuff = mean(mean(rfShuffRaw(:,:,rr,cc,:), 3), 4);

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
                    [~, neuronIdx] = ismember(sharedPhyIDs, phy_IDg);
                    gridSpikeRateSelected = gridSpikeRateSelected(:,:,neuronIdx,:,:);
                    gridShuffMean         = gridShuffMean(:,:,neuronIdx,:,:);
                end

            else
                % ----------------------------------------------------------
                % Standard path: use gridSpikeRate / gridSpikeRateShuff
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
                        gridSpikeRateSelected = gridSpikeRate(:,:,:,params.onOff,:,:);
                        gridShuffSelected     = gridSpikeRateShuff(:,:,:,:,params.onOff,:,:);

                        gridSpikeRateSelected = reshape(gridSpikeRateSelected, ...
                            [size(gridSpikeRateSelected,1), size(gridSpikeRateSelected,2), ...
                             size(gridSpikeRateSelected,3), size(gridSpikeRateSelected,5), ...
                             size(gridSpikeRateSelected,6)]);

                        gridShuffSelected = reshape(gridShuffSelected, ...
                            [size(gridShuffSelected,1), size(gridShuffSelected,2), ...
                             size(gridShuffSelected,3), size(gridShuffSelected,4), ...
                             size(gridShuffSelected,6), size(gridShuffSelected,7)]);

                    case "linearlyMovingBall"
                        gridSpikeRateSelected = gridSpikeRate;
                        gridShuffSelected     = gridSpikeRateShuff;
                end

                [~, neuronIdx] = ismember(sharedPhyIDs, phy_IDg(respU_all{s}));

                gridSpikeRateSelected = gridSpikeRateSelected(:,:,neuronIdx,:,:);
                gridShuffSelected     = gridShuffSelected(:,:,neuronIdx,:,:,:);

                % Average shuffle dimension (dim 4) to get baseline map
                gridShuffMean = mean(gridShuffSelected, 4);

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

            %Capture preferred-direction info for the current stim×neuron set.
            % These vectors will be written into the table inside the si/li loop.
            if stimType == "linearlyMovingBall" && params.prefDir && params.useRF
                % prefDirIdxShared: index (1:nDir) of each neuron's preferred direction
                storePrefDirIdx = prefDirIdxShared(:);          % [nShared, 1]
                % Convert preferred-direction index to degrees via the sorted unique directions
                storePrefDirDeg = rad2deg(uDirs(prefDirIdxShared(:)))';  % [nShared, 1]
            else
                % Non-MB or non-prefDir: fill with NaN so the table column exists
                storePrefDirIdx = nan(size(gridSpikeRateSelected, 3), 1);
                storePrefDirDeg = nan(size(gridSpikeRateSelected, 3), 1);
            end

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

                        rateVec      = rateFlat(:, u);
                        rateVecShuff = rateFlatShuff(:, u);

                        % Threshold for top-percent most active grid cells
                        threshold      = prctile(rateVec,      100 - params.topPercent);
                        thresholdShuff = prctile(rateVecShuff, 100 - params.topPercent);

                        topIdx       = find(rateVec      >= threshold);
                        topIdxShuff  = find(rateVecShuff >= thresholdShuff);
                        restIdx      = setdiff(1:nCells, topIdx);
                        restIdxShuff = setdiff(1:nCells, topIdxShuff);

                        % Mean rates in top and rest regions.
                        % Guard against empty restIdx (all cells above threshold
                        % when topPercent is large) — set meanRest=0 in that case
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
                        [rowIdx,      colIdx]      = ind2sub([nGrid nGrid], topIdx);
                        [rowIdxShuff, colIdxShuff] = ind2sub([nGrid nGrid], topIdxShuff); %#ok<ASGLU>

                        if size(rowIdx, 1) > 1
                            D = mean(pdist([rowIdx, colIdx], 'euclidean')) / maxDist;
                        else
                            D = 0;
                        end
                        if size(rowIdxShuff, 1) > 1
                            DShuff = mean(pdist([rowIdxShuff, colIdxShuff], 'euclidean')) / maxDist;
                        else
                            DShuff = 0;
                        end

                        % Shuffle-corrected geometric index
                        L_geometric(u) = (1 - D) - (1 - DShuff);

                        % L_combined: product of amplitude and geometric indices
                        L_combined(u)  = L_amplitude_diff(u) * L_geometric(u);

                    end % neuron loop

                    % Check for NaN indices and report for debugging
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
                    % phyID: sharedPhyIDs is stim-specific when allResponsive=true
                    rows.phyID             = sharedPhyIDs(:);
                    rows.onOff             = repmat(params.onOff, nN, 1);
                    rows.sizeIdx           = repmat(si, nN, 1);
                    rows.lumIdx            = repmat(li, nN, 1);
                    
                    % Preferred-direction index into the direction dimension of
                    % RFuSTDirSizeLum.  NaN for rectGrid or when prefDir=false.
                    rows.prefDirIdx = storePrefDirIdx;
                    % Preferred direction in degrees (0–360). NaN when not applicable.
                    rows.prefDirDeg = storePrefDirDeg;

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

idxCond       = tbl.onOff   == params.onOff   & ...
                tbl.sizeIdx == params.sizeIdx & ...
                tbl.lumIdx  == params.lumIdx;
tblCond       = tbl(idxCond, :);
tblCond.value = tblCond.(params.indexType);

sbLabel = 'rectGrid';
mbLabel = 'linearlyMovingBall';

for tt = 1:2
    if tt == 1
        stimLabel = sbLabel;  outField = 'topUnitsSB';
    else
        stimLabel = mbLabel;  outField = 'topUnitsMB';
    end

    if ~any(tblCond.stimulus == stimLabel)
        fprintf('  No data for %s — skipping top unit table.\n', stimLabel);
        results.(outField) = table();
        continue
    end

    tblStim         = tblCond(tblCond.stimulus == stimLabel, :);
    globalThreshold = prctile(tblStim.value, 80);
    topMask         = tblStim.value >= globalThreshold;
    tblTop          = sortrows(tblStim(topMask, :), 'value', 'descend');

    outTbl               = table();
    outTbl.animal        = tblTop.animal;
    outTbl.experimentNum = tblTop.experimentNum;
    outTbl.phyID         = tblTop.phyID;
    outTbl.indexValue    = tblTop.value;

    fprintf('  [%s] %d top units (top 20%%, threshold=%.4f).\n', ...
        stimLabel, height(outTbl), globalThreshold);
    results.(outField) = outTbl;
end

% =========================================================================
% PLOT
% =========================================================================
if params.plot

    % Filter table to the requested on/off, size, and lum condition
    idx     = tbl.onOff   == params.onOff   & ...
              tbl.sizeIdx == params.sizeIdx & ...
              tbl.lumIdx  == params.lumIdx;
    tblPlot = tbl(idx, :);
    tblPlot.value     = tblPlot.(params.indexType);
    tblPlot.insertion = tblPlot.experimentNum;  % rename for plotting compatibility

    pairs = {char(params.stimTypes(1)), char(params.stimTypes(2))};

    % ----------------------------------------------------------
    % Compute p-values:
    %   allResponsive=false: paired hierBoot on per-neuron differences
    %   allResponsive=true:  two-sample hierBoot on each group separately
    % ----------------------------------------------------------
    ps = zeros(size(pairs, 1), 1);
    j  = 1;

    for i = 1:size(pairs, 1)

        if ~params.allResponsive
            % ---------------------------------------------------------
            % Paired mode: per-neuron differences within each insertion,
            % then single hierBoot on the difference vector
            % ---------------------------------------------------------
            diffs   = [];
            insers  = [];
            animals = [];

            for ins = unique(tblPlot.insertion)'
                idx1 = tblPlot.insertion == categorical(ins) & tblPlot.stimulus == pairs{i,1};
                idx2 = tblPlot.insertion == categorical(ins) & tblPlot.stimulus == pairs{i,2};

                V1 = tblPlot.value(idx1);
                V2 = tblPlot.value(idx2);
                if isempty(V1) || isempty(V2), continue; end

                animal  = unique(tblPlot.animal(idx1));
                diffs   = [diffs;   V1 - V2];                                     %#ok<AGROW>
                insers  = [insers;  double(repmat(ins,    size(V1,1), 1))];       %#ok<AGROW>
                animals = [animals; double(repmat(animal, size(V1,1), 1))];       %#ok<AGROW>
            end

            if isempty(diffs)
                ps(j) = NaN;
            else
                % hierBoot on the paired difference, respecting nesting
                bootDiff = hierBoot(diffs, params.nBoot, insers, animals);
                ps(j)    = mean(bootDiff <= 0);  % P(stim1 <= stim2)
            end

        else
            % ---------------------------------------------------------
            % Unpaired (two-sample) mode:
            % Bootstrap each group separately, then compare distributions.
            % Nesting: neurons within insertions within animals.
            % ---------------------------------------------------------
            mask1 = tblPlot.stimulus == pairs{i,1};
            mask2 = tblPlot.stimulus == pairs{i,2};

            V1      = tblPlot.value(mask1);
            insers1 = double(tblPlot.insertion(mask1));
            anim1   = double(tblPlot.animal(mask1));

            V2      = tblPlot.value(mask2);
            insers2 = double(tblPlot.insertion(mask2));
            anim2   = double(tblPlot.animal(mask2));

            % Remove NaNs from each group independently
            valid1 = ~isnan(V1);
            valid2 = ~isnan(V2);

            if sum(valid1) < 3 || sum(valid2) < 3
                ps(j) = NaN;
                fprintf('  Not enough valid values for two-sample bootstrap (pair %d).\n', i);
            else
                % hierBoot on each group separately — same nesting structure
                % as the paired case but applied independently per group
                BootFirst = hierBoot(V1(valid1), params.nBoot, insers1(valid1), anim1(valid1));
                BootSec   = hierBoot(V2(valid2), params.nBoot, insers2(valid2), anim2(valid2));

                % One-tailed p: P(group2 >= group1), i.e. P(stim2 >= stim1)
                ps(j) = mean(BootSec >= BootFirst);
            end
        end

        j = j + 1;
    end

    % ----------------------------------------------------------
    % Plot swarm with bootstrap confidence intervals
    % ----------------------------------------------------------
    V1max = max(tblPlot.value, [], 'omitnan');

    fprintf('Length of ps: %d\n', numel(ps));
    fprintf('Size of pairs: %s\n', num2str(size(pairs)));

    % In allResponsive mode: suppress diff plot and inter-neuron lines
    % (neurons are unpaired so neither makes biological sense)
    useDiff  = false;  % diff is never shown directly here — handled inside plotSwarm
    useBoth =  ~params.allResponsive; 

    if isequal(pairs{1},'linearlyMovingBall'), pairs{1} = 'MB'; end

    if isequal(pairs{2},'rectGrid'), pairs{2} = 'SB'; end

    tblPlot.stimulus(tblPlot.stimulus == "linearlyMovingBall") = "MB";
    tblPlot.stimulus(tblPlot.stimulus == "rectGrid")           = "SB";

    [fig, ~] = plotSwarmBootstrapWithComparisons(tblPlot, pairs, ps, {'value'}, ...
        yLegend      = params.yLegend, ...
        yMaxVis      = max(params.yMaxVis, V1max), ...
        diff         = useDiff, ...
        Alpha        = params.Alpha, ...
        plotMeanSem  = true, ...
        drawLines    = false, ...
        showBothAndDiff = useBoth);
    
    % title(sprintf('%s — %s  (onOff=%d, size=%d, lum=%d, RF=%d, prefDir=%d, allResp=%d, union=%d)', ...
    %     params.indexType, strjoin(params.stimTypes, '/'), ...
    %     params.onOff, params.sizeIdx, params.lumIdx, ...
    %     params.useRF, params.prefDir, params.allResponsive, params.unionResponsive), ...
    %     'FontSize', 9);

    ax = gca;
    ax.YAxis.FontSize = 8; ax.YAxis.FontName = 'helvetica';
    ax.XAxis.FontSize = 8; ax.XAxis.FontName = 'helvetica';
    set(fig, 'Units', 'centimeters', 'Position', [20 20 5 4]);

    if params.PaperFig
        vs_first.printFig(fig, sprintf('SpatialTuningIndex-%s-%s-RF%d-prefDir%d-allResp%d', ...
            params.indexType, strjoin(params.stimTypes, '-'), ...
            params.useRF, params.prefDir, params.allResponsive), ...
            PaperFig = params.PaperFig);
    end

    results.fig = fig;
    results.ps  = ps;

end

% =========================================================================
% RECEPTIVE FIELD PAGES (multi-page PDF)
% Generates one PDF per stim type.  Each page shows a tile-layout of
% full-resolution RFs, sorted by descending tuning index.  Each tile is
% annotated with the neuron's phy ID, tuning-index value, and (for
% linearlyMovingBall) the preferred direction in degrees.
% =========================================================================
if params.plotRFs

    % ---- Guard: prefDir must be true for MB RF pages --------------------
    if ~params.prefDir
        % Without preferred-direction selection the RF slice to plot is
        % ambiguous for linearlyMovingBall, so we skip entirely.
        fprintf(['plotRFs requires prefDir=true for linearlyMovingBall.\n' ...
                 'Skipping RF page generation.\n']);
    else

        % ---- Page geometry ----------------------------------------------
        nCols        = 5;          % number of tile columns per page
        nRows        = 8;          % number of tile rows per page
        tilesPerPage = nCols * nRows;  % total tiles that fit on one page

        % ---- Filter the results table to the requested condition --------
        idxCond_rf       = tbl.onOff   == params.onOff   & ...
                           tbl.sizeIdx == params.sizeIdx & ...
                           tbl.lumIdx  == params.lumIdx;
        tblCondRF        = tbl(idxCond_rf, :);           % rows matching condition
        tblCondRF.value  = tblCondRF.(params.indexType);  % copy chosen index into 'value'

        % ---------------------------------------------------------
        % If plotRFunion: build a single neuron set from the union
        % of neurons responsive to either stim type, BEFORE entering
        % the stim-type loop.  Both PDFs will show the same neurons.
        % ---------------------------------------------------------
        if params.plotRFunion

            % Inline ternary helper for compact annotation formatting:
            % returns trueStr when cond is true, falseStr otherwise.
            ternaryStr = @(cond, trueStr, falseStr) ...
                subsref({falseStr; trueStr}, substruct('{}', {cond + 1}));

            % Separate table rows by stim type for per-stim-type lookup
            mbRows = tblCondRF(tblCondRF.stimulus == "linearlyMovingBall", :);
            sbRows = tblCondRF(tblCondRF.stimulus == "rectGrid", :);

            % Collect every unique (experimentNum, phyID) pair across both
            % stim types — this is the union neuron set
            keysMB  = mbRows(:, {'experimentNum', 'phyID'});
            keysSB  = sbRows(:, {'experimentNum', 'phyID'});
            allKeys = unique([keysMB; keysSB], 'rows');

            % Pre-allocate columns for each stim type's tuning index and
            % preferred-direction metadata (MB only)
            nUnion = height(allKeys);
            allKeys.valueMB    = nan(nUnion, 1);  % MB tuning index  (NaN = not responsive)
            allKeys.valueSB    = nan(nUnion, 1);  % SB tuning index  (NaN = not responsive)
            allKeys.prefDirIdx = nan(nUnion, 1);  % preferred direction index (MB only)
            allKeys.prefDirDeg = nan(nUnion, 1);  % preferred direction degrees (MB only)

            % Fill per-neuron values by matching on experiment + phyID
            for ri = 1:nUnion
                % Look up this neuron in the MB rows
                mMatch = mbRows.experimentNum == allKeys.experimentNum(ri) & ...
                    mbRows.phyID         == allKeys.phyID(ri);
                if any(mMatch)
                    mIdx = find(mMatch, 1);              % first (only) matching row
                    allKeys.valueMB(ri)    = mbRows.value(mIdx);
                    allKeys.prefDirIdx(ri) = mbRows.prefDirIdx(mIdx);
                    allKeys.prefDirDeg(ri) = mbRows.prefDirDeg(mIdx);
                end

                % Look up this neuron in the SB rows
                sMatch = sbRows.experimentNum == allKeys.experimentNum(ri) & ...
                    sbRows.phyID         == allKeys.phyID(ri);
                if any(sMatch)
                    allKeys.valueSB(ri) = sbRows.value(find(sMatch, 1));
                end
            end

            fprintf('  [plotRFunion] %d unique neurons in union set.\n', nUnion);
        end

        % ---- Iterate over each stimulus type ----------------------------
        for ss = 1:numel(params.stimTypes)

            stimType = params.stimTypes(ss);  % current stimulus type string

            % ---------------------------------------------------------
            % Select neurons: union set or per-stim-type set
            % ---------------------------------------------------------
            if params.plotRFunion

                % Use the pre-built union table — same neurons for both PDFs
                tblStim = allKeys;

                % Sort both PDFs by MB tuning index so neuron positions
                % match across the two PDFs for direct visual comparison.
                % Neurons not responsive to MB (NaN) sink to the bottom.
                tblStim.value = tblStim.valueMB;
                tblStim = sortrows(tblStim, 'value', 'descend', ...
                    'MissingPlacement', 'last');
            else
                % Default: only neurons responsive to this stim type
                stimMask = tblCondRF.stimulus == char(stimType);
                tblStim  = tblCondRF(stimMask, :);
            end

            % Skip if no data for this stim type
            if isempty(tblStim)
                fprintf('  [plotRFs] No neurons for %s — skipping.\n', char(stimType));
                continue
            end

            % Sort by tuning index descending (only for the non-union path;
            % the union path was already sorted above)
            if ~params.plotRFunion
                tblStim = sortrows(tblStim, 'value', 'descend');
            end

            nNeurons = height(tblStim);                % total neurons to plot
            nPages   = ceil(nNeurons / tilesPerPage);  % number of PDF pages

            % Build the output PDF path inside the combined-analysis directory
            if params.subtractShuffle
                shuffTag = '_shuffSub';   % tag indicating shuffle was subtracted
            else
                shuffTag = '_raw';        % tag indicating raw RF plotted
            end
            unionTag = '';  if params.plotRFunion, unionTag = '_union'; end
            pdfName = sprintf('RFpages_%s_%s%s%s.pdf', ...
                char(stimType), params.indexType, shuffTag, unionTag);
            pdfPath = fullfile(saveDir, pdfName);  % full path for the PDF

            % Cache loaded experiments so each experiment's heavy data is
            % read from disk only once (key = experiment number)
            cachedExp = containers.Map('KeyType', 'double', 'ValueType', 'any');

            % ---- Loop over pages ----------------------------------------
            for pg = 1:nPages

                % Create an invisible figure sized to A4 portrait (21 × 29.7 cm)
                fig_rf = figure('Visible', 'off', ...
                    'Units', 'centimeters', ...
                    'Position', [0 0 21 29.7]);

                % Tiled layout with compact spacing to maximise tile area
                tl = tiledlayout(nRows, nCols, ...
                    'TileSpacing', 'compact', ...
                    'Padding',     'compact');

                % Page-level title indicating stim type, mode, and page number
                if stimType == "linearlyMovingBall"
                    stimLabel_pg = 'Moving Ball RFs';
                else
                    stimLabel_pg = 'Rect Grid RFs';
                end
                if params.plotRFunion
                    stimLabel_pg = [stimLabel_pg ' (union)'];  %#ok<AGROW>
                end
                pageTitleStr = sprintf('%s — %s  (page %d/%d)', ...
                    stimLabel_pg, params.indexType, pg, nPages);
                title(tl, pageTitleStr, 'FontSize', 9, 'FontName', 'Helvetica');

                % Index range for the neurons that belong to this page
                startNeuron = (pg - 1) * tilesPerPage + 1;
                endNeuron   = min(pg * tilesPerPage, nNeurons);

                % ---- Loop over neurons on this page ---------------------
                for ni = startNeuron:endNeuron

                    nexttile;  % advance to the next tile in the layout

                    % Read metadata for this neuron from the sorted table
                    phyID = tblStim.phyID(ni);
                    tVal  = tblStim.value(ni);
                    ex    = str2double(string(tblStim.experimentNum(ni)));

                    % ----- Load experiment data (with caching) -----------
                    if ~cachedExp.isKey(ex)

                        NP_tmp = loadNPclassFromTable(ex);

                        % Always derive phy_IDg from linearlyMovingBallAnalysis
                        % to match the phyIDs stored in the results table
                        obj_lmb     = linearlyMovingBallAnalysis(NP_tmp);
                        p_s_tmp     = NP_tmp.convertPhySorting2tIc(obj_lmb.spikeSortingFolder);
                        phy_IDg_tmp = p_s_tmp.phy_ID(string(p_s_tmp.label') == 'good');

                        % Load RF data from the stim-specific analysis object
                        switch stimType
                            case "linearlyMovingBall"
                                obj_stim = obj_lmb;
                            case "rectGrid"
                                obj_stim = rectGridAnalysis(NP_tmp);
                        end
                        S_rf_tmp = obj_stim.CalculateReceptiveFields;

                        % -------------------------------------------------
                        % Compute preferred direction for ALL good units
                        % (needed for union mode, where some neurons may not
                        %  be MB-responsive and thus lack prefDirIdx in the
                        %  table).  Uses the same logic as the main computation.
                        % -------------------------------------------------
                        S_rf_lmb   = obj_lmb.CalculateReceptiveFields;
                        rfSpeed    = S_rf_lmb.params.speed;
                        rfField    = sprintf('Speed%d', rfSpeed);
                        NeuronResp = obj_lmb.ResponseWindow;
                        NeuronVals = NeuronResp.(rfField).NeuronVals;
                        % NeuronVals: [nGoodUnits, nConditions, nFeatures]

                        dirLabels  = NeuronVals(:,:,6);         % direction (radians)
                        spikeRates = NeuronVals(:,:,1);         % spike rate
                        uDirsLocal = unique(dirLabels(1,:));    % sorted unique directions

                        % Max spike rate per direction per good unit
                        nGU = size(NeuronVals, 1);
                        maxRPD = zeros(nGU, numel(uDirsLocal));
                        for dd = 1:numel(uDirsLocal)
                            dMask = dirLabels(1,:) == uDirsLocal(dd);
                            maxRPD(:,dd) = max(spikeRates(:, dMask), [], 2);
                        end
                        [~, prefDirIdxAll] = max(maxRPD, [], 2);        % [nGU,1]
                        prefDirDegAll = rad2deg(uDirsLocal(prefDirIdxAll))';  % [nGU,1]

                        % Store everything in the cache
                        cachedExp(ex) = struct( ...
                            'S_rf',           S_rf_tmp, ...
                            'S_rf_lmb',       S_rf_lmb, ...
                            'phy_IDg',        phy_IDg_tmp, ...
                            'prefDirIdxAll',  prefDirIdxAll, ...
                            'prefDirDegAll',  prefDirDegAll);
                    end

                    % Retrieve cached data for this experiment
                    cached = cachedExp(ex);

                    % Find this neuron's index within the good-unit array
                    [~, nIdx] = ismember(phyID, cached.phy_IDg);

                    % Guard: if phyID is not found, skip this tile
                    if nIdx == 0
                        text(0.5, 0.5, sprintf('phy%d\nnot found', phyID), ...
                            'HorizontalAlignment', 'center', ...
                            'FontSize', 5);
                        axis off;
                        continue
                    end

                    % ----- Extract the full-resolution RF slice -----------
                    switch stimType

                        case "linearlyMovingBall"
                            % Preferred direction: use table value if available
                            % (MB-responsive neuron), otherwise use the cached
                            % value computed for all good units
                            prefIdx = tblStim.prefDirIdx(ni);
                            prefDeg = tblStim.prefDirDeg(ni);

                            if isnan(prefIdx)
                                % Neuron was not MB-responsive — use cached
                                % preferred direction computed from NeuronVals
                                prefIdx = cached.prefDirIdxAll(nIdx);
                                prefDeg = cached.prefDirDegAll(nIdx);
                            end

                            % Slice RFuSTDirSizeLum: [nDir,nSize,nLum,rfY,rfX,nN]
                            % Use the MB-specific RF cache (S_rf_lmb) since
                            % S_rf may be rectGrid when this is the SB iteration
                            rfSlice = squeeze( ...
                                cached.S_rf_lmb.RFuSTDirSizeLum( ...
                                prefIdx, params.sizeIdx, params.lumIdx, :, :, nIdx));

                            if params.subtractShuffle
                                rfShuff = cached.S_rf_lmb.RFuShuffST(:, :, nIdx);
                                rfSlice = rfSlice - rfShuff;
                            end

                        case "rectGrid"
                            % Slice RFu: [2(onOff), nLums, nSize, sR, sR, nN]
                            rfSlice = squeeze( ...
                                cached.S_rf.RFu( ...
                                params.onOff, params.lumIdx, params.sizeIdx, :, :, nIdx));

                            if params.subtractShuffle
                                rfShuff = squeeze( ...
                                    cached.S_rf.RFuShuffMean( ...
                                    params.onOff, params.lumIdx, params.sizeIdx, :, :, nIdx));
                                rfSlice = rfSlice - rfShuff;
                            end
                    end

                    % ----- Plot the RF as a colour image -----------------
                    imagesc(rfSlice);
                    axis equal tight;
                    axis off;

                    if params.subtractShuffle
                        maxAbs = max(abs(rfSlice(:)));
                        if maxAbs > 0
                            clim([-maxAbs, maxAbs]);
                        end
                    end

                    cb = colorbar;
                    cb.FontSize  = 3.5;
                    cb.TickDirection = 'out';
                    cb.Ticks = linspace(cb.Limits(1), cb.Limits(2), 3);

                    % ----- Tile title annotation -------------------------
                    if params.plotRFunion
                        % Union mode: show both stim type indices on every tile
                        mbVal = tblStim.valueMB(ni);
                        sbVal = tblStim.valueSB(ni);
                        mbStr = ternaryStr(isnan(mbVal), 'n/r', sprintf('%.2f', mbVal));
                        sbStr = ternaryStr(isnan(sbVal), 'n/r', sprintf('%.2f', sbVal));

                        switch stimType
                            case "linearlyMovingBall"
                                % Also show preferred direction
                                tileTitle = sprintf('phy%d|MB %s|SB %s|%d°', ...
                                    phyID, mbStr, sbStr, round(prefDeg));
                            case "rectGrid"
                                tileTitle = sprintf('phy%d|MB %s|SB %s', ...
                                    phyID, mbStr, sbStr);
                        end
                    else
                        % Standard mode: show only this stim type's index
                        switch stimType
                            case "linearlyMovingBall"
                                tileTitle = sprintf('phy%d | %.2f | %d°', ...
                                    phyID, tVal, round(tblStim.prefDirDeg(ni)));
                            case "rectGrid"
                                tileTitle = sprintf('phy%d | %.2f', phyID, tVal);
                        end
                    end
                    title(tileTitle, 'FontSize', 5, 'FontName', 'Helvetica');

                end  % neuron loop (tiles on this page)

                % ----- Export this page to the PDF -----------------------
                if pg == 1
                    % First page: create the PDF file
                    exportgraphics(fig_rf, pdfPath, 'ContentType', 'vector');
                else
                    % Subsequent pages: append to the existing PDF
                    exportgraphics(fig_rf, pdfPath, 'ContentType', 'vector', 'Append', true);
                end

                close(fig_rf);  % close figure to free memory before next page

                fprintf('  [plotRFs] %s — page %d/%d exported.\n', ...
                    char(stimType), pg, nPages);

            end  % page loop

            fprintf('  [plotRFs] Saved %d pages to:\n    %s\n', nPages, pdfPath);

        end  % stim-type loop

    end  % prefDir guard

end  % plotRFs block

end