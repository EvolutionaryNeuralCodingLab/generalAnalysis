function results = SpatialTuningIndex(exList, params)

arguments
    exList   double
    params.stimTypes  (1,:) string  = ["rectGrid", "linearlyMovingBall"]
    params.topPercent double        = 10
    params.overwrite  logical       = false
    params.statType   string        = "BootstrapPerNeuron"
    params.speed      double        = 1
    params.plot       logical       = true
    params.indexType  string        = "L_amplitude"  % L_amplitude_diff,L_amplitude_ratio, L_geometric, L_combined
    params.onOff      double        = 1             % 1=on, 2=off (rectGrid only)
    params.sizeIdx    double        = 1
    params.lumIdx     double        = 1
    params.nBoot      double        = 10000
    params.yLegend    char          = 'Spatial Tuning Index'
    params.yMaxVis    double        = 1
    params.Alpha      double        = 0.4
    params.PaperFig   logical       = false
    params.useRF      logical       = false
end

% -------------------------------------------------------------------------
% Build save path
% -------------------------------------------------------------------------
NP_first = loadNPclassFromTable(exList(1));

switch params.stimTypes(1)
    case "rectGrid"
        vs_first = rectGridAnalysis(NP_first);
    case "linearlyMovingBall"
        vs_first = linearlyMovingBallAnalysis(NP_first);
end

p = extractBefore(vs_first.getAnalysisFileName, 'lizards');
p = [p 'lizards'];
if ~exist([p '\Combined_lizard_analysis'], 'dir')
    cd(p)
    mkdir Combined_lizard_analysis
end
saveDir = [p '\Combined_lizard_analysis'];

stimLabel  = strjoin(params.stimTypes, '-');
nameOfFile = sprintf('\\Ex_%d-%d_SpatialTuningIndex_%s.mat', ...
    exList(1), exList(end), stimLabel);

% -------------------------------------------------------------------------
% Decide whether to compute or load
% -------------------------------------------------------------------------
if exist([saveDir nameOfFile], 'file') == 2 && ~params.overwrite
    S = load([saveDir nameOfFile]);
    if isequal(S.expList, exList)
        fprintf('Loading saved SpatialTuningIndex from:\n  %s\n', [saveDir nameOfFile]);
        % Jump straight to table building
        tbl = S.tbl;
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

    tbl = table();

    for ei = 1:nExp

        ex = exList(ei);
        fprintf('\n=== Experiment %d ===\n', ex);

        try
            NP = loadNPclassFromTable(ex);
        catch ME
            warning('Could not load experiment %d: %s', ex, ME.message);
            continue
        end

        obj_s = linearlyMovingBallAnalysis(NP);

        nameParts  = split(NP.recordingName, '_');
        animalName = nameParts{1};

        % ----------------------------------------------------------
        % Find union of responsive neurons across ALL stim types
        % ----------------------------------------------------------

        % Get phy IDs once — same for all stim types
        p_s     = NP.convertPhySorting2tIc(obj_s.spikeSortingFolder);
        phy_IDg = p_s.phy_ID(string(p_s.label') == 'good');

        respPhyIDs_all = cell(1, nStim);
        respU_all      = cell(1, nStim);  % ADD — stores respU indices per stim

        for s = 1:nStim
            stimType = params.stimTypes(s);
            try
                switch stimType
                    case "rectGrid"
                        obj_s = rectGridAnalysis(NP);
                    case "linearlyMovingBall"
                        obj_s = linearlyMovingBallAnalysis(NP);
                end

                if params.statType == "BootstrapPerNeuron"
                    Stats = obj_s.BootstrapPerNeuron;
                else
                    Stats = obj_s.ShufflingAnalysis;
                end

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

                respU             = find(pvals < 0.05);
                respU_all{s}      = respU;           % ADD — index into gridSpikeRate dim 3
                respPhyIDs_all{s} = phy_IDg(respU);  % phy IDs of responsive neurons
                fprintf('  [%s] %d responsive neuron(s).\n', stimType, numel(respU));

            catch ME
                warning('Could not get pvals for %s exp %d: %s', stimType, ex, ME.message);
                respU_all{s}      = [];
                respPhyIDs_all{s} = [];
            end
        end

        % Intersection of responsive phy IDs across stim types
        sharedPhyIDs = respPhyIDs_all{1};
        for s = 2:nStim
            sharedPhyIDs = intersect(sharedPhyIDs, respPhyIDs_all{s});
        end

        if isempty(sharedPhyIDs)
            fprintf('  No neurons responsive to all stim types in exp %d — skipping.\n', ex);
            continue
        end

        fprintf('  %d neuron(s) responsive to all stim types in exp %d.\n', numel(sharedPhyIDs), ex);


        for s = 1:nStim

            stimType = params.stimTypes(s);

            % Build analysis object
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


            % ----------------------------------------------------------
            % Load grid results
            % ----------------------------------------------------------
            S_rf = obj.CalculateReceptiveFields;

            gridSpikeRate      = S_rf.gridSpikeRate;
            gridSpikeRateShuff = S_rf.gridSpikeRateShuff;

            switch stimType
                case "rectGrid"
                    gridSpikeRateSelected = gridSpikeRate(:,:,:,params.onOff,:,:);
                    gridShuffSelected     = gridSpikeRateShuff(:,:,:,:,params.onOff,:,:);

                    % Remove onOff singleton at dim 4 for rate: [9 9 nN 1 nSize nLum] -> [9 9 nN nSize nLum]
                    gridSpikeRateSelected = reshape(gridSpikeRateSelected, ...
                        [size(gridSpikeRateSelected,1), size(gridSpikeRateSelected,2), ...
                        size(gridSpikeRateSelected,3), size(gridSpikeRateSelected,5), ...
                        size(gridSpikeRateSelected,6)]);

                    % Remove onOff singleton at dim 5 for shuff: [9 9 nN nShuffle 1 nSize nLum] -> [9 9 nN nShuffle nSize nLum]
                    gridShuffSelected = reshape(gridShuffSelected, ...
                        [size(gridShuffSelected,1), size(gridShuffSelected,2), ...
                        size(gridShuffSelected,3), size(gridShuffSelected,4), ...
                        size(gridShuffSelected,6), size(gridShuffSelected,7)]);
                case "linearlyMovingBall"
                    gridSpikeRateSelected = gridSpikeRate;       % [nGrid nGrid nN nSize nLum]
                    gridShuffSelected     = gridSpikeRateShuff;  % [nGrid nGrid nN nShuffle nSize nLum]
            end

            % Find which indices of THIS stim's gridSpikeRate correspond to sharedPhyIDs
            [~, neuronIdx] = ismember(sharedPhyIDs, phy_IDg(respU_all{s}));

            gridSpikeRateSelected = gridSpikeRateSelected(:,:,neuronIdx,:,:);
            gridShuffSelected     = gridShuffSelected(:,:,neuronIdx,:,:,:);
           
            % Average over shuffles and reshape explicitly — no squeeze
            gridShuffMean = mean(gridShuffSelected, 4);  % [nGrid nGrid nN 1 nSize nLum]

            % Get dimensions explicitly
            nN    = size(gridSpikeRateSelected, 3);
            nSize = size(gridSpikeRateSelected, 4);
            nLum  = size(gridSpikeRateSelected, 5);
            nGrid = size(gridSpikeRateSelected, 1);

            fprintf('gridSpikeRateSelected size before reshape: %s\n', num2str(size(gridSpikeRateSelected)));
            fprintf('Expected: [%d %d %d %d %d]\n', nGrid, nGrid, nN, nSize, nLum);

            % Reshape both to clean [nGrid nGrid nN nSize nLum]
            gridSpikeRateSelected = reshape(gridSpikeRateSelected, [nGrid nGrid nN nSize nLum]);
            gridShuffMean         = reshape(gridShuffMean,         [nGrid nGrid nN nSize nLum]);

            nCells = nGrid * nGrid;
            maxDist  = sqrt(2) * (nGrid - 1);

            % Average over shuffles

            % ----------------------------------------------------------
            % Compute indices
            % ----------------------------------------------------------

            fprintf('gridSpikeRate size: %s\n', num2str(size(gridSpikeRate)));
            fprintf('gridSpikeRateShuff size: %s\n', num2str(size(gridSpikeRateShuff)));
            fprintf('gridShuffMean size: %s\n', num2str(size(gridShuffMean)));

            for si = 1:nSize
                for li = 1:nLum

                    rateFlat      = reshape(gridSpikeRateSelected(:,:,:,si,li), [nCells, nN]);
                    rateFlatShuff = reshape(gridShuffMean(:,:,:,si,li),          [nCells, nN]);

                    L_amplitude_diff = zeros(nN, 1);
                    L_amplitude_ratio = zeros(nN, 1);
                    L_geometric = zeros(nN, 1);
                    L_combined  = zeros(nN, 1);

                    for u = 1:nN

                        rateVec      = rateFlat(:, u);
                        rateVecShuff = rateFlatShuff(:, u);

                        % Top cells
                        threshold      = prctile(rateVec,      100 - params.topPercent);
                        thresholdShuff = prctile(rateVecShuff, 100 - params.topPercent);

                        topIdx      = find(rateVec      >= threshold);
                        topIdxShuff = find(rateVecShuff >= thresholdShuff);
                        restIdx      = setdiff(1:nCells, topIdx);
                        restIdxShuff = setdiff(1:nCells, topIdxShuff);

                        % Amplitude
                        meanTop       = mean(rateVec(topIdx));
                        meanRest      = mean(rateVec(restIdx));
                        meanAll       = mean(rateVec);
                        meanTopShuff  = mean(rateVecShuff(topIdxShuff));
                        meanRestShuff = mean(rateVecShuff(restIdxShuff));
                        meanAllShuff  = mean(rateVecShuff);

                        if meanAll      == 0, meanAll      = eps; end
                        if meanAllShuff == 0, meanAllShuff = eps; end

                        L_amplitude_diff(u) = ...
                            (meanTop - meanRest) / meanAll - ...
                            (meanTopShuff - meanRestShuff) / meanAllShuff;

                        shuffleNorm = (meanTopShuff - meanRestShuff) / meanAllShuff;
                        if shuffleNorm == 0, shuffleNorm = eps; end

                        L_amplitude_ratio(u) = ((meanTop - meanRest) / meanAll) / shuffleNorm;

                        % Geometric
                        [rowIdx,      colIdx]      = ind2sub([nGrid nGrid], topIdx);
                        [rowIdxShuff, colIdxShuff] = ind2sub([nGrid nGrid], topIdxShuff);

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

                        L_geometric(u) = (1 - D) - (1 - DShuff);
                        L_combined(u)  = L_amplitude_diff(u) * L_geometric(u);

                    end

                    % Build rows for this condition
                    rows             = table();
                    rows.L_amplitude_diff = L_amplitude_diff;
                    rows.L_amplitude_ratio = L_amplitude_ratio;
                    rows.L_geometric = L_geometric;
                    rows.L_combined  = L_combined;
                    rows.stimulus    = categorical(repmat({char(stimType)}, nN, 1));
                    rows.insertion   = categorical(repmat(ex,              nN, 1));
                    rows.animal      = categorical(repmat({animalName},    nN, 1));
                    rows.NeurID      = (1:nN)';
                    rows.onOff       = repmat(params.onOff, nN, 1); % params.onOff for rectGrid, meaningless but consistent for movingBall
                    rows.sizeIdx     = repmat(si, nN, 1);
                    rows.lumIdx      = repmat(li, nN, 1);

                    tbl = [tbl; rows];

                end
            end

            fprintf('  [%s] Indices computed. %d neurons.\n', stimType, nN);

        end % stim loop
    end % exp loop

    % Clean categories
    tbl.stimulus  = removecats(tbl.stimulus);
    tbl.animal    = removecats(tbl.animal);
    tbl.insertion = removecats(tbl.insertion);

    % Save
    S.expList = exList;
    S.tbl     = tbl;
    S.params  = params;
    save([saveDir nameOfFile], '-struct', 'S');
    fprintf('\nSaved to:\n  %s\n', [saveDir nameOfFile]);

end % compute block

results.tbl = tbl;

% =========================================================================
% PLOT
% =========================================================================
if params.plot

    % Filter table to requested condition
    idx = tbl.onOff   == params.onOff   & ...
        tbl.sizeIdx == params.sizeIdx & ...
        tbl.lumIdx  == params.lumIdx;

    tblPlot = tbl(idx, :);
    tblPlot.value = tblPlot.(params.indexType);  % select which index to plot

    % ----------------------------------------------------------
    % Compute p-values using hierBoot
    % ----------------------------------------------------------
    ps = [];

    pairs = {char(params.stimTypes(1)), char(params.stimTypes(2))};


    ps = zeros(size(pairs, 1), 1);
    j  = 1;

    for i = 1:size(pairs, 1)
        diffs   = [];
        insers  = [];
        animals = [];

        for ins = unique(tblPlot.insertion)'
            idx1 = tblPlot.insertion == categorical(ins) & tblPlot.stimulus == pairs{i,1};
            idx2 = tblPlot.insertion == categorical(ins) & tblPlot.stimulus == pairs{i,2};

            V1 = tblPlot.value(idx1);
            V2 = tblPlot.value(idx2);

            if isempty(V1) || isempty(V2)
                continue
            end

            animal  = unique(tblPlot.animal(idx1));
            diffs   = [diffs;   V1 - V2];
            insers  = [insers;  double(repmat(ins,    size(V1,1), 1))];
            animals = [animals; double(repmat(animal, size(V1,1), 1))];
        end

        if isempty(diffs)
            ps(j) = NaN;
        else
            bootDiff = hierBoot(diffs, params.nBoot, insers, animals);
            ps(j)    = mean(bootDiff <= 0);
        end
        j = j + 1;
    end


    % ----------------------------------------------------------
    % Plot
    % ----------------------------------------------------------
    V1max = max(tblPlot.value, [], 'omitnan');

    fprintf('Length of ps: %d\n', numel(ps));
    fprintf('Size of pairs: %s\n', num2str(size(pairs)));

    [fig, ~] = plotSwarmBootstrapWithComparisons(tblPlot, pairs, ps, {'value'}, ...
        yLegend     = params.yLegend, ...
        yMaxVis     = max(params.yMaxVis, V1max), ...
        diff        = true, ...
        Alpha       = params.Alpha, ...
        plotMeanSem = true);

    title(sprintf('%s — %s  (onOff=%d, size=%d, lum=%d)', ...
        params.indexType, strjoin(params.stimTypes, '/'), ...
        params.onOff, params.sizeIdx, params.lumIdx), ...
        'FontSize', 9);

    if params.PaperFig
        vs_first.printFig(fig, sprintf('SpatialTuningIndex-%s-%s', ...
            params.indexType, strjoin(params.stimTypes, '-')), ...
            PaperFig = params.PaperFig);
    end

    results.fig = fig;
    results.ps  = ps;

end

end