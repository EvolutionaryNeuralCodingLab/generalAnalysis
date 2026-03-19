function results = SpatialTuningIndex(exList, params)

arguments
    exList double
    params.stimTypes  (1,:) string  = ["rectGrid", "linearlyMovingBall"]
    params.topPercent double        = 10
    params.overwrite  logical       = false
    params.statType   string        = "BootstrapPerNeuron"
    params.speed      double        = 1
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
        results = S;
        return
    else
        fprintf('Experiment list mismatch — recomputing.\n');
    end
end

% -------------------------------------------------------------------------
% EXPERIMENT LOOP
% -------------------------------------------------------------------------
nExp  = numel(exList);
nStim = numel(params.stimTypes);

% Will grow as we discover dimensions from first valid experiment
L_amplitude_all = cell(nStim, nExp);
L_geometric_all = cell(nStim, nExp);
L_combined_all  = cell(nStim, nExp);

for ei = 1:nExp

    ex = exList(ei);
    fprintf('\n=== Experiment %d ===\n', ex);

    try
        NP = loadNPclassFromTable(ex);
    catch ME
        warning('Could not load experiment %d: %s', ex, ME.message);
        continue
    end

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
        % Check for responsive neurons
        % ----------------------------------------------------------
        try
            if params.statType == "BootstrapPerNeuron"
                Stats = obj.BootstrapPerNeuron;
            else
                Stats = obj.ShufflingAnalysis;
            end

            p_sort = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
            label  = string(p_sort.label');
            goodU  = p_sort.ic(:, label == 'good');

            % Resolve field name depending on stim type
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

            respU = find(pvals < 0.05);

        catch ME
            warning('Could not load stats for %s exp %d: %s', stimType, ex, ME.message);
            L_amplitude_all{s, ei} = [];
            L_geometric_all{s, ei} = [];
            L_combined_all{s, ei}  = [];
            continue
        end

        if isempty(respU)
            fprintf('  [%s] No responsive neurons in exp %d — skipping.\n', stimType, ex);
            L_amplitude_all{s, ei} = [];
            L_geometric_all{s, ei} = [];
            L_combined_all{s, ei}  = [];
            continue
        end

        fprintf('  [%s] %d responsive neuron(s) in exp %d.\n', stimType, ex, numel(respU));

        % Load grid results

        S_rf = obj.CalculateReceptiveFields;

        gridSpikeRate      = S_rf.gridSpikeRate;       % [nGrid x nGrid x nN x onOff x nSize x nLum]
        gridSpikeRateShuff = S_rf.gridSpikeRateShuff;  % [nGrid x nGrid x nN x nShuffle x nSize x nLum]

        [nGrid, ~, nN, nOnOff, nSize, nLum] = size(gridSpikeRate);
        nShuffle = size(gridSpikeRateShuff, 4);
        nCells   = nGrid * nGrid;

        % Average over shuffles
        gridShuffMean = mean(gridSpikeRateShuff, 4); % [nGrid x nGrid x nN x nSize x nLum]

        L_amplitude = zeros(nN, nOnOff, nSize, nLum);
        L_geometric = zeros(nN, nOnOff, nSize, nLum);
        L_combined  = zeros(nN, nOnOff, nSize, nLum);

        maxDist = sqrt(2) * (nGrid - 1);

        for oi = 1:nOnOff
            for si = 1:nSize
                for li = 1:nLum

                    rateFlat      = reshape(gridSpikeRate(:,:,:,oi,si,li),  [nCells, nN]);
                    rateFlatShuff = reshape(gridShuffMean(:,:,:,si,li),      [nCells, nN]);

                    for u = 1:nN

                        rateVec      = rateFlat(:, u);
                        rateVecShuff = rateFlatShuff(:, u);

                        %% ---- Shared: top cells ----
                        threshold      = prctile(rateVec,      100 - params.topPercent);
                        thresholdShuff = prctile(rateVecShuff, 100 - params.topPercent);

                        topIdx      = find(rateVec      >= threshold);
                        topIdxShuff = find(rateVecShuff >= thresholdShuff);

                        restIdx      = setdiff(1:nCells, topIdx);
                        restIdxShuff = setdiff(1:nCells, topIdxShuff);

                        %% ---- 1. Amplitude index ----
                        meanTop      = mean(rateVec(topIdx));
                        meanRest     = mean(rateVec(restIdx));
                        meanAll      = mean(rateVec);

                        meanTopShuff  = mean(rateVecShuff(topIdxShuff));
                        meanRestShuff = mean(rateVecShuff(restIdxShuff));
                        meanAllShuff  = mean(rateVecShuff);

                        if meanAll      == 0, meanAll      = eps; end
                        if meanAllShuff == 0, meanAllShuff = eps; end

                        L_amplitude(u, oi, si, li) = ...
                            (meanTop - meanRest) / meanAll - ...
                            (meanTopShuff - meanRestShuff) / meanAllShuff;

                        %% ---- 2. Geometric index ----
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

                        L_geometric(u, oi, si, li) = (1 - D) - (1 - DShuff);

                        %% ---- 3. Combined index ----
                        L_combined(u, oi, si, li) = L_amplitude(u, oi, si, li) * L_geometric(u, oi, si, li);

                    end
                end
            end
        end

        L_amplitude_all{s, ei} = L_amplitude;  % [nN x nOnOff x nSize x nLum]
        L_geometric_all{s, ei} = L_geometric;
        L_combined_all{s, ei}  = L_combined;

        fprintf('  [%s] Done. %d neurons.\n', stimType, nN);

    end % stim loop
end % experiment loop

% -------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------
S.expList         = exList;
S.L_amplitude_all = L_amplitude_all;  % {nStim x nExp} cell, each [nN x nOnOff x nSize x nLum]
S.L_geometric_all = L_geometric_all;
S.L_combined_all  = L_combined_all;
S.params          = params;

save([saveDir nameOfFile], '-struct', 'S');
fprintf('\nSaved SpatialTuningIndex to:\n  %s\n', [saveDir nameOfFile]);

results = S;

end