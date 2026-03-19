function [fig, tbl] = plotSpatialTuningIndex(exList, pairs, params)

arguments
    exList   double
    pairs    cell   = {}
    params.stimTypes    (1,:) string = ["rectGrid", "linearlyMovingBall"]
    params.indexType    string       = "L_combined"  % L_amplitude, L_geometric, L_combined
    params.onOff        double       = 1             % 1=on, 2=off (rectGrid only)
    params.sizeIdx      double       = 1
    params.lumIdx       double       = 1
    params.nBoot        double       = 10000
    params.overwrite    logical      = false
    params.yLegend      char         = 'Spatial Tuning Index'
    params.yMaxVis      double       = 1
    params.Alpha        double       = 0.4
    params.PaperFig     logical      = false
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
saveDir = [p '\Combined_lizard_analysis'];

stimLabel  = strjoin(params.stimTypes, '-');
nameOfFile = sprintf('\\Ex_%d-%d_SpatialTuningIndex_%s.mat', ...
    exList(1), exList(end), stimLabel);

% -------------------------------------------------------------------------
% Load SpatialTuningIndex results
% -------------------------------------------------------------------------
if ~exist([saveDir nameOfFile], 'file')
    error('SpatialTuningIndex results not found. Run SpatialTuningIndex first.');
end

S = load([saveDir nameOfFile]);

% -------------------------------------------------------------------------
% Build long table
% -------------------------------------------------------------------------
tbl = table();

nExp  = numel(exList);
nStim = numel(params.stimTypes);

for ei = 1:nExp
    ex = exList(ei);

    % Get animal/insertion info
    try
        NP = loadNPclassFromTable(ex);
    catch
        warning('Could not load experiment %d — skipping.', ex);
        continue
    end

    for s = 1:nStim

        stimType = params.stimTypes(s);

        % Get the right index matrix for this stim/exp
        switch params.indexType
            case "L_amplitude"
                idxMat = S.L_amplitude_all{s, ei};
            case "L_geometric"
                idxMat = S.L_geometric_all{s, ei};
            case "L_combined"
                idxMat = S.L_combined_all{s, ei};
        end

        if isempty(idxMat)
            continue
        end

        % idxMat is [nN x nOnOff x nSize x nLum]
        % for linearlyMovingBall there is no onOff dimension — handle both
        if ndims(idxMat) == 3
            % [nN x nSize x nLum] — no onOff
            vals = idxMat(:, params.sizeIdx, params.lumIdx);
            oi   = 1;
        else
            vals = idxMat(:, params.onOff, params.sizeIdx, params.lumIdx);
            oi   = params.onOff;
        end

        nN = numel(vals);

        % Build rows for this experiment/stim
        rows            = table();
        rows.value      = vals;
        rows.stimulus   = categorical(repmat({char(stimType)}, nN, 1));
        rows.insertion  = categorical(repmat(ex,              nN, 1));
        rows.animal     = categorical(repmat({NP.animalName}, nN, 1));
        rows.NeurID     = (1:nN)';
        rows.onOff      = repmat(oi,                nN, 1);
        rows.sizeIdx    = repmat(params.sizeIdx,    nN, 1);
        rows.lumIdx     = repmat(params.lumIdx,     nN, 1);
        rows.indexType  = categorical(repmat({params.indexType}, nN, 1));

        tbl = [tbl; rows];
    end
end

if isempty(tbl)
    warning('No data found — table is empty.');
    fig = [];
    return
end

% Clean up categories
tbl.stimulus  = removecats(tbl.stimulus);
tbl.animal    = removecats(tbl.animal);
tbl.insertion = removecats(tbl.insertion);

% -------------------------------------------------------------------------
% Compute p-values using hierBoot
% -------------------------------------------------------------------------
ps = [];

if ~isempty(pairs)
    ps  = zeros(size(pairs, 1), 1);
    j   = 1;

    for i = 1:size(pairs, 1)
        diffs   = [];
        insers  = [];
        animals = [];

        for ins = unique(tbl.insertion)'
            idx1 = tbl.insertion == categorical(ins) & tbl.stimulus == pairs{i,1};
            idx2 = tbl.insertion == categorical(ins) & tbl.stimulus == pairs{i,2};

            V1 = tbl.value(idx1);
            V2 = tbl.value(idx2);

            if isempty(V1) || isempty(V2)
                continue
            end

            animal = unique(tbl.animal(idx1));
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
end

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
V1max = max(tbl.value, [], 'omitnan');

[fig, ~] = plotSwarmBootstrapWithComparisons(tbl, pairs, ps, {'value'}, ...
    yLegend  = params.yLegend, ...
    yMaxVis  = max(params.yMaxVis, V1max), ...
    diff     = false, ...
    Alpha    = params.Alpha, ...
    plotMeanSem = true);

title(sprintf('%s — %s  (size=%d, lum=%d)', ...
    params.indexType, strjoin(params.stimTypes,'/'), ...
    params.sizeIdx, params.lumIdx), ...
    'FontSize', 9);

if params.PaperFig
    vs_first.printFig(fig, sprintf('SpatialTuningIndex-%s-%s', ...
        params.indexType, strjoin(params.stimTypes, '-')), ...
        PaperFig = params.PaperFig);
end

end