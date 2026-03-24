function [result] = getNeuronDepths(exList)
% getNeuronDepths  Returns cortical depths of good units across all experiments,
%                 and computes 3 globally-defined equal depth bins.
%
% Inputs:
%   exList  - vector of experiment numbers (same as used in plotPSTH_MultiExp)
%
% Outputs:
%   result  - struct with fields:
%               .depthTable   - table with columns: Experiment, Unit, Depth_um
%               .depthBinEdges - 1x4 vector [min, t1, t2, max] in um
%               .perExp       - struct array with per-experiment data:
%                                 .exNum, .goodU, .p_sort

% ------------------------------------------------------------------
% Load Excel once
% ------------------------------------------------------------------
excelPath = '\\sil3\data\Large_scale_mapping_NP\Experiment_Excel.xlsx';
T = readtable(excelPath);

% ------------------------------------------------------------------
% Preallocate collections
% ------------------------------------------------------------------
expCol   = [];   % experiment number per unit
unitCol  = [];   % unit index (1-based) per unit
depthCol = [];   % depth in um per unit

result.perExp(numel(exList)) = struct('exNum', [], 'goodU', [], 'p_sort', []);

% ------------------------------------------------------------------
% Loop over experiments
% ------------------------------------------------------------------
for ei = 1:numel(exList)

    ex = exList(ei);
    fprintf('Loading experiment %d ...\n', ex);

    try
        NP     = loadNPclassFromTable(ex);
        obj    = linearlyMovingBallAnalysis(NP);
        p_sort = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
    catch ME
        warning('Could not load experiment %d: %s', ex, ME.message);
        result.perExp(ei).exNum  = ex;
        result.perExp(ei).goodU  = [];
        result.perExp(ei).p_sort = [];
        continue
    end

    % coor_Z for this experiment
    coor_Z = T.coor_Z(ex);

    % Good units
    label = string(p_sort.label');
    goodU = p_sort.ic(:, label == 'good');   % nTimePoints x nGoodUnits
    nGood = size(goodU, 2);

    % Channel IDs (0-based) → Y positions → real depths
    channelIDs   = goodU(1, :);                               % 1 x nGoodUnits, 0-based
    yPos         = NP.chLayoutPositions(2, channelIDs + 1);   % 1 x nGoodUnits
    neuronDepths = coor_Z - yPos;                             % 1 x nGoodUnits, in um

    % Accumulate table columns
    expCol   = [expCol,   repmat(ex,    1, nGood)];
    unitCol  = [unitCol,  1:nGood                ];
    depthCol = [depthCol, neuronDepths           ];

    % Store per-experiment data
    result.perExp(ei).exNum  = ex;
    result.perExp(ei).goodU  = goodU;
    result.perExp(ei).p_sort = p_sort;

    fprintf('  coor_Z = %.0f um | Good units: %d | Depth range: %.0f - %.0f um\n', ...
        coor_Z, nGood, min(neuronDepths), max(neuronDepths));

end

% ------------------------------------------------------------------
% Build table
% ------------------------------------------------------------------
result.depthTable = table(expCol(:), unitCol(:), depthCol(:), ...
    'VariableNames', {'Experiment', 'Unit', 'Depth_um'});

% ------------------------------------------------------------------
% Global depth bins
% ------------------------------------------------------------------
dMin = min(depthCol);
dMax = max(depthCol);
step = (dMax - dMin) / 3;

result.depthBinEdges = [dMin, dMin+step, dMin+2*step, dMax];

fprintf('\nGlobal depth range: %.0f - %.0f um\n', dMin, dMax);
fprintf('Depth bins:\n');
fprintf('  Bin 1 (shallow) : %.0f - %.0f um\n', result.depthBinEdges(1), result.depthBinEdges(2));
fprintf('  Bin 2 (middle)  : %.0f - %.0f um\n', result.depthBinEdges(2), result.depthBinEdges(3));
fprintf('  Bin 3 (deep)    : %.0f - %.0f um\n', result.depthBinEdges(3), result.depthBinEdges(4));

% ------------------------------------------------------------------
% Save to disk
% ------------------------------------------------------------------

n = extractBefore(obj.getAnalysisFileName,'lizards');
saveName = [n 'lizards' filesep 'Combined_lizard_analysis' filesep 'NeuronDepths.mat'];
save(saveName, '-struct', 'result');
fprintf('\nSaved to: %s\n', saveName);
end