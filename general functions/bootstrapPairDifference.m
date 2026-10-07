function pVal = bootstrapPairDifference(tbl, pair, nBoot, metric)
% bootstrapPairDifference  Hierarchical bootstrap test for a single pair.
%
%   BUG FIX (was Bug #2): Pairs neurons explicitly by NeurID within each
%   insertion using innerjoin, instead of relying on row order.  The
%   previous row-order approach silently produced wrong differences if the
%   table was ever sorted or filtered asymmetrically.
%
%   INPUTS
%     tbl    table   Long-format with columns: insertion, stimulus, NeurID,
%                    animal, and the metric column.
%     pair   {1x2}  Cell pair of stimulus labels.
%     nBoot  double  Number of bootstrap iterations.
%     metric char    Column name to test ('Z-score' or 'SpkR').
%
%   OUTPUT
%     pVal   double  Two-tailed p-value from hierarchical bootstrap.

    diffs   = [];       % per-neuron paired differences
    insers  = [];       % insertion label for each difference (hierBoot level 1)
    animals = [];       % animal label for each difference   (hierBoot level 2)

    for ins = unique(tbl.insertion)'

        % Logical masks: rows for this insertion x each stimulus
        mask1 = tbl.insertion == ins & tbl.stimulus == pair{1};
        mask2 = tbl.insertion == ins & tbl.stimulus == pair{2};

        if ~any(mask1) || ~any(mask2), continue; end    % skip if either is absent

        % Extract sub-tables with the metric, NeurID, and animal columns
        sub1 = tbl(mask1, {metric, 'NeurID', 'animal'});
        sub2 = tbl(mask2, {metric, 'NeurID'});

        % Inner-join on NeurID ensures correct neuron-to-neuron pairing
        merged = innerjoin(sub1, sub2, 'Keys', 'NeurID', ...
                           'LeftVariables',  {metric, 'NeurID', 'animal'}, ...
                           'RightVariables', {metric});

        % MATLAB auto-suffixes duplicated variable names after innerjoin;
        % find both metric columns by prefix matching
        mergedVarNames = merged.Properties.VariableNames;
        metricCols     = mergedVarNames(startsWith(mergedVarNames, metric));

        % Paired difference: stimulus 1 minus stimulus 2
        d = merged.(metricCols{1}) - merged.(metricCols{2});

        % Animal identity (constant within an insertion)
        animal = merged.animal(1);

        % Append to accumulators
        nPaired = numel(d);
        diffs   = [diffs;   d];                                             %#ok<AGROW>
        insers  = [insers;  double(repmat(ins,    nPaired, 1))];            %#ok<AGROW>
        animals = [animals; double(repmat(animal, nPaired, 1))];            %#ok<AGROW>
    end

    % Level-2 key must be unique ACROSS animals. hierBootMatchFreq selects its
    % level-2 pool with data(lvl2==id) -- keyed on lvl2 alone, not on the
    % (animal, insertion) pair -- so per-animal insertion numbering would merge
    % two animals into one unit and destroy the between-animal independence the
    % hierarchy exists to preserve. findgroups gives one id per distinct
    % (animal, insertion) combination. Matches computeSwarmStats and the
    % fraction-responsive bootstrap in AllExpAnalysis.
    insKey = findgroups(animals, insers);   % composite within-animal insertion id

    % Hierarchical bootstrap: resample animals -> insertions -> neurons
    bootMeans = hierBootMatchFreq(diffs, nBoot, animals, insKey);
    

    % Two-tailed p-value: probability that |bootstrap mean| >= |observed|
    pLeft  = mean(bootMeans <= 0);          % fraction of bootstrap means <= 0
    pRight = mean(bootMeans >= 0);          % fraction of bootstrap means >= 0
    pVal   = 2 * min(pLeft, pRight);        % double the smaller tail
    pVal   = min(pVal, 1);                  % cap at 1
end
