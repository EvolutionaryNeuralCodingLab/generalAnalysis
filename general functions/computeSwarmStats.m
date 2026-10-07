function statsTbl = computeSwarmStats(tblPlot, stimuli, params)
% computeSwarmStats  Per-group descriptives + hierarchical-bootstrap intervals.
%
%   statsTbl = computeSwarmStats(tblPlot, stimuli, params)
%
%   The canonical statistics routine behind every swarm panel. Returns exactly
%   one row per entry of `stimuli`, in that same order, so the row index doubles
%   as the x-position of the corresponding swarm column. Groups that are absent,
%   or have n < 3, still occupy a row (see emptyStatsRow).
%
%   Call this directly to produce a schema-compatible stats block WITHOUT
%   plotting — that is how analyses whose bootstraps do not route through
%   plotSwarmBootstrapWithComparisons (e.g. the depth and stripe-direction
%   comparisons in analyzeStripeNeurons) get rows in the supplementary table.
%
%   INPUTS
%     tblPlot - table with a `stimulus` column, a `value` column, and the
%               columns named in params.bootGroupVars (normally animal,
%               insertion).
%     stimuli - cellstr of group labels, in the order rows should appear.
%     params  - struct with fields:
%                 nBoot         bootstrap replicates
%                 ciMethod      'percentile' | 'sem'
%                 bootGroupVars cell, normally {'animal','insertion'};
%                               {} forces a flat (non-hierarchical) bootstrap
%
%   The point estimate drawn as the bar is the OBSERVED mean, not bootMean;
%   see statsSchema for what each column means and why.
%
%   See also: statsSchema, emptyStatsRow, hierBootEffect, hierBootMatchFreq,
%             plotSwarmBootstrapWithComparisons, buildStatsTable

statsTbl = statsSchema();                              % typed accumulator (stays 0-row if stimuli is empty)

for i = 1:numel(stimuli)
    idx = tblPlot.stimulus == stimuli{i};

    % Absent group: emit a placeholder so row i still maps to category i
    if ~any(idx)
        statsTbl = [statsTbl; emptyStatsRow(stimuli{i})];   %#ok<AGROW>
        continue
    end

    % Row indices for this group, with NaN values dropped. Indexing by row
    % number (rather than two successive logical masks) keeps the value vector
    % and every grouping column trivially aligned.
    rows = find(idx);                                  % all rows of this group
    vals = tblPlot.value(rows);                        % their values
    keep = ~isnan(vals);                               % NaN mask
    rows = rows(keep);                                 % rows that survive
    vals = vals(keep);                                 % values that survive

    % Descriptive counts: how many animals/insertions actually contribute
    if ismember('animal', tblPlot.Properties.VariableNames)
        nAnimals = numel(unique(tblPlot.animal(rows)));    % distinct animals in this group
    else
        nAnimals = NaN;                                    % column absent
    end
    if ismember('insertion', tblPlot.Properties.VariableNames)
        nInsertions = numel(unique(tblPlot.insertion(rows)));  % distinct insertions
    else
        nInsertions = NaN;                                 % column absent
    end

    % Plain pooled descriptives over the observations actually plotted.
    n         = numel(vals);
    rawMean   = mean(vals);
    rawMedian = median(vals);

    % Observed effect sizes, both computed on the plotted values so they carry
    % the meaning of whichever table this is:
    %   groups -> one-sample Cohen's d against zero, and P(value > 0)
    %   diffs  -> Cohen's d_z (mean of paired differences / their SD), and the
    %             common-language effect size P(stimA > stimB)
    rawSD = std(vals);                                 % SD of the plotted values
    if rawSD > 0
        rawEffect = rawMean / rawSD;                   % Cohen's d / d_z
    else
        rawEffect = NaN;                               % degenerate (all values identical)
    end
    rawPropPos = mean(vals > 0);                       % P(>0); Cliff's delta = 2*propPos - 1 when there are no ties

    % Too few points for a meaningful bootstrap: keep the observed statistics,
    % leave every CI column NaN, and draw no bar.
    if n < 3
        fprintf('Stimulus %s: n=%d < 3; skipping mean/SE bar.\n', char(stimuli{i}), n);
        statsTbl = [statsTbl; {string(stimuli{i}), n, nAnimals, nInsertions, ...
            rawMean, rawMedian, NaN, NaN, NaN, ...
            rawEffect, NaN, NaN, rawPropPos, NaN, NaN}];   %#ok<AGROW>
        continue
    end

    % Pull each grouping column aligned with NaN drop
    groupVars = params.bootGroupVars;
    groupVals = cell(1, numel(groupVars));
    for g = 1:numel(groupVars)
        col = tblPlot.(groupVars{g})(rows);
        if iscategorical(col)
            col = double(col);                         % hierBoot requires numeric
        end
        groupVals{g} = col;
    end

    % Make the level-2 (insertion) key unique *within* animal. Insertion
    % numbering is reset per animal, so double(insertion) gives the same index
    % to "insertion 1" of every animal; findgroups over (animal, insertion)
    % yields one ID per distinct pair, so two animals' same-labelled insertions
    % can never be merged when the bootstrap builds its animal->insertion
    % nesting. No-op if the resampler already enumerates insertions within each
    % animal; a silent-bug fix if it keys insertion globally.
    if numel(groupVars) == 2 && all(strcmp(groupVars(:)', {'animal','insertion'}))
        groupVals{2} = findgroups(groupVals{1}, groupVals{2});   % composite within-animal insertion ID
    end


    % Hierarchical (animal -> insertion) or flat bootstrap of the mean.
    % hierBootMatchFreq is hard-coded for exactly two grouping levels and
    % expects them coarsest-first (lvl1 = animal, lvl2 = insertion). The
    % default bootGroupVars ordering {'animal','insertion'} satisfies this.
    if isempty(groupVars)
        bootMean   = bootstrp(params.nBoot, @mean, vals);        % no hierarchy: ordinary bootstrap
        bootEffect = bootstrp(params.nBoot, ...                  % same resamples, effect statistics
            @(x) [mean(x)/std(x), mean(x > 0)], vals);
    else
        assert(numel(groupVars) == 2, ...                        % enforce the 2-level contract
            ['computeSwarmStats: hierBootMatchFreq needs exactly two grouping ', ...
            'levels (animal, insertion); got %d. Adjust bootGroupVars, or ', ...
            'pass {} for a flat bootstrap.'], numel(groupVars));
        bootMean = hierBootMatchFreq(vals, params.nBoot, ...     % lvl1 = animal (coarsest)
            groupVals{1}, groupVals{2});                         % lvl2 = insertion
        % Effect-size CIs need a statistic other than the mean, which
        % hierBootMatchFreq cannot return, so a second pass over the SAME
        % resampling scheme is run here (see hierBootEffect).
        bootEffect = hierBootEffect(vals, params.nBoot, ...
            groupVals{1}, groupVals{2});
    end

    % Mean of the bootstrap distribution. This is NOT the plotted point
    % estimate — it is kept as a diagnostic. It is also NOT an equal-weight-
    % per-animal mean: hierBootMatchFreq pools every resampled observation and
    % takes one flat mean, so the resampled statistic is a ratio with a random
    % denominator (animals contribute different numbers of neurons). Its
    % expectation therefore drifts away from the observed pooled mean whenever
    % neuron counts are unbalanced across animals, landing between the pooled
    % mean and the mean-of-animal-means. The gap (rawMean - mu) is exactly the
    % bootstrap's bias estimate for this statistic.
    mu = mean(bootMean);

    % Uncertainty bounds. Following the standard bootstrap convention (Efron &
    % Tibshirani), the POINT ESTIMATE is the observed statistic (rawMean, drawn
    % as the bar) and the resampling supplies only the interval.
    switch lower(params.ciMethod)
        case 'sem'
            % "estimate +/- bootstrap SE": the interval must be centred on the
            % same value the bar draws, or the bar could fall outside its own
            % error bar when the bootstrap is biased.
            se  = std(bootMean);                        % bootstrap standard error
            yLo = rawMean - se;                         % centred on the OBSERVED estimate
            yHi = rawMean + se;
        case 'percentile'
            % A percentile interval is defined by the bootstrap distribution
            % itself, so it is left as-is and may sit asymmetrically about the
            % observed estimate. That asymmetry is informative, not an error.
            yLo = prctile(bootMean,  2.5);
            yHi = prctile(bootMean, 97.5);
        otherwise
            error('Unknown ciMethod: %s. Use ''sem'' or ''percentile''.', params.ciMethod);
    end

    % Effect-size intervals are ALWAYS percentile, regardless of ciMethod: a
    % symmetric +/-SE interval is inappropriate for a bounded quantity like
    % P(>0) (it can run outside [0,1]) and for a ratio like Cohen's d.
    effLo  = prctile(bootEffect(:,1),  2.5);           % Cohen's d lower bound
    effHi  = prctile(bootEffect(:,1), 97.5);           % Cohen's d upper bound
    propLo = prctile(bootEffect(:,2),  2.5);           % P(>0) lower bound
    propHi = prctile(bootEffect(:,2), 97.5);           % P(>0) upper bound

    % Store the row: observed statistics (mean/median/effect/propPos) plus the
    % bootstrap intervals, and bootMean as a bias diagnostic.
    statsTbl = [statsTbl; {string(stimuli{i}), n, nAnimals, nInsertions, ...
        rawMean, rawMedian, mu, yLo, yHi, ...
        rawEffect, effLo, effHi, rawPropPos, propLo, propHi}];   %#ok<AGROW>
end

end % computeSwarmStats
