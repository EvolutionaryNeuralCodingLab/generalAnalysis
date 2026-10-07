function [pVal, diffs, animLabels, insLabels] = bootstrapFractionDifference(tbl, pair, nBoot)
% bootstrapFractionDifference  Hierarchical bootstrap test on the per-insertion
%   difference in fraction-responsive units between two stimuli.
%
%   [pVal, diffs, animLabels, insLabels] = bootstrapFractionDifference(tbl, pair, nBoot)
%
%   tbl   - table with columns: insertion, stimulus, animal, respNeur,
%           totalSomaticN. One row per (insertion, stimulus).
%   pair  - 1x2 cell of stimulus labels; the difference is pair{1} - pair{2}.
%   nBoot - bootstrap replicates.
%
%   Returns the two-tailed p-value plus the underlying per-insertion
%   differences and their grouping labels, so callers can reuse them.
%
%   IMPORTANT — what this bootstrap can and cannot weight:
%   Each insertion contributes EXACTLY ONE value (its paired fraction
%   difference), because the units are collapsed into a ratio before the
%   bootstrap sees them. The level-3 (within-insertion) resampling is
%   therefore vacuous, and every insertion carries equal weight regardless of
%   how many somatic units it held. An insertion with 2 units — whose fraction
%   can only be 0, 0.5 or 1 — counts as much as one with 113. That
%   quantisation noise inflates the SD and attenuates the effect size; it
%   cannot be fixed by the bootstrap. See the minSomaticN sensitivity analysis
%   in AllExpAnalysis. This contrasts with the unit-level comparisons
%   (bootstrapPairDifference), where an insertion contributes one value per
%   neuron and is therefore weighted by its neuron count.
%
%   Based on Saravanan et al. (2020) via hierBootMatchFreq.

% Accumulators: one entry per insertion holding both stimuli of the pair
diffs      = [];    % per-insertion fraction difference
insLabels  = [];    % insertion ids -> hierBootMatchFreq level 2 (finest)
animLabels = [];    % animal ids    -> hierBootMatchFreq level 1 (coarsest)

for ins = unique(tbl.insertion)'

    % Rows for this insertion under each stimulus of the pair
    idx1 = tbl.insertion == ins & tbl.stimulus == pair{1};
    idx2 = tbl.insertion == ins & tbl.stimulus == pair{2};

    % Both stimuli must be present for a valid paired comparison
    if ~any(idx1) || ~any(idx2)
        continue
    end

    total = tbl.totalSomaticN(idx1);            % somatic units recorded (shared denominator)
    f1    = tbl.respNeur(idx1) / total;         % responsive fraction, stimulus 1
    f2    = tbl.respNeur(idx2) / total;         % responsive fraction, stimulus 2

    diffs      = [diffs;      f1 - f2];                        %#ok<AGROW>
    insLabels  = [insLabels;  double(ins)];                    %#ok<AGROW>
    animLabels = [animLabels; double(tbl.animal(idx1))];       %#ok<AGROW>
end

% Nothing to test (no insertion carried both stimuli)
if numel(diffs) < 2
    warning('bootstrapFractionDifference:tooFewInsertions', ...
        'Only %d insertion(s) hold both %s and %s; returning NaN.', ...
        numel(diffs), pair{1}, pair{2});
    pVal = NaN;
    return
end

% Level-2 key must be unique ACROSS animals. hierBootMatchFreq selects its
% level-2 pool with data(lvl2==id) -- keyed on lvl2 alone, not on the
% (animal, insertion) pair -- so per-animal insertion numbering would merge two
% animals into one unit and destroy the between-animal independence the
% hierarchy exists to preserve.
insKey = findgroups(animLabels, insLabels);     % composite within-animal insertion id

% Hierarchical bootstrap: resample animals -> insertions -> fractions
bootDiff = hierBootMatchFreq(diffs, nBoot, animLabels, insKey);

% Two-tailed p-value: double the smaller tail, capped at 1
pVal = min(2 * min(mean(bootDiff <= 0), mean(bootDiff >= 0)), 1);

end
