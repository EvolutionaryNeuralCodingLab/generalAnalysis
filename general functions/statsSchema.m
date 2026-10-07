function tblOut = statsSchema()
% statsSchema  Canonical 0-row summary table for swarm-panel statistics.
%
%   tblOut = statsSchema()
%
%   Defined once, in one file, so every producer of stats rows
%   (computeSwarmStats, emptyStatsRow, and any analysis assembling its own
%   blocks) emits identical column names and types. buildStatsTable relies on
%   this schema being stable; changing it changes the supplementary table.
%
%   Columns:
%     group        - label ('MB' for a group, 'MB-SB' for a difference)
%     n            - non-NaN observations
%     nAnimals     - distinct animals contributing
%     nInsertions  - distinct insertions contributing
%     mean, median - plain pooled sample statistics; `mean` is what the swarm
%                    bar draws
%     bootMean     - mean of the hierarchical-bootstrap distribution
%                    (DIAGNOSTIC only, not plotted; a large gap from `mean`
%                    signals unbalanced per-animal counts)
%     ciLo, ciHi   - bootstrap interval on the mean
%     effect       - Cohen's d (d_z when the values are paired differences)
%     effectLo/Hi  - 95% percentile CI on effect
%     propPos      - fraction of values > 0
%     propPosLo/Hi - 95% percentile CI on propPos
%
%   See also: computeSwarmStats, emptyStatsRow, buildStatsTable

names = {'group','n','nAnimals','nInsertions', ...       % identity and sample sizes
         'mean','median','bootMean','ciLo','ciHi', ...   % central tendency + its CI
         'effect','effectLo','effectHi', ...             % Cohen's d (d_z for differences) + CI
         'propPos','propPosLo','propPosHi'};             % P(value > 0) + CI

types      = repmat({'double'}, 1, numel(names));        % every column numeric...
types{1}   = 'string';                                   % ...except the group label

tblOut = table('Size', [0 numel(names)], ...
    'VariableTypes', types, 'VariableNames', names);

end
