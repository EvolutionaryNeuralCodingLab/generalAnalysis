function row = emptyStatsRow(groupName)
% emptyStatsRow  One all-NaN summary row for a group or pair with no usable data.
%
%   row = emptyStatsRow(groupName)
%
%   Keeps the one-row-per-group and one-row-per-pair invariants that callers
%   assert on: a group that is absent, or a pair for which no unit could be
%   matched, still occupies its slot so row index continues to correspond to
%   x-position (and, in plotAllPairDiffs, to the entry of `pairs`).
%
%   n is set to 0 and every statistic to NaN.
%
%   See also: statsSchema, computeSwarmStats

row  = statsSchema();                                  % 0-row typed table
vals = num2cell(NaN(1, width(row) - 4));               % all statistic columns unknown
row  = [row; [{string(groupName), 0, 0, 0}, vals]];    % single placeholder row (n = 0)

end
