function bt = hierBootEffect(vals, nrep, lvl1, lvl2)
% hierBootEffect  Hierarchical bootstrap of two effect-size statistics.
%
%   bt = hierBootEffect(vals, nrep, lvl1, lvl2)
%
%   Uses the SAME resampling scheme as hierBootMatchFreq (Saravanan et al.
%   2020): resample level-1 units (animals) with replacement, then level-2
%   units (insertions) with replacement within each drawn level-1 unit, then
%   observations with replacement within each drawn level-2 unit, pooling every
%   draw into one vector. The only difference is that hierBootMatchFreq
%   hard-codes mean(tmp) as the statistic; this returns Cohen's d and P(>0)
%   instead. Running it with mean() reproduces hierBootMatchFreq's distribution
%   to Monte Carlo error, so the two agree on the resampling itself.
%
%   lvl2 must be unique ACROSS lvl1 units (use findgroups(animal, insertion));
%   this is asserted, because a shared id would merge two animals into one
%   level-2 pool and destroy the between-animal independence.
%
%   Returns an nrep x 2 matrix: column 1 = mean/std, column 2 = fraction > 0.
%
%   See also: computeSwarmStats, hierBootMatchFreq

grp1   = unique(lvl1);                                 % level-1 (animal) ids
ngrps1 = numel(grp1);

% Pre-index the nesting once: pools{k}{m} holds the observations of the m-th
% level-2 unit belonging to the k-th level-1 unit. Selecting on lvl1 AND lvl2
% (not lvl2 alone) means a level-2 label reused across animals can never merge
% two animals' data.
pools    = cell(ngrps1, 1);
maxTotal = 0;                                          % largest possible resample size
for k = 1:ngrps1
    inThisGrp = (lvl1 == grp1(k));                     % rows of this level-1 unit
    grp2      = unique(lvl2(inThisGrp));               % its level-2 units
    pools{k}  = arrayfun(@(g) vals(inThisGrp & lvl2 == g), grp2, 'UniformOutput', false);

    % Every level-2 id must belong to exactly one level-1 unit. A shared id
    % means insertion numbering restarted per animal and two animals' data
    % would be merged; that is the documented pitfall, so fail loudly.
    assert(~any(ismember(grp2, lvl2(~inThisGrp))), ...
        ['hierBootEffect: level-2 id shared across level-1 units; ids must be ', ...
         'globally unique (use findgroups(animal, insertion)).']);

    maxTotal = max(maxTotal, sum(cellfun(@numel, pools{k})));  % obs in the largest level-1 unit
end

% A resample draws ngrps1 level-1 units, each contributing at most maxTotal
% observations, so one buffer of this size can hold any replicate.
buf = zeros(ngrps1 * maxTotal, 1);
bt  = nan(nrep, 2);                                    % [Cohen's d, P(>0)] per replicate

for j = 1:nrep
    res1 = randi(ngrps1, ngrps1, 1);                   % level-1 units, with replacement
    ptr  = 0;                                          % write position in buf

    for k = 1:ngrps1
        p    = pools{res1(k)};                         % level-2 pools of the drawn unit
        n2   = numel(p);
        res2 = randi(n2, n2, 1);                       % level-2 units, with replacement

        for m = 1:n2
            pool = p{res2(m)};                         % observations of the drawn level-2 unit
            np   = numel(pool);
            buf(ptr + (1:np)) = pool(randi(np, np, 1));% observations, with replacement
            ptr = ptr + np;
        end
    end

    tmp = buf(1:ptr);                                  % this replicate's pooled sample
    sd  = std(tmp);
    if sd > 0
        bt(j,1) = mean(tmp) / sd;                      % Cohen's d of the replicate
    end                                                % else leave NaN (degenerate replicate)
    bt(j,2) = mean(tmp > 0);                           % P(>0) of the replicate
end

end % hierBootEffect
