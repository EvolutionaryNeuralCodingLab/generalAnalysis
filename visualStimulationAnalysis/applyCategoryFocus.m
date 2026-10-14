function [keepMask, catID, catTable, sepNames] = applyCategoryFocus(C, colNames, levelFactors, lockSpec)
% applyCategoryFocus  Build effective response categories from a condition
%   matrix under a focus spec: keep named factors SEPARATE (their crossing
%   defines the levels), POOL every unnamed factor, and optionally LOCK
%   factors to a fixed value (pre-filter). Shared by StatisticsPerNeuron and
%   the PSTH/raster plotters so statistics and displays use identical trials.
%
%   INPUTS
%     C            [nTrials × (1+nParams)]  C(:,1)=onset, C(:,k+1)=parameter k
%     colNames     [1×nParams] names; colNames{k} <-> C(:,k+1)
%     levelFactors factors kept separate ("all"/{} -> every parameter separate)
%     lockSpec     {name1,val1,...} locked factors ({} -> none)
%   OUTPUTS
%     keepMask  [nTrials×1] logical  trials surviving the lock filter
%     catID     [nKept×1]            effective category (1..nCats) per kept trial
%     catTable  [nCats×nSep]         level values per category (cols = sepNames)
%     sepNames  [1×nSep] string      separating-factor names

    tol      = 1e-6;                                   % float tolerance for level matching
    nTrials  = size(C, 1);                             % total trials
    colNames = string(colNames(:))';                   % normalise to a string row
    nameToCol = @(nm) findParamCol(colNames, nm);      % name -> column in C (errors if absent)

    % --- Lock filter (pre-select trials) ----------------------------------
    keepMask = true(nTrials, 1);                        % start with all trials
    if ~isempty(lockSpec)
        assert(mod(numel(lockSpec),2)==0, ...          % must be name/value pairs
            'Lock must be name/value pairs, e.g. {''luminosities'',255}.');
        for q = 1:2:numel(lockSpec)                     % iterate pairs
            col      = nameToCol(lockSpec{q});          % column of this factor
            keepMask = keepMask & abs(C(:,col) - lockSpec{q+1}) < tol;  % keep matching trials
        end
        assert(any(keepMask), 'Lock specification matched no trials.');  % must leave something
    end

    % --- Resolve separating factors ---------------------------------------
    if isempty(levelFactors) || (isscalar(string(levelFactors)) && string(levelFactors)=="all")
        sepNames = colNames;                            % default: every parameter separate
    else
        sepNames = string(levelFactors(:))';            % requested subset
    end
    sepCols = arrayfun(nameToCol, sepNames);            % their columns in C (errors on any miss)

    % --- Group kept trials by the separating-factor crossing --------------
    sub        = C(keepMask, sepCols);                  % [nKept × nSep] levels per kept trial
    subRounded = round(sub / tol) * tol;                % snap to tolerance grid (no float-split)
    [catTable, ~, catID] = unique(subRounded, 'rows');  % unique combos -> catTable; catID in 1..nCats
    % Unnamed factors are absent from sepCols, so their trials share a catID
    % -> pooled automatically. catTable/catID are deterministic (sorted rows).
end

function col = findParamCol(colNames, nm)
% findParamCol  Case-insensitive name -> C-column lookup (hard error on miss).
    idx = find(strcmpi(colNames, string(nm)), 1);       % match parameter name
    assert(~isempty(idx), 'Factor "%s" not found. Available: %s', ...
        string(nm), strjoin(colNames, ', '));
    col = idx + 1;                                       % +1 because C(:,1) is onset
end