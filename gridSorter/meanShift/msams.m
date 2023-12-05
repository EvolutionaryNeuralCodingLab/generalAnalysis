function [clust,clustinfo] = msams(varargin)
%NOTE OLIVIER: x has been replaced by varargin{1} everywhere to spare some
%memory.

% MSAMS: clustering by Minimally-Significant Adaptive Mean Shift
%
% MSAMS is a variant of a mean shift clustering algorithm that adaptively
% chooses the size of the region to average over. It also has a number of
% implementation details that make it quite efficient.
%
% Syntax:
%
%   clust = msams(x)
%
% Here x is a d-by-N matrix of data points, where each column of x is one
% point in d dimensions (N points total). On output, clust will be a
% 1-by-N vector giving the cluster # assigned to each point.
%
%   clust = msams(x,lminfo)
%
% Perform MSAMS using landmark-based optimizations. lminfo is a structure
% of the type returned by CHOOSE_LANDMARKS. As a rough guideline, you might
% want to choose the # of landmarks as sqrt(N). This syntax allows MSAMS to run
% somewhat faster, because only the landmark points are used; however, the
% cost is somewhat lower "resolution" on your clustering. In this case,
% clust is assigned for the landmark points; you can obtain assignments for
% all the points in x from clust(lminfo.landmarkAssignment).
%
%   clust = msams(x,y)
%
% An alternative syntax in which clustering is performed for "probe points"
% y, moving on a fixed background of "data points" x. This is useful in
% cases where you don't want to use landmark points as probe points.
%
%   [clust,clustinfo] = msams(...)
%
% This provides extra diagnostic information in the structure clustinfo.
% The precise form of this structure depends on the version of the
% algorithm being run (MEX or non-MEX). The most detailed description can
% be found in the help for MSAMS_CONVERGE_MEX.
%
%   clust = msams(...,options)
%
% This syntax allows you to control the algorithm with some options. See
% MSAMS_CONVERGE_MEX for the details in the "fast" case (when you
% supply lminfo), see MSAMS_CONVERGE and MSAMS_STEP1 for other cases. The
% field "compute_lminfo_if_missing" (default true) allows you to control
% whether the lminfo is calculated for calls like msams(x,y).
%
% See also: MSAMS_CONVERGE_MEX, CHOOSE_LANDMARKS, MSAMS_RESAMPLE, MSAMS_CONVERGE, MSAMS_STEP1, VMAMS_STEP1.

% Copyright 2006-2008 by Timothy E. Holy


if (length(varargin) < 1)
    error('Must supply at least one input argument');
end
%  x = varargin{1};
N = size(varargin{1},2);
have_lminfo = true;
probePointsAreDataPoints = false;
options = struct;
if (length(varargin) == 1)
    probePointsAreDataPoints = true;
    lminfo = choose_landmarks(varargin{1},ceil(sqrt(N)));
else
    if isnumeric(varargin{2})
        y = varargin{2};
        have_lminfo = false;
        if (length(varargin) > 2)
            options = varargin{3};
        end
    elseif isstruct(varargin{2})
        % is this a landmark structure or an options structure?
        if isfield(varargin{2},'landmarks')
            lminfo = varargin{2};
            if (length(varargin) > 2)
                options = varargin{3};
            end
        else
            lminfo = choose_landmarks(varargin{1},N);
            options = varargin{2};
        end
    else
        error('Input argument #2 not recognized');
    end
end

if ~isfield(options,'leapfrog')
    options.leapfrog = true;
end
if ~isfield(options,'consolidate')
    options.consolidate = true;  % change this to true once debugged
end
if ~isfield(options,'use_mex')
    options.use_mex = (exist('msams_converge_mex','file') == 3);
end
if ~isfield(options,'compute_lminfo_if_missing')
    options.compute_lminfo_if_missing = true;
end

if have_lminfo
    if probePointsAreDataPoints
        y = varargin{1};
    else
        y = lminfo.landmarks;
    end
elseif options.compute_lminfo_if_missing
    lminfo = choose_landmarks(varargin{1},size(y,2),struct('seed_landmarks',y));
    have_lminfo = true;
end

%
% Do the MSAMS steps
%


if options.use_mex && have_lminfo
    if probePointsAreDataPoints
        clustinfo = msams_converge_mex(varargin{1},lminfo,[],options);
    else
        time_passing = now;
        clustinfo = msams_converge_mex(varargin{1},lminfo,y,options);
    end
    yf = clustinfo.yf;
    map = clustinfo.closestLandmark;
    n = clustinfo.n;
else
    [yf,map,n,clustinfo] = msams_converge(varargin{1},y,options);
    clustinfo.yf = yf;
end
clustinfo.y0 = y;
clustinfo.map0 = map;
clustinfo.n = n;
% Wait to set variable_metric here so that msams_converge can handle
% defaults correctly
options.variable_metric = isfield(clustinfo,'scale');

%
% Aggregate points
%
if options.leapfrog
    % Do leapfrog convergence
    time_passing = now;
    
    map = leapfrog(map,n);
    %fprintf('Done in %.3f s\n',24*3600*(now-time_passing));
else
    % We have to aggregate the points some other way
    error('Aggregation by means other than leapfrog not yet implemented');
end
clustinfo.map1 = map;

%
% Now each cluster is defined by landmarks that map to themselves; these
% should correspond in some sense to the peaks (or at least flattest
% regions) of density. We need to check these to make sure they're not
% just some "local cycle" by testing to see whether a plurality of the
% points within the local neighborhood are a member of this cluster.
% Essentially, each point identified as belonging to the local
% neighborhood of a self-mapping landmark will "vote" for the assignment
% of the cluster overall.
%
time_passing = now;

if (options.consolidate && have_lminfo)
    mapOld = [];
    clustinfo.n_consolidated = 0;
    % For all that start out mapping to self, determine the list of
    % neighbors. We only need to do this once.
    maps_to_self0 = find(map == 1:length(map));
    maps_to_self = maps_to_self0;
    neighbor_list = cell(1,length(maps_to_self));
    for clusterIndex = 1:length(maps_to_self)
        tLandmark = maps_to_self(clusterIndex);
        sd = sqrdist(varargin{1},yf(:,tLandmark));
        [sd,neighbor_list_tmp] = sort(sd);
        neighbor_list{clusterIndex} = neighbor_list_tmp(1:n(tLandmark));
    end
    while ~isequal(mapOld,map)
        mapOld = map;
        n_maps_to_self = length(maps_to_self);
        votes = zeros(n_maps_to_self,n_maps_to_self);
        index0 = findainb(maps_to_self,maps_to_self0);
        for clusterIndex = 1:length(maps_to_self)
            % Recreate the local set of points, as the list "indexContrib"
            % We have to recalculate this (rather than saving it from before),
            % because cycling might mean that the last iteration was not the
            % one kept, and it's too memory intensive to store the full list on
            % each iteration.
            tLandmark = maps_to_self(clusterIndex);
            indexContrib = neighbor_list{index0(clusterIndex)};
            % Of the points in the neighborhood, what is their cluster
            % assignment?
            if probePointsAreDataPoints
                nbrmap = map(indexContrib);
            else
                nbrmap = map(lminfo.landmarkAssignment(indexContrib));
            end
            % Calculate the number of points assigned to each cluster within the
            % local neighborhood---this tallies the "votes"
            [unbrmap,tmp,relabelled_nbrmap] = unique(nbrmap);
            [tmp,votes_tmp] = agglabel(relabelled_nbrmap);
            votes(clusterIndex,findainb(unbrmap,maps_to_self)) = votes_tmp;
        end
        % Re-assign each self-mapping landmark to the one that receives the
        % most "votes" in its local neighborhood (note we do this only after
        % all votes have been tallied for all "candidates," so that there can't
        % be any order effects.  Not that those seem likely...)
        [max_votes,maxIndex] = max(votes,[],2);
        map(maps_to_self) = maps_to_self(maxIndex);
        clustinfo.n_consolidated = clustinfo.n_consolidated + sum(maxIndex' ~= 1:n_maps_to_self);
        % Iterate map to get final assignments
        map = leapfrog(map);
        maps_to_self = find(map == 1:length(map));
    end
end

clustinfo.map = map;
% Re-label the clusters 1,2,3,...
[umap,tmp,clust] = unique(map);
%fprintf('Done in %.3f s\n',24*3600*(now-time_passing));
