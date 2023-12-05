function lminfo = choose_landmarks(x,k,options)
% choose_landmarks: split a group of points up into local regions
% Syntax:
%   lminfo = choose_landmarks(x,k)
%   lminfo = choose_landmarks(x,k,options)
% where
%   x is a d-by-N matrix of data points, each column a separate data
%     point in d dimensions;
%   k is the number of landmarks desired;
%   options is a structure with the following fields:
%     frac_random (default 0): the fraction of landmarks chosen randomly
%       from the data points. The remainder of landmarks will be chosen
%       iteratively as the point farthest from the existing landmarks.
%       A value of 0 will choose the first landmark randomly, and the rest
%       according to distance. frac_random = 0 will tend to make landmark
%       groupings of very similar physical size; frac_random = 1 will tend
%       to make landmark groupings containing relatively similar #s of
%       points.
%     metric (default @sqrdist): a handle to a function that computes the
%       distance between points.  This function must be able to handle a
%       matrix of data points for its first argument. 'sqrdist' is an
%       example of such a function.
%     seed_landmarks: if present, starts the process of landmark creation
%       using the supplied landmarks, and then creates the rest as needed.
%     center: if true, causes the landmarks to move inward from the "edge"
%       of their assigned group towards the mean (it's like doing a single
%       round of k-means).
%     progress (default true): if true, will provide an update about
%       progress every five seconds.
% and
%   lminfo is an output structure with the following fields:
%     landmarks: a d-by-k matrix, giving the positions of the landmarks;
%     landmarkAssignment: a 1-by-N vector, landmarkAssignment(i) is the
%       index of the landmark closest to the ith data point, x(:,i);
%     landmarkList: a 1-by-k cell array, each element giving the indices of
%       the points in x closest to this landmark;
%     dist_to_closest_landmark: an N-by-1 vector, containing the distance
%       to the closest landmark for each data point;
%     radius_of_neighborhood: a k-by-1 vector, giving the distance to the
%       farthest point assigned to each landmark;
%     landmark_xIndex: a 1-by-k vector, giving the indices of the points in
%       x used as landmarks (NaNs are used for seed_landmarks);
%     metric: the function handle used to compute distances.
%
% See also: REINDEX_LANDMARKS.
  
% Copyright 2006 by Timothy E. Holy
    
  if (nargin < 3)
    options = struct;
  end
  if ~isfield(options,'metric')
    options.metric = @sqrdist;
  end
  if ~isfield(options,'frac_random')
    options.frac_random = 0;
  end
  if ~isfield(options,'center')
    options.center = false;
  end
  if ~isfield(options,'progress')
    options.progress = true;
  end

  % See if we're using sqrdist (that will allow us to use mindist to save memory)
  funcstr = func2str(options.metric);
  if isequal(funcstr,'sqrdist')
    using_sqrdist = true;
  end
  
  [d,N] = size(x);
  
  if options.progress
    first_update = true;
    try
      % Try to avoid resetting the timer if possible
      tprev = toc;
    catch
      tic;
      tprev = 0;
    end
  end
  
  if(k>=N)
     lminfo.landmarks=x;
     lminfo.landmarkAssignment=1:N;
     lminfo.landmarkList=agglabel(1:N);
     lminfo.dist_to_closest_landmark = zeros(1,N);
     lminfo.radius_of_neighborhood = zeros(1,N);
     lminfo.landmark_xIndex = 1:N;
     lminfo.metric = options.metric;
     return
  end
  
  y = nan([d k],class(x));
  
  if isfield(options,'seed_landmarks')
    % Use the seed landmarks
    n_start = min(k,size(options.seed_landmarks,2));
    y(:,1:n_start) = options.seed_landmarks(:,1:n_start);
    xIndex(1:n_start) = nan;
  else
    % Start with random point(s)
    n_start = 1;
    if (n_start < 1)
      n_start = 1;
    end
    if (n_start > k)
      n_start = k;
    end
    if (n_start == 1)
      xIndex = round(N*rand(1,1)+0.5);
      xIndex=1;
    else
      xIndex = randperm(N);
    end
    xIndex = xIndex(1:n_start);
    y(:,1:n_start) = x(:,xIndex);
  end
  % Calculate distances to "starter" landmarks
  % Do this in a way that will work even if memory is tight
  if using_sqrdist
    if (options.progress && d*N*n_start > 1e7)
      fprintf('Choosing %d landmarks: ',k);
      first_update = false;
    end
    
    [R2x,landmarkIndex] = mindist(x,y(:,1:n_start));
  else
    try
      cdist = options.metric(x,y(:,1:n_start));
      % For each data point, calculate the distance to the closest landmark (&
      % keep track of which one it is)
      [R2x,landmarkIndex] = min(cdist,[],2);
      landmarkIndex = landmarkIndex';
    catch
      % We go here if memory is insufficient for doing them all at once
      R2x = zeros(1,N);
      landmarkIndex = zeros(1,N);
      for i = 1:N
        cdist = options.metric(x(:,i),y(:,1:n_start));
        [R2x(i),landmarkIndex(i)] = min(cdist);
      end
    end
  end
  % In case the data set contains duplicated points, we need to make sure
  % that every one of the starter landmarks has points assigned, and toss
  % any that have no points
  [xlabel,nlabel] = agglabel_mex(landmarkIndex);
  keepFlag = (nlabel > 0);
  xlabel = xlabel(keepFlag);
  n_start = sum(keepFlag);
  xIndex = xIndex(keepFlag);
  y(:,1:n_start) = x(:,xIndex);
  for i = 1:n_start
    landmarkIndex(xlabel{i}) = i;
  end
  % Now choose the rest of the landmarks iteratively by finding the point
  % that is most distance from the current landmarks
  for yIndex = n_start+1:k
    % Find the point farthest from the existing landmarks
    [maxdist,maxIndex] = max(R2x);
    xIndex(yIndex) = maxIndex;  % this is the new landmark
    % Compute the distance of each point to this new landmark
    cdist = options.metric(x(:,maxIndex),x);
    % Determine the set of points for which this landmark is the closest yet
    is_closest = (cdist < R2x);
    landmarkIndex(is_closest) = yIndex;
    R2x(is_closest) = cdist(is_closest);
    y(:,yIndex) = x(:,maxIndex);
    if options.progress
      if (toc - tprev > 5)
        tprev = toc;
        if first_update
          fprintf('Choosing %d landmarks: ',k);
          first_update = false;
        end
        fprintf('%d...',yIndex);
      end
    end
  end

  if (options.progress && ~first_update)
    fprintf('done\n');
  end

  [xlabel,nlabel] = agglabel(landmarkIndex);
  if (options.center && using_sqrdist)
    % Move each landmark to the center of the points that are currently
    % assigned to it
    for yIndex = 1:k
      y(:,yIndex) = mean(x(:,xlabel{yIndex}),2);
    end
    [R2x,landmarkIndex] = mindist(x,y);
    [xlabel,nlabel] = agglabel(landmarkIndex);
    % Toss any that have no points assigned
    keepflag = (nlabel > 0);
    y = y(:,keepflag);
    xlabel = xlabel(keepflag);
    k = sum(keepflag);
  end
  %compute radii
  for yIndex = 1:k
    R2y(yIndex) = max(R2x(xlabel{yIndex}));
  end
  
  lminfo.landmarks = y;
  lminfo.landmarkAssignment = landmarkIndex;
  lminfo.landmarkList = xlabel;
  lminfo.dist_to_closest_landmark = sqrt(R2x);
  lminfo.radius_of_neighborhood = sqrt(R2y(:));
  lminfo.landmark_xIndex = xIndex;
  lminfo.metric = options.metric;
