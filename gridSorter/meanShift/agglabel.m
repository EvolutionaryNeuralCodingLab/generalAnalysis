function [clabel,clabelsize] = agglabel(label,maxlabel)
% AGGLABEL: aggregate entries with the same label
% This is basically a substitute for looping over find(label == i)
% statements, which is an O(nlabels*N) operation.
%
% Syntax:
%   clabel = agglabel(label)
%   [clabel,nlabel] = agglabel(label)
% where
%   label is a set of integer labels (must not be less than 1)
% and
%   clabel is a cell array, where clabel{i} = find(label == i);
%   nlabel is a vector containing the number of each label type.
%
% When the number of labels is large, this algorithm is rather more
% efficient than looping over find(label == i) statements.  However, it
% doesn't get down to O(N logN) as would be possible with better control over
% memory management.
%
% Update: this now calls a MEX file, which does yield O(N logN). In fact,
% the performance is almost as good as a single call to find(label == i).
% However, it only works with labels of type double.
  
% Copyright 2004 Timothy E. Holy
  
%fprintf('Using the matlab version')

  if ~all(label >= 1)
    warning('Some labels were less than 1 and will be discarded');
  end
  [clabel,clabelsize] = agglabel_mex(label);
  return
  
  
  
  
  % Here's the old algorithm
  llabel = length(label);
  if (nargin < 2)
    maxlabel = max(label);
  end
  if (maxlabel < 100)
    for i = 1:maxlabel
      clabel{i} = find(label == i);
    end
  else
    m = sparse(1:llabel,label,ones(1,llabel),llabel,maxlabel);
    %  m = sparse(label,1:llabel,ones(1,llabel),maxlabel,length(label));
    for i = 1:maxlabel
      clabel{i} = find(m(:,i))';
      %    clabel{i} = find(m(i,:));
    end
  end
  if (nargout > 1)
    for i = 1:length(clabel)
      clabelsize(i) = length(clabel{i});
    end
  end