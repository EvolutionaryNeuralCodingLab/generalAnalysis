function [M]=normZeroOne(M,scale)
%[M]=normZeroOne(M,scale)
%normalize each column between 0-1 (or scale) separately
if nargin==1
    scale=[0 1];
end
S=size(M);
if S(1)==1 || S(2)==1
    M=(M-min(M))./(max(M)-min(M));
elseif numel(S)==2
    M=bsxfun(@rdivide,bsxfun(@minus,M,min(M)),max(M)-min(M));
end
M=M.*(scale(2)-scale(1))+scale(1);
