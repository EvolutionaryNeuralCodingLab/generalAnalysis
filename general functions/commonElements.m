function [Common,UnCommon]=commonElements(Vec1,Vec2,MaxDiff)
%[Common,UnCommon]=commonElements(Vec1,Vec2,MaxDiff)
%finds element places in vector 1 which values are common to vector 2
%different from intersect(A, B) in that the function returns all the indices and not only first appearances
Vec2=sort(Vec2);
Common=[];
for i=1:length(Vec2)
    Common=[Common find(Vec1>=(Vec2(i)-MaxDiff) & Vec1<=(Vec2(i)+MaxDiff))];
end
Common=unique(Common);
UnCommon=1:length(Vec1);
UnCommon(Common)=[];