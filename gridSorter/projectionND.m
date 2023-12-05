function p=projectionND(v,d)
%Calculate projection between a vector and a set of dots in multi dimensional space
%v = [1 x N] - vector
%d = [M X N] - M dot locations

%calculate the cos angle between vector and dots
cosAng=v*d'./(sqrt(sum(v.^2))*sqrt(sum(d'.^2)));
p=cosAng.*sqrt(sum(d'.^2));

