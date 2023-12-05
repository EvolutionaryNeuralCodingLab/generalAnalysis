% [ic_sep]=Separate_networks(ic);
function [ic_sep]=Separate_networks(ic);

side_a=[1 2 3 7 8 9 10 15 16 17 18 23 24 25 26 31 32 33 34 39 40 41 42 47 48 49 50 55 56 57];
side_b=[4 5 6 11 12 13 14 19 20 21 22 27 28 29 30 35 36 37 38 43 44 45 46 51 52 53 54 58 59 60];

place=[];
for i=1:size(ic,2)
    if size(find(side_a==ic(1,i)),2)~=0
        place=[place 1];
    else
        place=[place 2];
    end
end
[sorted,order]=sort(place);
ic_sep=ic(:,order);
