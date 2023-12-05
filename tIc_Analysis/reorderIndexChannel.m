%[tNew,icNew]=reorderIndexChannel(t,ic)
% Function purpose : order channels in t,ic in increasing order
%
% Function recives :    t [ms] - firing timings
%                       ic - indexc channel
%
% Function give back :  tNew [ms] - firing timings after ordering
%                       icNew - firing timings after ordering
%
% Last updated : 27/06/12
function [tNew,icNew]=reorderIndexChannel(t,ic)

tNew=[];icNew=[];
ch=unique(ic(1,:));
c=0;
for i=1:numel(ch)
    p=find(ic(1,:)==ch(i));
    p2=unique(ic(2,p));
    for j=1:numel(p2)
        pOld=p(p2(j));
        tNew=[tNew t( ic(3,pOld):ic(4,pOld) )];
        L=ic(4,pOld)-ic(3,pOld)+1;
        icNew=[icNew [ic(1:2,pOld);c+1;c+L] ];
        c=c+L;
    end
end