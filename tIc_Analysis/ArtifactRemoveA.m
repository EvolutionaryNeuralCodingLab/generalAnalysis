% [t_clean,ic_clean]=ArtifactRemove(t,ic,Fraction,ArtifWidth);
% Function purpose : removes artifact of 1ms synchoronization in many channels
%
% Recives :    I - activity intensity
%                           t [ms only] - firing timings
%                           ic - index chanel
%                           Fraction - the minmal fraction of neurons with artifact to remove artifact (best with ~0.3)
%                           ArtifWidth [ms]-the maximal delay of artifact between different neurons. 
%                                                                                                                                                        
% Function give back :   I_clean - activity intensity
%                                                   t_clean [ms only] - firing timings
%                                                   ic_clean - index chanel
% Recomended usage  : [I_clean,t_clean,ic_clean]=ArtifactRemove(I,t,ic,0.3,0.2);
% Last updated : 01/08/09
function [I_clean, t_clean,ic_clean]=ArtifactRemoveA(I,t,ic,Fraction,ArtifWidth);

%move to resolution where spikes are alligned.
t=round(t./ArtifWidth);
NConsecutiveSpike_Threshold=round(size(ic,2)*Fraction);
I_clean=I;t_clean=t;ic_clean=ic;
SortedT=sort(t);
[UniqueTimes LastAppearace n]=unique(SortedT);
ToDelete=UniqueTimes(logical([0 diff(LastAppearace)>NConsecutiveSpike_Threshold]));

t_clean=t_clean*ArtifWidth;
ToDelete=ToDelete*ArtifWidth;
for i=1:length(ToDelete)
    [I_clean,t_clean,ic_clean]=RemoveSegmentA(I_clean,t_clean,ic_clean,ToDelete(i)-(6/ArtifWidth),ToDelete(i)+(6/ArtifWidth));
end
