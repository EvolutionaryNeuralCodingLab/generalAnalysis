% [t_clean,ic_clean]=ArtifactRemove(t,ic,Fraction,ArtifWidth);
% Function purpose : removes artifact of 1ms synchoronization in many channels
%
% Recives :    t [1/12 ms only] - firing timings
%              ic - index chanel
%              Fraction - the minmal fraction of neurons with artifact to remove artifact (best with ~0.3)
%               ArtifWidth [1/12 ms]-the maximal delay of artifact between different neurons. 
%                                                                                                                                                        
% Function give back :  t_clean [1/12 ms only] - firing timings
%                       ic_clean - index chanel
% Recomended usage  : [t_clean,ic_clean]=ArtifactRemove(t,ic,0.3,3);
% Last updated : 07/01/09
function [t_clean,ic_clean]=ArtifactRemove(t,ic,Fraction,ArtifWidth);

%move to resolution where spikes are alligned.
t=round(t./ArtifWidth);
NConsecutiveSpike_Threshold=round(size(ic,2)*Fraction);
t_clean=t;ic_clean=ic;
SortedT=sort(t);
[UniqueTimes LastAppearace n]=unique(SortedT);
ToDelete=UniqueTimes(logical([0 diff(LastAppearace)>NConsecutiveSpike_Threshold]));

t_clean=t_clean*ArtifWidth;
ToDelete=ToDelete*ArtifWidth;
for i=1:length(ToDelete)
    [t_clean,ic_clean]=RemoveSegment(t_clean,ic_clean,ToDelete(i)-(6/ArtifWidth),ToDelete(i)+(6/ArtifWidth));
end
