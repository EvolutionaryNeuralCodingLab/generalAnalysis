function [IB,tB,icB]=Bursts2Spikes(BP,BI,icA)
%[IB,tB,icB]=Bursts2Spikes(BP,BI,icA);
% Function purpose : Coverts burst timing (Burst peak) and intensity (burst intensity) cell arrays to index channel representation
%
% Function recives :    BI [ms]- cell array of the number of clusters length with the total activty intensity of every burst
%                                                 BP [ms]- cell array of the number of clusters length with the time (peak location) of every burst
%                                                 icA - the indexc channel of the recording
%
% Function give back :  IB - total activity intesities of the bursts
%                                                   tI [ms] - timings of the bursts
%                                                   icB - index channel of the bursts
%
% Last updated : 18/06/10
icB(3,1)=1;
tB=[];
IB=[];
c=0;
for i=1:size(icA,2)
    if ~isempty(BP{i})
        c=c+1;
        tB=[tB BP{i}];
        IB=[IB BI{i}];
        icB(1:2,c)=icA(1:2,i);
        icB(4,c)=icB(3,c)+length(BP{i})-1;
        icB(3,c+1)=icB(4,c)+1;
    end
end
icB(:,end)=[];
        