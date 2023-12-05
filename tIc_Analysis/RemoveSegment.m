% [tNew,icNew]=RemoveSegment(t,ic,tStart,tEnd);
%
% Function purpose : removes a time segment from t and sort channel
%
% Function recives :    t  - firing timings
%                       ic - index channel
%                       tStart - begining of segment time
%                       tEnd - ending of segment time
%                       
% Function give back :  tNew - new firing timings
%                       icNew - new index channel
% Last updated : 28/06/05

function [tNew,icNew]=RemoveSegment(t,ic,tStart,tEnd);

icNew=zeros(size(ic));
icNew(1:2,:)=ic(1:2,:);
tNew=[];
for i=1:length(ic),
    tChannel=t(ic(3,i):ic(4,i));
    tChannel(find(tChannel>tStart & tChannel<tEnd))=[];
    if isempty(tChannel),
        tChannel=NaN;
        fprintf('\nChannel %d is empty after removing segment',ic(1,i));
    end
    icNew(3,i)=length(tNew)+1;
    tNew=[tNew tChannel];
    icNew(4,i)=length(tNew);
end

    