% [INew,tNew,icNew]=RemoveSegmentA(I,t,ic,tStart,tEnd);
%
% Function purpose : removes a time segment from t and sort channel
%
% Function recives :    I - activity intensity
%                                                 t  - firing timings
%                                                 ic - index channel
%                                                  tStart - begining of segment time
%                                                  tEnd - ending of segment time
%                       
% Function give back :  INew - activity intensity
%                                                   tNew - new firing timings
%                                                   icNew - new index channel
% Last updated : 1/07/09

function [INew,tNew,icNew]=RemoveSegmentA(I,t,ic,tStart,tEnd)

icNew=zeros(size(ic));
icNew(1:2,:)=ic(1:2,:);
tNew=[];
INew=[];
for i=1:length(ic),
    tChannel=t(ic(3,i):ic(4,i));
    IChannel=I(ic(3,i):ic(4,i));
    Places=find(tChannel>tStart & tChannel<tEnd);
    tChannel(Places)=[];
    IChannel(Places)=[];
    if isempty(tChannel),
        tChannel=NaN;
        IChannel=NaN;
        fprintf('\nChannel %d is empty after removing segment',ic(1,i));
    end
    icNew(3,i)=length(tNew)+1;
    tNew=[tNew tChannel];
    INew=[INew IChannel];
    icNew(4,i)=length(tNew);
end