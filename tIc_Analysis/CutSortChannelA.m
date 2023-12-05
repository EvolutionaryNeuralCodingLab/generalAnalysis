% [INew,tNew,icNew]=CutSortChannelA(I,t,ic,tStart,tEnd);
% Function purpose : cuts out a peace of a Activity recording between timings tStart
%                    if a neuron does not fire in this piece it is deleted
%
% Function recives :    A - activity intensity file
%                       t [ms] - firing timings
%                       ic - indexc channel
%                       tStart - starting time [ms]
%                       tEnd - ending time [ms]
%                
% Function give back :  tNew [ms]- new firing timings
%                       icNew - new indexc channel
%                       INew - new activity intensity
%
% To show paster plot use : plotraster(tNew,icNew);
% Last updated : 14/05/09

function [INew,tNew,icNew]=CutSortChannelA(I,t,ic,tStart,tEnd);

to_delete=[];
icNew=zeros(size(ic));
icNew(1:2,:)=ic(1:2,:);
tNew=[];INew=[];
for i=1:size(ic,2),
    tChannel=t(ic(3,i):ic(4,i));
    IChannel=I(ic(3,i):ic(4,i));
    P=find(tChannel>tStart & tChannel<tEnd);
    tCut=tChannel(P);
    ICut=IChannel(P);
    if isempty(tCut),
        fprintf('\nNeuron %d %d does not fire during selected interval and therfore was deleted',ic(1:2,i));
        to_delete=[to_delete i];
    else        
        icNew(3,i)=length(tNew)+1;
        tNew=[tNew, tCut];
        INew=[INew, ICut];
        icNew(4,i)=length(tNew);
    end
end
icNew(:,to_delete)=[];
fprintf('\nPlease remember that cut segment started at time %d in the whole sort Channel\n',tStart);
tNew=tNew-tStart;

    