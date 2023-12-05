% [tNew,icNew,tStart,to_delete]=CutSortChannel(t,ic,tStart,tEnd,outMessage,cutEmptyTime,cutEmptyNeurons);
% Function purpose : cuts out a peace of a recording between timings tStart
%                    if a neuron does not fire in this piece it is deleted
%
% Function recives :    t - firing timings
%                       ic - indexc channel
%                       tStart - starting time 
%                       tEnd - ending time
%                       cutEmptyTime - change times so that the first spike will be in time 1
%                       cutEmptyNeurons - remove neurons not firing in the time segment
%                       outMessage - send out messages
%                
% Function give back :  tNew - new firing timings
%                       icNew - new indexc channel
%                        tStart - the start time of the new sort channel
%                       to_delete - number of empty channels after cutting
%
% To show paster plot use : plotraster(tNew,icNew);
% Last updated : 02/03/09

function [tNew,icNew,tStart,to_delete]=CutSortChannel(t,ic,tStart,tEnd,outMessage,cutEmptyTime,cutEmptyNeurons)
if nargin<4
    error('Not enough arguments');
end
if nargin<=4
    outMessage=true;
end
if nargin<=5
    cutEmptyTime=false;
end
if nargin<=6
    cutEmptyNeurons=true;
end

to_delete=[];
icNew=zeros(size(ic));
icNew(1:2,:)=ic(1:2,:);
tNew=[];
for i=1:size(ic,2),
    tChannel=t(ic(3,i):ic(4,i));
    tCut=tChannel(find(tChannel>tStart & tChannel<tEnd));
    if isempty(tCut) & cutEmptyNeurons,
        if outMessage
            fprintf('Neuron %d %d does not fire during selected interval and therfore was deleted\n',ic(1:2,i));
        end
        to_delete=[to_delete i];
    else
        icNew(3,i)=length(tNew)+1;
        tNew=[tNew, tCut];
        icNew(4,i)=length(tNew);
    end
end
icNew(:,to_delete)=[];
if cutEmptyTime
    tNew=tNew-tStart;
    if outMessage
        fprintf('Please remember that cut segment started at time %d in the whole sort Channel\n',tStart);
    end
end

    