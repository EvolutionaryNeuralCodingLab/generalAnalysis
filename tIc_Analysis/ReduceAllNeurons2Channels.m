%[ic]=ReduceAllNeurons2Channels(ic);
%
% Function purpose : Reduces sort channel to one neurons per channel, by
%                    merging all the neurons at each channel
%
% Function recives :    t - firing timings
%                       ic - index channel
%
% Function give back : t , ic - after reducing to a sigle neurons at each channel
%
% Last updated : 27/12/16

function [t_new,ic_new]=ReduceAllNeurons2Channels(t,ic)
%function assumes that ic is sorted
channels=unique(ic(1,:));
t_new=zeros(size(t));
for i=1:numel(channels)
    pNeurons=find(ic(1,:)==channels(i));
    startIdx=ic(3,pNeurons(1));
    endIdx=ic(4,pNeurons(end));
    ic_new(:,i)=[channels(i);1;startIdx;endIdx];
    t_new(startIdx:endIdx)=sort(t(startIdx:endIdx));
end
