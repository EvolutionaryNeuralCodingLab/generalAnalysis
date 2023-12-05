%[t,ic]=ReduceSortChannel(t,ic,neurons_not_to_reduce);
%
% Function purpose : Reduces sort channel to one neurons per channel, by
%                    merging all the neurons at each channel
%
% Function recives :    t - firing timings
%                       ic - index channel
%                       neurons_not_to_reduce - ignores reducing these neurons
%                       form of vector : [26 56...] . Defult empty
%
% Function give back : t , ic - after reducing to a sigle neurons at each channel
%
% Last updated : 27/12/05

function [t,ic]=ReduceSortChannel(t,ic,neurons_not_to_reduce);

if nargin<3
    neurons_not_to_reduce=[];
end
channels=unique(ic(1,:));

if ~isempty(neurons_not_to_reduce)
    for i=1:length(neurons_not_to_reduce)
        channels(find(channels==neurons_not_to_reduce(i)))=[];
    end
end
        
for i=1:size(channels,2)
    neurons=ic(2,find(ic(1,:)==channels(i)));
    if size(neurons,2)>1
        for j=neurons(2:end)
            [t,ic]=MergeNeurons(t,ic,[channels(i);neurons(1)],[channels(i);j]);
        end
    end
end
