%[t_new,ic_new]=MergeNeurons(t,ic,varargin);
%
% Function purpose : merges several neurons in sort_channel
%
% Function recives :    t - firing timings
%                       ic - index channel
%
% Function give back : t_new , ic_new - t and ic after merging neurons 
%
% Last updated : 13/07/14

function [t_new,ic_new]=MergeNeurons(t,ic,varargin)
%sort neurons in ascending order
ch2Merge=cell2mat(varargin);

% get the spike times of the merged channels and remove them from the original time list (t_new)
t_new=t;
ic_new=ic;

nNeurons=size(ch2Merge,2);
neuron=zeros(1,nNeurons);
t_neuron=cell(1,nNeurons);
for i=1:nNeurons
    pNeuron=find(ic_new(1,:)==ch2Merge(1,i) & ic_new(2,:)==ch2Merge(2,i));
    if isempty(pNeuron)
        error(['Neuron  ' num2str(ch2Merge([1 2],i)') '  does not appear in ic!']);
    else
        neuron(i)=pNeuron;
        t_neuron{i}=t_new(ic_new(3,neuron(i)):ic_new(4,neuron(i)));
        [t_new,ic_new]=removeNeurons(t_new,ic_new,ch2Merge(:,i));
    end
end

t_neuron=sort(cell2mat(t_neuron));
thresh=1; %Threshold for minimum distance between spikes (in units of t)
double_spikes=find(diff(t_neuron)<thresh);
% t_merged(double_spikes+1)=[];
fprintf('The number of correlated spikes is %d out of %d.\n',numel(double_spikes),numel(t_neuron));

[t_new,ic_new]=addNeurons(t_new,ic_new,ch2Merge(:,1),t_neuron);