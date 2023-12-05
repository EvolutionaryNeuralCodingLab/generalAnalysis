function [t,indexchannel]=MergeSortChannel(t1,ic1,t2,ic2);
% [t,indexchannel]=MergeSortChannel(t1,ic1,t2,ic2);
%
% Function purpose : takes two sort chanels and merges them into one each
% channel is taken from both recordings sorted and combined to the rest.
%
% Function recives :    t1,2  - firing timings of input timings
%                       ic1,2 - index channel of t1,2
%                       
% Function give back :  t - firing timings
%                       indexchannel - index channel
% Last updated : 07/06/05

channels=sort(unique([ic1(1,:) ic2(1,:)]));

%Making a vector all channels with all channels and neurons in merged recoding
all_channels=[];
for i=1:length(channels)
    max_neurons=max([length((find(ic1(1,:)==channels(i)))) length((find(ic2(1,:)==channels(i))))]);
    all_channels=[all_channels [(ones(1,max_neurons).*channels(i));(1:max_neurons)]];
end

%going over on all channels combines both recordings and build indexchannel
t=[];
indexchannel=[];
for i=1:size(all_channels,2)
    tmp_ind1=find(ic1(1,:)==all_channels(1,i) & ic1(2,:)==all_channels(2,i));
    t_tmp1=t1(ic1(3,tmp_ind1):ic1(4,tmp_ind1));
    tmp_ind2=find(ic2(1,:)==all_channels(1,i) & ic2(2,:)==all_channels(2,i));
    t_tmp2=t2(ic2(3,tmp_ind2):ic2(4,tmp_ind2));
    if max(t_tmp1)>min(t_tmp2)
        fprintf('\nBoth recordings are one same times or the order is swiched\nRemember the first recording chronologically comes first');
    end
    indexchannel=[indexchannel [all_channels(1:2,i);(length(t)+1);(length(t)+length([t_tmp1 t_tmp2]))]];
    t=[t sort([t_tmp1 t_tmp2])];
end