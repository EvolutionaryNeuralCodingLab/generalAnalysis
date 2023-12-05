%[t_new,ic_new]=addNeurons(t,ic,varargin);
%
% Function purpose : adds neurons to t and ic
%
% Function recives :    t - firing timings
%                       ic - index channel
%                       varargin - [ch1;neu1],t1,[ch2;neu2],t2,...
%
% Function give back : t_new , ic_new - t and ic after merging neurons 
%
% Last updated : 13/07/14

function [t_new,ic_new]=MergeNeurons(t,ic,varargin)
%sort neurons in ascending order
times=varargin(2:2:end);
chNames=cell2mat(varargin(1:2:end));

nNeurons=size(chNames,2);
% get the spike times of the merged channels and remove them from the original time list (t_new)
t_new=t;
ic_new=ic;
for i=1:nNeurons
    nNewSpikes=numel(times{i});
    pCh=find(ic_new(1,:)==chNames(1,i));
    if isempty(pCh)
        pInIc=find(ic_new(1,:)>=chNames(1,i),1,'first');
    else
        neuronsInCh=ic_new(2,pCh);
        [~,pSort]=sort([neuronsInCh chNames(2,i)]);
        pInIc=pCh(1)+find(pSort==numel(neuronsInCh)+1)-1;
    end
    pInt=ic_new(3,pInIc);
    t_new=[t_new(1:(pInt-1)) times{i} t_new(pInt:end)];
    ic_new(3:4,pInIc:end)=ic_new(3:4,pInIc:end)+nNewSpikes;
    ic_new=[ic_new(:,1:pInIc-1),...
        [chNames([1 2],i);pInt;pInt+nNewSpikes-1],...
        ic_new(:,pInIc:end)];
end
