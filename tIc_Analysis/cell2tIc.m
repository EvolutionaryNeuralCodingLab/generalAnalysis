function [t,ic]=cell2tIc(c,neuronNames)
% [t,ic]=cellArray2tIc(c,neuronNames)
% Convert cell array with neuron names to t+ic format
if numel(c)~=size(neuronNames,2)
    error('Size of input vectors does not match');
end
nSpk=cellfun(@(x) numel(x),c);
cumSpk=cumsum(nSpk);
ic=[neuronNames;[1 cumSpk(1:end-1)+1];cumSpk];
if size(c{1},2) ~= size(c,2)
    c=c';
end
t=cell2mat(c);
t=t(:)';