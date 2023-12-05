function [newIndexchannel]=renameNeuron(indexchannel,varargin);
%function [newIndexchannel]=renameNeuron(indexchannel,oldname,newname);
% a function that reduce size of indexchannel. takes varagin{2:end} and
% adds them to varargin{1}
%example:  [newT,newIndexchannel]=reduceNeurons(t,indexchannel,[2;1],[2;3],[2;4]) -
%takes neurons [2;3] and [2;4] in indexchannel and adds them to index [2;1].

for i=1:length(varargin)
    x(:,i)=varargin{i};
end

if length(x)~=2,
    error('improper inputs');
end

fprintf('rename  neuron [%d;%d] to [%d,%d]\n',x(1,1),x(2,1),x(1,2),x(2,2));
channel=find(indexchannel(1,:)==x(1,1) & indexchannel(2,:)==x(2,1));
indexchannel(1:2,channel)=x(:,2);
newIndexchannel=indexchannel;
    

