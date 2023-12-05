function [newIndexchannel]=flipnameNeurons(indexchannel,varargin);
%function [newIndexchannel]=flipnameNeurons(indexchannel,varargin);
% a function that flip names in indexchannel. takes varagin{2:end} and
% flips with varargin{1}
%example:  [newIndexchannel]=flipnameNeurons(indexchannel,[2;1],[2;3]) -

for i=1:length(varargin)
    x(:,i)=varargin{i};
end

if length(x)~=2,
    error('improper inputs');
end

fprintf('flip names of  neurons [%d;%d],[%d,%d]\n',x(1,1),x(2,1),x(1,2),x(2,2));
channel1=find(indexchannel(1,:)==x(1,1) & indexchannel(2,:)==x(2,1));
channel2=find(indexchannel(1,:)==x(1,2) & indexchannel(2,:)==x(2,2));
indexchannel(1:2,channel1)=x(:,2);
indexchannel(1:2,channel2)=x(:,1);

newIndexchannel=indexchannel;
    

