function [newT,newIndexchannel]=reduceNeurons(t,indexchannel,varargin);
%function [newT,newIndexchannel]=reduceNeurons(t,indexchannel,varargin);
% a function that reduce size of indexchannel. takes varagin{2:end} and
% adds them to varargin{1}
%example:  [newT,newIndexchannel]=reduceNeurons(t,indexchannel,[2;1],[2;3],[2;4]) -
%takes neurons [2;3] and [2;4] in indexchannel and adds them to index [2;1].

for i=1:length(varargin)
    x(:,i)=varargin{i};
end

if length(x)<2,
    error('not enough inputs');
end
%if ~isempty(find(x(1,:)>length(indexchannel))),
%    error('channel does not exist');
%end
if length(x)>length(indexchannel),
    error('too many inputs');
end

channel=find(indexchannel(1,:)==x(1,1) & indexchannel(2,:)==x(2,1));
for i=2:length(x),
    fprintf('reducing neuron [%d;%d] into neuron [%d;%d]\n',x(1,i),x(2,i),x(1,1),x(2,1));
    index=find(indexchannel(1,:)==x(1,i) & indexchannel(2,:)==x(2,i));
    %cut out a peace of t:
    cutT=t(indexchannel(3,index):indexchannel(4,index));
    if indexchannel(3,index)~=1,
        tBefore=t(1:indexchannel(3,index)-1);
    else
        tBefore=[];
    end
    if indexchannel(4,index)~=length(t),
        tAfter=t(indexchannel(4,index)+1:end);
    else
        tAfter=[];
    end
    t=[tBefore,tAfter];
    %update indexes:
    indexchannel(3,find(indexchannel(3,:)>indexchannel(3,index)))=indexchannel(3,find(indexchannel(3,:)>indexchannel(3,index)))-length(cutT);
    indexchannel(4,find(indexchannel(4,:)>indexchannel(4,index)))=indexchannel(4,find(indexchannel(4,:)>indexchannel(4,index)))-length(cutT);
    %and paste it in a new location:
    tBefore=t(1:indexchannel(4,channel));
    if indexchannel(3,channel)~=length(t),
        tAfter=t(indexchannel(4,channel)+1:end);
    else
        tAfter=[];
    end
    t=[tBefore,cutT,tAfter];
    %update indexes:
    indexchannel(3,find(indexchannel(3,:)>indexchannel(3,channel)))=indexchannel(3,find(indexchannel(3,:)>indexchannel(3,channel)))+length(cutT);
    indexchannel(4,find(indexchannel(4,:)>indexchannel(4,channel)))=indexchannel(4,find(indexchannel(4,:)>indexchannel(4,channel)))+length(cutT);
    indexchannel(4,channel)=indexchannel(4,channel)+length(cutT);
end

%remove channels:
for i=2:length(x),
    indexes(i-1)=find(indexchannel(1,:)==x(1,i) & indexchannel(2,:)==x(2,i));
end
indexchannel(:,indexes)=[];

newT=t;
newIndexchannel=indexchannel;
    

