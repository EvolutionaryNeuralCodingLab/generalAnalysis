function [newt,newic]=removeNeurons(t,ic,varargin)
% [newt,newic]=RemoveNeurons(t,ic,varargin);
% Function purpose : deletes all neurons except the ones mentioned.
%
% Recives :     t [units of choice] - firing timings
%                           ic - index chanel
%                           varargin - a) the list of channels to be deleted: e.g. [12;1],[15;3],....
%                                      b) [1XN] vector (assuming every channel has only one neuron) : [2 13 14 16 25...] 
%                                      c) [2XN] vector with channel number and neuron number
%                                                                                                                                                                   
% Function give back :  newt [units of choice] - new firing timings
%                                                   newic - new index chanel
% Recomended usage  :   [newt,newic]=RemoveNeurons(t,ic,[2;1],[2;3],[2;4]);
% Last updated : 16/09/09

outputMessages=false;

%get neurons to keep
x=[];
if numel(varargin)>0 %checks if there are enough inputs
    %if channels are given in a format of [ch;neuron]
    if size(varargin{1},1)==2
        if size(varargin{1},2)==1
            for i=1:length(varargin)
                x(:,i)=varargin{i};
            end
        else
            x=varargin{1};
        end
    else
        %if channels are given as a vector of channels
        for i=1:length(varargin{1})
            p=find(ic(1,:)==varargin{1}(i));
            v=ic(1:2,p);
            x=[x v];
        end
    end
else
    error('No channels to keep in input');
end
x(3,:)=0;


if isempty(varargin{1}),
    disp('not enough inputs!!!!!!!!!!');
    newt=t;newic=ic;
    return;
end
%if ~isempty(find(x(1,:)>size(ic,2))),
%    error('channel does not exist');
%end
if size(x,2)>length(ic),
    error('too many inputs');
end

for i=1:size(x,2),
    index=find(ic(1,:)==x(1,i) & ic(2,:)==x(2,i));
    if ~isempty(index)
        if outputMessages
            fprintf('removing neuron [%d;%d]\n',x(1,i),x(2,i));
        end
        %cut out a peace of t:
        cutt=t(ic(3,index):ic(4,index));
        if ic(3,index)~=1,
            tBefore=t(1:ic(3,index)-1);
        else
            tBefore=[];
        end
        if ic(4,index)~=length(t),
            tAfter=t(ic(4,index)+1:end);
        else
            tAfter=[];
        end
        t=[tBefore,tAfter];
        %update indexes:
        ic(3,find(ic(3,:)>ic(3,index)))=ic(3,find(ic(3,:)>ic(3,index)))-length(cutt);
        ic(4,find(ic(4,:)>ic(4,index)))=ic(4,find(ic(4,:)>ic(4,index)))-length(cutt);
    else
        if outputMessages
            fprintf('Neuron [%d;%d] was not removed since it does not exist\n',x(1,i),x(2,i));
        end
        x(3,i)=1;
    end
end
%remove channels:
x(:,x(3,:)==1)=[];
for i=1:size(x,2),
    indexes(i)=find(ic(1,:)==x(1,i) & ic(2,:)==x(2,i));
end
ic(:,indexes)=[];

newt=t;
newic=ic;

%%%%%%% To remove neurons that dont fire (instead of the nan assigned to them);
for i=1:size(newic,2)
    if isnan(newt(newic(3,i)))
        newt(newic(3,i))=[];
        newic(:,i)=NaN;
        newic(3:4,(i+1):end)=newic(3:4,(i+1):end)-1;        
    end
end
newic(:,isnan(newic(1,:)))=[];        