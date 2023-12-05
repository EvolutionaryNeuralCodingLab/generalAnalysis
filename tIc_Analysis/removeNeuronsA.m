function [newI,newt,newic]=RemoveNeuronsA(I,t,ic,varargin)
% [newI,newt,newic]=RemoveNeuronsA(I,t,ic,varargin);
% Function purpose : deletes all neurons except the ones mentioned.
%
% Recives :     I - Activity intensity function
%                           t [units of choice] - firing timings
%                           ic - index chanel
%                           varargin - the list of channels to be deleted: e.g. [12;1],[15;3],....
%                                                    can be entered also as [1XN] vector (assuming every channel has only one neuron) : [2 13 14 16 25...] 
%                                                                                                                                                                   
% Function give back :  newI - new activity intensity function
%                                                   newt [units of choice] - new firing timings
%                                                   newic - new index chanel
% Recomended usage  :   [newI,newt,newic]=RemoveNeuronsA(I,t,ic,[2;1],[2;3],[2;4]);
% Last updated : 16/09/09

if length(varargin)>=1
    if size(varargin{1},1)==2
        for i=1:length(varargin)
            x(:,i)=varargin{i};
        end
    else
        for i=1:length(varargin{1})
            x(:,i)=[varargin{1}(i);1];
        end
    end
end
x(3,:)=0;

if isempty(varargin{1}),
    disp('not enough inputs!!!!!!!!!!');
    newI=I;newt=t;,newic=ic;
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
        fprintf('removing neuron [%d;%d]\n',x(1,i),x(2,i));
        %cut out a peace of t:
        cutt=t(ic(3,index):ic(4,index));
        if ic(3,index)~=1,
            tBefore=t(1:ic(3,index)-1);
            IBefore=I(1:ic(3,index)-1);
        else
            tBefore=[];
            IBefore=[];
        end
        if ic(4,index)~=length(t),
            tAfter=t(ic(4,index)+1:end);
            IAfter=I(ic(4,index)+1:end);
        else
            tAfter=[];
            IAfter=[];
        end
        t=[tBefore,tAfter];
        I=[IBefore,IAfter];
        %update indexes:
        ic(3,find(ic(3,:)>ic(3,index)))=ic(3,find(ic(3,:)>ic(3,index)))-length(cutt);
        ic(4,find(ic(4,:)>ic(4,index)))=ic(4,find(ic(4,:)>ic(4,index)))-length(cutt);
    else
        fprintf('Neuron [%d;%d] was not removed since it does not exist\n',x(1,i),x(2,i));
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
newI=I;
newic=ic;

%%%%%%% To remove neurons that dont fire (instead of the nan assigned to them);
for i=1:size(newic,2)
    if isnan(newt(newic(3,i)))
        newI(newic(3,i))=[];
        newt(newic(3,i))=[];
        newic(:,i)=NaN;
        newic(3:4,(i+1):end)=newic(3:4,(i+1):end)-1;        
    end
end
newic(:,isnan(newic(1,:)))=[];        
    

