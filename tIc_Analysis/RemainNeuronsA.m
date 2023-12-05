function [I,t,ic]=RemainNeuronsA(I,t,ic,varargin)
% [I,t,ic]=RemainNeuronsA(I,t,ic,varargin);
% Function purpose : deletes all neurons except the ones mentioned.
%
% Recives :     I - Activity intensity function
%                           t [units of choice] - firing timings
%                           ic - index chanel
%                           varargin - the list of channels to not be deleted: e.g. [12;1],[15;3],....
%                                                    can be entered also as [1XN] vector (assuming every channel has only one neuron) : [2 13 14 16 25...] 
%                                                                                                                                                                   
% Function give back :  newI - new activity intensity function
%                                                   newt [units of choice] - new firing timings
%                                                   newic - new index chanel
% Recomended usage  :   [newI,newt,newic]=RemainNeuronsA(I,t,ic,[2;1],[2;3],[2;4]);
% Last updated : 25/01/10

if length(varargin)>=1
    if size(varargin{1})==2
        for i=1:length(varargin)
            tmp(:,i)=varargin{i};
        end
    else
        for i=1:length(varargin{1})
            tmp(:,i)=[varargin{1}(i);1];
        end
    end
else
    error('No channels in input');
end
x=[];
for i=1:size(ic,2)
    index=find(ic(1,i)==tmp(1,:) & ic(2,i)==tmp(2,:));
    if isempty(index)
        x=[x [ic(1,i);ic(2,i)]];
    end
end

if isempty(varargin),
    error('not enough inputs');
end
%if ~isempty(find(x(1,:)>size(ic,2))),
%    error('channel does not exist');
%end
if size(x,2)>length(ic),
    error('too many inputs');
end

ic_old=ic;
Idx=[];indexes=[];
for i=1:size(x,2),
    %fprintf('removing neuron [%d;%d]\n',x(1,i),x(2,i));
    indexes(i)=find(ic(1,:)==x(1,i) & ic(2,:)==x(2,i));
    %cut out a peace of t:
    Idx=[Idx ic_old(3,indexes(i)):ic_old(4,indexes(i))];
    
    %update indexes:
     Lt=ic(4,indexes(i))-ic(3,indexes(i))+1;
    ic(3,find(ic(3,:)>ic(3,indexes(i))))=ic(3,find(ic(3,:)>ic(3,indexes(i))))-Lt;
    ic(4,find(ic(4,:)>ic(4,indexes(i))))=ic(4,find(ic(4,:)>ic(4,indexes(i))))-Lt;

end

%remove channels:
ic(:,indexes)=[];

t(Idx)=[];
I(Idx)=[];

%%%%%%% To remove neurons that dont fire (instead of the nan assigned to them);
for i=1:size(ic,2)
    if isnan(t(ic(3,i)))
        I(ic(3,i))=[];
        t(ic(3,i))=[];
        ic(:,i)=NaN;
        ic(3:4,(i+1):end)=ic(3:4,(i+1):end)-1;        
    end
end
ic(:,isnan(ic(1,:)))=[];

