function [tNew,icNew]=RemainNeurons(t,ic,varargin)
% [t,ic]=RemainNeurons(t,ic,varargin);
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
% Recomended usage  :   [newt,newic]=RemainNeurons(t,ic,[2;1],[2;3],[2;4]);
%                                                 [newt,newic]=RemainNeurons(t,ic,[1 4 23]);
% Last updated : 17/12/14

%get neurons to keep
tmp=[];
if numel(varargin)>0 %checks if there are enough inputs
    %if channels are given in a format of [ch;neuron]
    if size(varargin{1},1)==2
        if size(varargin{1},2)==1
            for i=1:length(varargin)
                tmp(:,i)=varargin{i};
            end
        else
            tmp=varargin{1};
        end
    else
        %if channels are given as a vector of channels
        for i=1:length(varargin{1})
            p=find(ic(1,:)==varargin{1}(i));
            v=ic(1:2,p);
            tmp=[tmp v];
        end
    end
else
    error('No channels to keep in input');
end

%calculate output channels
x=[];counter=0;
nNeuorns=size(tmp,2);
tNew=cell(1,nNeuorns);
icNew=cell(1,nNeuorns);
for i=1:nNeuorns
    idx2Keep=find(ic(1,:)==tmp(1,i) & ic(2,:)==tmp(2,i));
    tNew{i}=t(ic(3,idx2Keep):ic(4,idx2Keep));
    icNew{i}=[ic(1:2,idx2Keep);counter+1;counter+numel(tNew{i})];
    counter=counter+numel(tNew{i});
end
tNew=cell2mat(tNew);
icNew=cell2mat(icNew);