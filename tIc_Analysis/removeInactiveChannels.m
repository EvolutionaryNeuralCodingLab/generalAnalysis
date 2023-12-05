% [channelsToRemove,channelEventFraction]=removeInactiveChannels(t,ic,res,minimalEventFraction,varargin);
% Function purpose : Finds inactive neurons in a recording
%
% Recives :t [ms] - firing timings
%                       ic - index chanel
%                       res [ms] - resolution of finding bursts (100);
%                       minimalEventFraction - the minimal fraction of time the channel was active. Below this fraction the channel is removed
%                       varargin - 'Plot' - if entered the results are plotted
%                                                                                                                                                        
% Function give back :  channelsToRemove - the channel with activity below the predefined threshold
%                                                   channelEventFraction - the fraction of time every channel was active
% Recomended usage  : [channelsToRemove]=removeInactiveChannels(t,ic,res,5e-5,'Plot');
% Last updated : 20/09/10
function [channelsToRemove,channelEventFraction]=removeInactiveChannels(t,ic,res,minimalEventFraction,varargin)
for i=1:nargin-4
    eval([varargin{i} '=1;']);
end
NChannels=size(ic,2);
Duration=max(t)-min(t);
minimalEventNumber=(minimalEventFraction*res)/Duration;
for i=1:NChannels
    channelEventFraction(i)=(ic(4,i)-ic(3,i))*res/Duration;
end
if exist('Plot','var')
    figure;hold on;
    bar(channelEventFraction)
    set(gca,'XTick',1:NChannels);
    set(gca,'XTickLabel',ic(1,:),'FontSize',8);
    xl=xlim;
    line([xl],[minimalEventFraction minimalEventFraction],'color','r');
    hold off;
    ylim([0 minimalEventFraction*3]);
end
channelsToRemove=ic(1:2,channelEventFraction<minimalEventFraction);
   