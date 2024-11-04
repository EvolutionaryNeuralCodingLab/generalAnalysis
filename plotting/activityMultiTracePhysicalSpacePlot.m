%[hPlot]=activityMultiTracePhysicalSpacePlot(h,Channels,Waveform,En,varargin)
% Function purpose : Plots the raw voltage data on the real space of the electrode
%
% Function recives :    Channels - the channels vector
%                                                 Waveform [NChannels,Trials,Time] - the raw voltage samples of all channels
%                                                 En - Electrode setup matrix
%                                                 h - axis handle
%
% Function give back :  h - axis handle
%
% Last updated : 23/04/13
function [hPlot]=activityMultiTracePhysicalSpacePlot(h,Channels,Waveform,En,varargin)
DrawGrid=1;
DrawElectrodeNumbers=1;
scaling='stdMedian';
lineWidth=1;
scaleFac=0.8;
%Scaling parameter in case the data should be super imposed on a figure. see ConnectivityPhysicalSpacePlot
aX=1;
aY=1;
bX=0;
bY=0;
hPlot1=[];
hPlot2=[];
hPlot3=[];
hPlot4=[];
hPlot5=[];
yScale=[];

%print out default arguments and values if no inputs are given
if nargin==0
    defaultArguments=who;
    for i=1:numel(defaultArguments)
        eval(['defaultArgumentValue=' defaultArguments{i} ';']);
        disp([defaultArguments{i} ' = ' num2str(defaultArgumentValue)]);
    end
    return;
end

%Collects all options
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end

if ~ishandle(h) %check that a valid axes handle was entered
    error('no valid axes handle provided!');
end

if ~isempty(yScale)
    scaling='manual';
end

if isempty(h)
    h=axes;
end
%creating translation between electrodes and locations
[nCh,nTrials,nSamples]=size(Waveform);

%cmap=colormap(jet(nTrials));
cmap=lines(nTrials);

translation=NaN(nCh,3);
for i=1:nCh
    [n,m]=find(En==Channels(i));
    if ~isempty(n)
        translation(i,:)=[m n Channels(i)];
    else
        translation(i,:)=Channels(i);
    end
end
[max_grid_y,max_grid_x]=size(En);

XL=[0.5 max_grid_x+0.5];
YL=[0.5 max_grid_y+0.5];

xlim(h,XL);
ylim(h,YL);

W=0.01*(XL(2)-XL(1)+YL(2)-YL(1))/2;

if DrawGrid
    xM=(1.5:(max_grid_y-0.5))';
    yM=(1.5:(max_grid_x-0.5))';
    
    hPlot1=line([XL(1) XL(2)],[xM*aY+bY xM*aY+bY]','LineWidth',2,'Color',[0.9 0.8 0.7],'Parent',h);
    hPlot2=line([yM *aX+bX yM*aX+bX]',[YL(1) YL(2)],'LineWidth',2,'Color',[0.9 0.8 0.7],'Parent',h);
end
%xlim([0 max_grid_x]);ylim([0 max_grid_y]);
set(h,'Box','on');
set(h,'xtick',[]);
set(h,'ytick',[]);
set(h,'XColor',[.8 .8 .8]);
set(h,'YColor',[.8 .8 .8]);

% Plotting Activity points
WaveformShiftX=(1/nSamples):(1/nSamples):1;
if strcmp(scaling,'std')
    waveformStd=10*std(Waveform(:));
    waveformMean=mean(Waveform(:));
    Waveform=(Waveform-waveformMean)/waveformStd/nTrials/scaleFac;
    Waveform=bsxfun(@plus,Waveform,(1:nTrials)/(nTrials+1));
elseif strcmp(scaling,'stdMedian')
    waveformStd=10*median(abs(  Waveform(:)-median(Waveform(:))  )) / 0.6745;
    waveformMean=median(Waveform(:));
    Waveform=(Waveform-waveformMean)/waveformStd/nTrials/scaleFac;
    Waveform=bsxfun(@plus,Waveform,(1:nTrials)/(nTrials+1));        
elseif strcmp(scaling,'noOverlap')
    Waveform=0.5+Waveform/(max(Waveform(:))-min(Waveform(:)));
elseif strcmp(scaling,'noOverlapLocal')
    %Waveform=bsxfun(@rdivide,Waveform,(max(Waveform,[],2)-min(Waveform,[],2)));
elseif strcmp(scaling,'manual')
    Waveform=0.5+Waveform/yScale;
end

x=aX*translation(:,1)+bX;
y=aY*translation(:,2)+bY;
hold on;
hPlot3=[];
for i=1:nTrials
    if nCh==1
        hPlotTmp=plot(ones(nSamples,1)*x'-0.5+WaveformShiftX'*ones(1,nCh),ones(nSamples,1)*y'+squeeze(Waveform(:,i,:))-0.5,'Parent',h,'color',cmap(i,:),'lineWidth',lineWidth);
    else
        hPlotTmp=plot(ones(nSamples,1)*x'-0.5+WaveformShiftX'*ones(1,nCh),ones(nSamples,1)*y'+squeeze(Waveform(:,i,:))'-0.5,'Parent',h,'color',cmap(i,:),'lineWidth',lineWidth);
    end
    hPlot3=[hPlot3;hPlotTmp];
end
%Plotting numbers on top of the propagation points
%Plotting numbers on top of the propagation points
if size(Channels,2)>1
    Channels=Channels';
end
if DrawElectrodeNumbers
    hPlot4=text(aX*(translation(:,1)-0.40)+bX,0.35+aY*(translation(:,2)-0.6)+bY,num2str(Channels),'fontsize',8,'Parent',h,'horizontalAlignment','left');
end
hPlot=[hPlot1;hPlot2;hPlot3;hPlot4];
