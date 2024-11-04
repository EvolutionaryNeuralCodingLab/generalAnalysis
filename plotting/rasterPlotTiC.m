%[hPlot]=rasterPlotTiC(h,t,ic,neu,tStart,tEnd,varargin)
% Function purpose : Plots the raw voltage data on the real space of the electrode
%
% Function recives :    h - axes handle
%                       M - raster Matrix
%                       varargin - format 'property','value',...
%
% Function give back :  h - axis handle
%
% Last updated : 22/03/15
function [hPlot,scaleFac]=activityTracePhysicalSpacePlot(h,Channels,Waveform,En,varargin)
%% default variables


drawElectrodeCenters=false;
electrodeColor=[0.9 0.1 0.1];
electrodeSize=5;
%% Output list of default variables
%print out default arguments and values if no inputs are given
if nargin==0
    defaultArguments=who;
    for i=1:numel(defaultArguments)
        eval(['defaultArgumentValue=' defaultArguments{i} ';']);
        if numel(defaultArgumentValue)==1
            disp([defaultArguments{i} ' = ' num2str(defaultArgumentValue)]);
        else
            fprintf([defaultArguments{i} ' = ']);
            disp(defaultArgumentValue);
        end
    end
    return;
end
%% Collects all input variables
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end
%% Main code
hPlot1=[];hPlot2=[];hPlot3=[];hPlot4=[];
if isempty(h)
    h=axes;
else
    axes(h);
end
hold(h,'on');

%Scaling parameter in case the data should be super imposed on a figure. see ConnectivityPhysicalSpacePlot
aX=1;
aY=1;
bX=0;
bY=0;
%En2=fliplr(En);

nCh=numel(Channels);
if nCh~=size(Waveform,1)
    error('The length of channel names does not match the Waveform matrix!!!');
end

%creating translation between electrodes and locations
nSamples=size(Waveform,2);
translation=NaN(nCh,3);
for i=1:nCh
    [n,m]=find(En==Channels(i));
    if ~isempty(n)
        translation(i,:)=[m n Channels(i)];
    else
        translation(i,:)=-Channels(i);
    end
end
[max_grid_y,max_grid_x]=size(En);

XL=[0.5 max_grid_x+0.5];
YL=[0.5 max_grid_y+0.5];

W=0.01*(XL(2)-XL(1)+YL(2)-YL(1))/2;

if DrawGrid
    xM=(1.5:(max_grid_y-0.5))';
    yM=(1.5:(max_grid_x-0.5))';
    
    if ~isempty(xM) & plotYGrid
        hPlot1=line([XL(1) XL(2)],[xM*aY+bY xM*aY+bY]','LineWidth',gridLineWidth,'Color',gridColor,'Parent',h);
    end
    if ~isempty(yM) & plotXGrid
        hPlot2=line([yM *aX+bX yM*aX+bX]',[YL(1) YL(2)],'LineWidth',gridLineWidth,'Color',gridColor,'Parent',h);
    end
end
%xlim([0 max_grid_x]);ylim([0 max_grid_y]);
set(h,'Box','on');
set(h,'xtick',[]);
set(h,'ytick',[]);
set(h,'XColor',[.8 .8 .8]);
set(h,'YColor',[.8 .8 .8]);

% Plotting Activity points
NW=size(Waveform,2);
WaveformShiftX=(1/NW):(1/NW):1;
if strcmp(scaling,'std')
    if isempty(scaleFac)
        scaleFac=[nanmean(Waveform(:)) 10*nanstd(Waveform(:))];
    end
    Waveform=(Waveform-scaleFac(1))./scaleFac(2);
elseif strcmp(scaling,'std0Mean')
    if isempty(scaleFac)
        scaleFac=10*nanstd(Waveform(:));
    end
    Waveform=Waveform./scaleFac(1);
elseif strcmp(scaling,'noOverlap')
    if isempty(scaleFac)
        scaleFac=max(Waveform(:))-min(Waveform(:));
    end
    Waveform=Waveform/scaleFac;
elseif strcmp(scaling,'noOverlapLocal')
    if isempty(scaleFac)
        scaleFac=max(Waveform,[],2)-min(Waveform,[],2);
    end
    Waveform=bsxfun(@rdivide,Waveform,scaleFac);
elseif strcmp(scaling,'manual')
    Waveform=bsxfun(@rdivide,Waveform-(scaleFac(1)+scaleFac(2))/2,scaleFac(2)-scaleFac(1));
    if isempty(scaleFac)
        error('Manual scaling mode requires a scaling factor (scaleFac) input [min max]');
    end
elseif strcmp(scaling,'none')
    scaleFac=1;
end

x=aX*translation(:,1)+bX;
y=aY*translation(:,2)+bY;
hPlot3=plot(ones(nSamples,1)*x'-0.5+WaveformShiftX'*ones(1,nCh),ones(nSamples,1)*y'+Waveform','Parent',h,'color',traceColor,'LineWidth',LineWidth);

if drawElectrodeCenters
    plot(x,y,'.','color',electrodeColor,'MarkerSize',electrodeSize);
end

if ~plotTracesAboveGrid %move grid to be above traces
    uistack(hPlot1,'top');
    uistack(hPlot2,'top');
end
%Plotting numbers on top of the propagation points
if size(Channels,2)>1
    Channels=Channels';
end
if DrawElectrodeNumbers
    hPlot4=text(aX*(translation(:,1)-0.40)+bX,0.35+aY*(translation(:,2)-0.6)+bY,num2str(Channels),'fontsize',8,'Parent',h,'horizontalAlignment','center');
end

if lockXYRatio
    axis equal;
end
xlim(h,XL);
ylim(h,YL);

hPlot=[hPlot1;hPlot2;hPlot3;hPlot4];
hold(h,'off');
