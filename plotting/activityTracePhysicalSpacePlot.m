%[hPlot]=activityTracePhysicalSpacePlot(h,Channels,Waveform,En,varargin)
% Function purpose : Plots the raw voltage data on the real space of the electrode
%
% Function recives :    h - axes handle
%                       Channels - the channels vector
%                       Waveform [NChannels,Time] - the raw voltage samples of all channels
%                       En - Electrode setup matrix
%                       varargin - format 'property','value',...
%
% Function give back :  h - axis handle
%
% Last updated : 22/03/15
function [hPlot,scaleFac,En]=activityTracePhysicalSpacePlot(h,Channels,Waveform,En,varargin)
%% default variables
DrawGrid=1;
DrawElectrodeNumbers=0;
scaling='std'; %'noOverlap','noOverlapLocal','none','manual'
traceColor=[0.1 0.1 0.9];
gridColor=[0.8 0.8 0.8];
LineWidth=1;
scaleFac=[];
gridLineWidth=2;
plotYGrid=true;
plotXGrid=true;
plotTracesAboveGrid=true;
lockXYRatio=false;
averageSubstruction=true; %
averageSubstructionSamples=[]; % if empty uses all samples for substraction

backgroundImgData=[]; %

drawElectrodeCenters=false;
electrodeColor=[0.9 0.1 0.1];
electrodeSize=5;

%Scaling parameter in case the data should be super imposed on a figure. see ConnectivityPhysicalSpacePlot
aX=1;
aY=1;
bX=0;
bY=0;
rot=0; %rotation in degrees
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

%remove empty lines in En
p1=all(isnan(En));
p2=all(isnan(En'));
En(:,p1)=[];En(p2,:)=[];

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


XL=[0.5 max_grid_x+0.5]*aX+bX;
YL=[0.5 max_grid_y+0.5]*aY+bY;

if rot~=0
    rotMat=[cos(rot/180*pi) -sin(rot/180*pi);sin(rot/180*pi) cos(rot/180*pi)];
    [XL, YL] = rotMat*[XL;YL];
end

if DrawGrid
    xM=(1.5:(max_grid_y-0.5))'*aY+bY;
    yM=(1.5:(max_grid_x-0.5))'*aY+bY;
    
    if rot~=0
        [xM, yM] = tformfwd(tForm, xM, yM);
    end
    
    if ~isempty(xM) & plotYGrid
        hPlot1=line([XL(1) XL(2)],[xM xM]','LineWidth',gridLineWidth,'Color',gridColor,'Parent',h);
    end
    if ~isempty(yM) & plotXGrid
        hPlot2=line([yM yM]',[YL(1) YL(2)],'LineWidth',gridLineWidth,'Color',gridColor,'Parent',h);
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
WaveformShiftX=(((1/NW):(1/NW):1)-0.5)*aX+bX;

if averageSubstruction
    if isempty(averageSubstructionSamples)
        averageSubstructionSamples=1:NW;
    end
    Waveform=bsxfun(@minus,Waveform,mean(Waveform(:,averageSubstructionSamples),2));
end

if strcmp(scaling,'std')
    if isempty(scaleFac)
        scaleFac=[nanmedian(Waveform(:)) 10*nanstd(Waveform(:))];
    end
    Waveform=(Waveform-scaleFac(1))./scaleFac(2);
elseif strcmp(scaling,'std0Mean')
    if isempty(scaleFac)
        scaleFac(2)=10*nanstd(Waveform(:));
        scaleFac(1)=0;
    end
    Waveform=Waveform./scaleFac(2);
elseif strcmp(scaling,'noOverlap')
    if isempty(scaleFac)
        scaleFac(2)=max(Waveform(:))-min(Waveform(:));
        scaleFac(1)=min(Waveform(:))+scaleFac(2)/2;
    end
    Waveform=(Waveform-scaleFac(1))/scaleFac(2);
elseif strcmp(scaling,'noOverlapLocal')
    if isempty(scaleFac)
        scaleFac=max(Waveform,[],2)-min(Waveform,[],2);
    end
    Waveform=bsxfun(@rdivide,Waveform,scaleFac);
elseif strcmp(scaling,'manual')
    if isempty(scaleFac)
        error('Manual scaling mode requires a scaling factor (scaleFac) input [min max]');
    end
    Waveform=bsxfun(@rdivide,Waveform-(scaleFac(1)+scaleFac(2))/2,scaleFac(2)-scaleFac(1));
elseif strcmp(scaling,'none')
    scaleFac=1;
end

Waveform=bY+aY*Waveform';
x=aX*translation(:,1)+bX;
y=aY*translation(:,2)+bY;

if rot~=0
    [x, y] = tformfwd(tForm, x, y);
    Waveform(:)=tformfwd(tForm, Waveform(:));
end

hPlot3=plot(ones(nSamples,1)*x'+WaveformShiftX'*ones(1,nCh),ones(nSamples,1)*y'+Waveform,'Parent',h,'color',traceColor,'LineWidth',LineWidth);

if drawElectrodeCenters
    plot(x,y,'.','color',electrodeColor,'MarkerSize',electrodeSize);
end

if ~isempty(backgroundImgData)
    tForm=maketform('composite',flipud(cell2mat(backgroundImgData.transforms)));
    transImg = imtransform(backgroundImgData.img,tForm,'FillValues',0,'XData',[0 max(x)+100], 'YData',[0 max(y)+100],'XYScale',1);
    hImg=imshow(transImg,'Parent',h);
    uistack(hImg,'bottom');
    h.YDir='normal';
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
    hPlot4=text(aX*(translation(:,1)-0.40)+bX,aY*(0.35+translation(:,2)-0.6)+bY,num2str(Channels),'fontsize',8,'Parent',h,'horizontalAlignment','center');
end

if lockXYRatio
    axis equal;
end
xlim(h,XL);
ylim(h,YL);

hPlot=[hPlot1;hPlot2;hPlot3;hPlot4];
hold(h,'off');
