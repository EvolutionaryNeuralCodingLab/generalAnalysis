% [hCbar]=IntensityPhysicalSpacePlot(Channels,Intensity,En)
% Function purpose : Plots an intensity channel map on real space
%
% Function recives :    Channels - the channels vector
%                       Intensity - intensity vector
%                       En - Electrode setup matrix
%                       varargin - 'property',value
%                           h - axis handle
%                           markerSize
%                           txtHorizontalAlignment
%                           fontSize
%                           plotGridLines
%                           plotElectrodeNumbers
%                           plotColorBar
%                           Ilim
%                           saturations
%
% Function give back :  Intensity plot on real space
%
% Last updated : 11/7/09
function [hCbar,h,hSbar]=IntensityPhysicalSpacePlot(Channels,Intensity,En,varargin)
hCbar=[];
markerSize=50; %if vector, plots sizes in addition to color
txtHorizontalAlignment='center';
fontSize=12;
plotGridLines=1;
plotElectrodeNumbers=1;
plotColorBar=1;

plotSizeBar=1;
logMarkerSize=false;
nSizes2Show=3; %number of markers to show in legened if plotSizeBar is true
maxLogMarkerSizeOnPlot=50; %this increases the maximal size of the marker on the plot and scales all linearly

Ilim=0;
cMap=[];

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

if ~exist('h','var')
    h=axes;
else
    axes(h);
end
if Ilim==0
    Act_norm=(Intensity-min(Intensity))./max(Intensity-min(Intensity));
    Ilim=[min(Intensity) max(Intensity)];
else
    Intensity(Intensity<Ilim(1))=Ilim(1);
    Intensity(Intensity>Ilim(2))=Ilim(2);
    Act_norm=(Intensity-Ilim(1))./(Ilim(2)-Ilim(1));
end

%creating translation between electrodes and locations
translation=[];
En2=fliplr(En);
Elecs=sort(En(:));
Elecs=Elecs(~isnan(Elecs));
for i=1:length(Elecs)
    [n,m]=find(En==Elecs(i));
    translation=[translation;m n Elecs(i)];
end
[max_grid_y,max_grid_x]=size(En);
if plotGridLines
    for i=1:(max_grid_y-1)
        line([0 max_grid_x+1],[i i],'LineWidth',2,'Color',[.8 .8 .8]);
    end
    for i=1:(max_grid_x-1)
        line([i i],[0 max_grid_y+1],'LineWidth',2,'Color',[.8 .8 .8]);
    end
end
xlim([0 max_grid_x]);ylim([0 max_grid_y]);
set(h,'Box','on');
set(h,'xtick',[]);
set(h,'ytick',[]);
set(h,'XColor',[.8 .8 .8]);
set(h,'YColor',[.8 .8 .8]);
axis(h,'square');

%check for colormap
if isempty(cMap)
    cMap=colormap;
end
nColors=size(cMap,1);

if plotColorBar
    IlimTxt=num2str([Ilim(1) round(100*Ilim(2))/100]',3);
    hCbar=colorbar;set(hCbar,'ytick',[0 1]);set(hCbar,'yticklabel',IlimTxt);
    colormap(hCbar,cMap);
    hCbar.Position=[0.8360    0.6550    0.0280    0.2707];
end
hold on;

if exist('saturations','var')
    saturations=saturations/max(saturations);
else
   saturations=ones(1,length(Channels));
end

if numel(markerSize)==1
    markerSize=markerSize*ones(1,length(Channels));
    plotSizeBar=false;
end

logMultiplier=10;
originalSize=markerSize;
if logMarkerSize
    markerSize=logMultiplier*log10(1+markerSize/logMultiplier)+eps;
    markerScaling=max(markerSize)/maxLogMarkerSizeOnPlot;
    markerSize=markerSize/markerScaling;
end

sizMultiplier=10;
if plotSizeBar
    xl=xlim;yl=ylim;
    minSiz=min(originalSize);
    maxSiz=max(originalSize);
    barSizes=[(sizMultiplier*minSiz+maxSiz)/(sizMultiplier+1) (minSiz+sizMultiplier*maxSiz)/(sizMultiplier+1)];
    barSizes=barSizes(1):(barSizes(2)-barSizes(1))/(nSizes2Show-1):barSizes(2);
    markerTxt=num2str(barSizes',2);
    if logMarkerSize
        barSizes=logMultiplier*log10(1+barSizes/logMultiplier)+eps;
        barSizes=barSizes/markerScaling; %point should have the same normalization
    end
    sx=xl(2)+(xl(2)-xl(1))*0.05;
    sxTxt=ones(1,nSizes2Show)*(xl(2)+(xl(2)-xl(1))*0.1);
    sy=(yl(1)+yl(2))/2+(yl(2)-yl(1))*(-0.45:(0.45/nSizes2Show+eps):0);
    for i=1:nSizes2Show,hSbar(i)=plot(sx,sy(i),'.k','markersize',barSizes(i));end;
    text(sxTxt,sy,markerTxt)
    set(h,'clipping','off');
end


% Plotting Activity points
for i=1:length(Channels)
    x=0.5+translation(find(translation(:,3)==Channels(i)),1)-1;
    y=0.5+translation(find(translation(:,3)==Channels(i)),2)-1;
    %check if possible to use scatter instead (make sure that ball sizes are correctly plotted with scatter
    if ~isnan(Act_norm(i))
        colorHSV=rgb2hsv(cMap(ceil(1+Act_norm(i)*(nColors-1)),:));
        colorHSV(2)=colorHSV(2)*saturations(i);
        colorRGB=hsv2rgb(colorHSV);
        plot(x,y,'.','color',colorRGB,'markersize',markerSize(i));
    end
end

%Plotting numbers on top of the propagation points
if plotElectrodeNumbers
    for i=1:length(Elecs)
        if Elecs(i)<10
            text(translation(i,1)-0.5,translation(i,2)-0.5,num2str(Elecs(i)),'fontsize',fontSize,'HorizontalAlignment',txtHorizontalAlignment);
        else
            text(translation(i,1)-0.5,translation(i,2)-0.5,num2str(Elecs(i)),'fontsize',fontSize,'HorizontalAlignment',txtHorizontalAlignment);
        end
    end
end

hold off;