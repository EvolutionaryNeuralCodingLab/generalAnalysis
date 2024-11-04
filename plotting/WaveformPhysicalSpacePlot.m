%[h]=WaveformPhysicalSpacePlot(Channels,Waveform,En,DrawWhat,h)
% Function purpose : Plots the raw voltage data on the real space of the electrode
%
% Function recives :    Channels - the channels vector
%                                                 Waveform [NChannels,Time] - the raw voltage samples of all channels
%                                                 En - Electrode setup matrix
%                                                 DrawWhat - 0-Grid and electrode numbers
%                                                                          1-Grid only
%                                                                          2-Electrode numbers only
%                                                                          3-nothing
%                                                 h - axis handle
%
% Function give back :  h - axis handle
%
% Last updated : 05/01/10
function [h]=WaveformPhysicalSpacePlot(Channels,Waveform,En,DrawWhat,h)

if isempty(DrawWhat) || DrawWhat==0
    DrawElectrodeNumbers=1;
    DrawGrid=1;
elseif DrawWhat==1
    DrawElectrodeNumbers=0;
    DrawGrid=1;
elseif DrawWhat==2
    DrawElectrodeNumbers=1;
    DrawGrid=0;
elseif DrawWhat==3
    DrawElectrodeNumbers=0;
    DrawGrid=0;
else
    error('DrawWhat parameter has a invalid value');
end

%Scaling parameter in case the data should be super imposed on a figure. see ConnectivityPhysicalSpacePlot
aX=1;
aY=1;
bX=0;
bY=0;

if ~exist('h','var')
    h=axes;
else
    axes(h);
end

%creating translation between electrodes and locations
translation=[];
%En2=fliplr(En);
Elecs=sort(En(~isnan(En)));
for i=1:length(Elecs)
    [n,m]=find(En==Elecs(i));
    translation=[translation;m n Elecs(i)];
end
[max_grid_y,max_grid_x]=size(En);


xlim([0.5 max_grid_x+0.5]);
ylim([0.5 max_grid_y+0.5]);

XL=xlim;
YL=ylim;
W=0.01*(XL(2)-XL(1)+YL(2)-YL(1))/2;

if DrawGrid
    for i=1.5:(max_grid_y-0.5)
        TransparentLine([XL(1) XL(2)],[i*aY+bY i*aY+bY],[.9 .5 .2],W,0.4);
        %line([XL(1) XL(2)],[i*aY+bY i*aY+bY],'LineWidth',2,'Color',[.1 .5 .7]);
    end
    for i=1.5:(max_grid_x-0.5)
        TransparentLine([i *aX+bX i *aX+bX],[YL(1) YL(2)],[.9 .5 .2],W,0.4);
        %line([i *aX+bX i *aX+bX],[YL(1) YL(2)],'LineWidth',2,'Color',[.1 .5 .7]);
    end
end
%xlim([0 max_grid_x]);ylim([0 max_grid_y]);
set(gca,'Box','on');
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'XColor',[.8 .8 .8]);
set(gca,'YColor',[.8 .8 .8]);
hold on;

% Plotting Activity points
NW=size(Waveform,2);
WaveformShiftX=(1/NW):(1/NW):1;
YScaling=1/(max(max(Waveform))-min(min(Waveform)));
if ~isempty(Waveform)
    for i=1:length(Channels)
        x(i)=aX*translation(find(translation(:,3)==Channels(i)),1)+bX;
        y(i)=aY*translation(find(translation(:,3)==Channels(i)),2)+bY;
        plot(x(i)-0.5+WaveformShiftX,y(i)+YScaling*Waveform(i,:))
    end
end

%Plotting numbers on top of the propagation points
if DrawElectrodeNumbers
    for i=1:length(Elecs)
        if Elecs(i)<10
            text(aX*(translation(i,1)-0.20)+bX,0.35+aY*(translation(i,2))+bY,num2str(Elecs(i)),'fontsize',8,'FontWeight','Bold');
        else
            text(aX*(translation(i,1)-0.44)+bX,0.35+aY*(translation(i,2))+bY,num2str(Elecs(i)),'fontsize',8,'FontWeight','Bold');
        end
    end
end

