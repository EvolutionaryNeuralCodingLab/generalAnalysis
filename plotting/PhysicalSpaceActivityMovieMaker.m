function [fps]=PhysicalSpaceActivityMovieMaker(MovieName,Fs,NodeIntensity,Channels,En,varargin)
% [fps]=PhysicalSpaceActivityMovieMaker(MovieName,NodeIntensity,Channels,En);
% Function purpose : Makes movies of nodes and links on physical space
%
% Recives :     MovieName - the name of the movie
%               Fs - sampling frequency
%               NodeIntensity -  [M X N] with the intensity of every channel point N, in every frame M (if Intensity is empty on channel marker are added).
%               Channels - the channels vector
%               En - Electrode position matrix
%               varargin - list of arguments in ('property name',property value) format
%                   Compression [0/1] - if 1 uses compression.
%                   TimeBin - the real time of every colomn in ActM [sec]
%                   TimeFactor - the time factor used for slowing or increasing the movie speed,
%                       if 0.1 is chosen the movie will be 10 times slower than in reality, if 5 is chosen the movie will
%                       be 5 times faster than in reality
%                   cMap - color map
%                   RImage - the real network image, image should be adjusted to fit the electrode grid             
%                   Allignment [2 X 3] - two allignment markers of electrode coordinates on the overlayed figure (if one is added)
%                   	every marker includes [ElectrodeNumber X-coordinate Y-coordinate]
%                   	if no allignment is needed put empty value []
%                   DrawWhat - see "ConnectivityPhysicalSpacePlot"
%
% Function give back :fps - frames per second in movie
% Recomended usage  :
% Last updated :

%default parameters
cMapName='MyJetColorMap';
DrawWhat=2;
TimeFactor=0.2;
fps=60;
RImage=[]; %the real network image, image should be adjusted to fit the electrode grid
Allignment=[]; %two allignment markers of electrode coordinates on the overlayed figure (if one is added) every marker includes [ElectrodeNumber X-coordinate Y-coordinate]if no allignment is needed put empty value []
mode='elevationPlot';%optional modes: 'elevationPlot','colorCoded'
normalizeActivity=true;
hFig=[];

%print out default arguments and values if no inputs are given
if nargin==0
    defaultArguments=who;
    for i=1:numel(defaultArguments)
        eval(['defaultArgumentValue=' defaultArguments{i} ';']);
        disp([defaultArguments{i} ' = ' num2str(defaultArgumentValue)]);
    end
    return;
end

[nChannels,nSamples]=size(NodeIntensity);

%Collects all options
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end

load(cMapName);

%Normalization
if normalizeActivity
    NodeIntensity=(NodeIntensity-min(min(NodeIntensity)))/(max(max(NodeIntensity))-min(min(NodeIntensity)));
end

%Setting figure properties and calculating first frame
scrsz = get(0,'ScreenSize');
%f=figure('Position',[20 70 scrsz(3)/1.4 scrsz(4)/1.2]);
if isempty(hFig)
    f=figure('Position',[100 100 400 400]);
else
    figure(hFig);
end

h=axes;
if ~isempty(RImage)
    image(RImage);axis equal;
    set(h,'XTick',[]);
    set(h,'YTick',[]);
    whitebg(f,'k');
    set(f,'color','k');
end

axis tight;

%set(gca,'CLimMode','Manual');
set(h,'nextplot','replacechildren');
set(h,'position',get(h,'OuterPosition'));
set(gcf,'Renderer','zbuffer');

%Startfrom the first frame before grabing the first image to correct initial mismatch
IntensityPhysicalSpacePlot(Channels,NodeIntensity(:,1),En,'plotElectrodeNumbers',0,'plotGridLines',0,'plotColorBar',0,'h',h);

writerObj = VideoWriter([MovieName '.avi']);
writerObj.FrameRate = fps;
open(writerObj);

% Record the movie
maxFrame=floor((nSamples/Fs)*fps/TimeFactor);
for i = 1:maxFrame
    ip=ceil((i/fps)*Fs*TimeFactor); %the location in the array
    disp(['Frame ' num2str(i)]);
    IntensityPhysicalSpacePlot(Channels,NodeIntensity(:,ip),En,'plotElectrodeNumbers',0,'plotGridLines',0,'plotColorBar',0,'h',h,'Ilim',[0 1]);

    frame=getframe;
    writeVideo(writerObj,frame)
end
close(writerObj);

fprintf('\nFinished\n');

