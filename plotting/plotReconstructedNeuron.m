function [hand,outTree]=plotReconstructedNeuron(anatomicalData,varargin)
%% default variables
IE_identity='E';
reconstructedFrom='Ventricle'; %'Ventricle'/'Pia'
rotation=[0 90];
dZ=10;
addScaleBarToFig=false;
colorReferencePlane=[];
ZcolorLim=[];
plotContour=false;
plotNoregion=false;
plotPlane=false;
colorCodeZ=true;
colorE=[0.9 0.1 0.1];
colorI=[0.1 0.1 0.9];
colorCell=[0.5 0.5 0.5];
colorAxon=[0.5 0.5 0.5];
colorContour=[0 0 0];
colorDendrites=[];
positionFileName=[];
savePositionData=false;
zoomFactorXY=1; %sometimes a zoom out/in is needed for the plots (one is no zoom)
hAxis=[];
resampleTree=false;
resamplingInterval=2; %in um


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
NSKToolBoxMainDir=fileparts(which('identifierOfMainDir4NSKToolBox'));
NSKToolBoxMainDir=regexp(NSKToolBoxMainDir,filesep,'split');
NSKToolBoxMainDir=fullfile(NSKToolBoxMainDir{1:end-1});

addpath(genpath([NSKToolBoxMainDir filesep 'TREES2']));
%addpath(genpath('D:\Mark\matlabFunctions\TREES1.15'));

if ~isempty(positionFileName) %save position data if file name for saving was entered
    savePositionData=true;
    if strcmp(positionFileName(end-3:end),'.mat')
        positionFileName=positionFileName(1:end-4);
    end
end
    
%check if ASC file was provided or mat file
[dirName,name,ext]=fileparts(anatomicalData);
if strcmp(ext,'.ASC')
    [treeInitTmp, coords, contoursInit] = neurolucida_tree(anatomicalData);
    if iscell(treeInitTmp)
        treeInit=treeInitTmp{1};
        for i=2:numel(treeInitTmp)
            treeInit = cat_tree(treeInit,treeInitTmp{i});
        end
    else
        treeInit=treeInitTmp;
    end
    if resampleTree
        disp('Resampling tree...');
        treeInit=resample_tree(treeInit, resamplingInterval);
    end
    if strcmp(reconstructedFrom,'Ventricle')
        treeInit = rot_tree(treeInit,[0 180 0]); % rotate along x-axis
    else
        disp('Assuming neuron was reconstructed from Pia');
    end
    
elseif strcmp(ext,'.mat')
    load(anatomicalData);
elseif isempty(ext)
    error('Extension name should be provided in the file name');
end

[outTree]=sep_tree(treeInit);

if isempty(hAxis)
    hand.hFigure=figure;
    hand.hAxis=axes;
else
    hand.hAxis=hAxis;
    hand.hFigure=hand.hAxis.Parent;
end

if isempty(colorDendrites)
    if IE_identity=='I'
        hand.dendrites=plot_tree(outTree.Dendrite,colorE,[],[],[],'-p -3l');
    elseif IE_identity=='E'
        hand.dendrites=plot_tree(outTree.Dendrite,colorI,[],[],[],'-p -3l');
    else
        error('No IE identity');
    end
else
    hand.dendrites=plot_tree(outTree.Dendrite,colorDendrites,[],[],[],'-p -3l');
end

hand.cell=plot_tree(outTree.CellBody,colorCell,[],[],[],'-p -thick -3l');
if isfield(outTree,'Closed') & plotContour
    hand.contour=plot_tree(outTree.Closed,colorContour,[],[],[],'-p -thick -3l');
end
if isfield(outTree,'noregion') & plotNoregion
    hand.contour=plot_tree(outTree.noregion,colorContour,[],[],[],'-p -3l');
end

if colorCodeZ
    if ~isempty(colorReferencePlane)
        ZDist=abs(colorReferencePlane(1:3)'*[outTree.Axon.X';outTree.Axon.Y';outTree.Axon.Z']+colorReferencePlane(4))./norm(colorReferencePlane(1:3));
        if plotPlane
            xp = 500*[1 -1 -1 1]; % Generate data for x vertices
            yp = 500*[1 1 -1 -1]; % Generate data for y vertices
            zp = -1/colorReferencePlane(3)*(colorReferencePlane(1)*xp + colorReferencePlane(2)*yp + colorReferencePlane(4)); % Solve for z vertices data
            patch(xp, yp, zp,'k','facealpha',0.05,'edgeColor','none');
        end
    else
        ZDist=outTree.Axon.Z;
    end
    if isempty(ZcolorLim)
        ZcolorLim=[min(ZDist) max(ZDist)];
    end
    Zlims=round(ZcolorLim(1):dZ:ZcolorLim(2));
    cMap=flipud(jet(numel(Zlims)));
    
    hand.axon=plot_tree(outTree.Axon,[0.5 0.5 0.5],[],[],[],'-3l');
    for k=1:numel(Zlims)-1
        p=find(ZDist>=Zlims(k) & ZDist<Zlims(k+1));
        set(hand.axon(p),'color',cMap(k,:));
    end
    
    colormap(hand.hFigure,cMap);
    hand.cb=colorbar;
    clim=hand.cb.YLim;
    set(hand.cb,'Position',[ 0.8900    0.77    0.013    0.19],'YTick',clim,'YTickLabel',round(ZcolorLim));
    
    %{
    if  strcmp(reconstructedFrom,'Pia')
        set(hand.cb,'YTickLabel',round(ZcolorLim([2 1])));
    end
    %}
    ylabel(hand.cb,'Z [\mum]');

else
    hand.axon=plot_tree(outTree.Axon,colorAxon,[],[],[],'-3l');hold on;
end

view(rotation(1),rotation(2));

if zoomFactorXY~=1
    xl=xlim;yl=ylim;
    xlim([mean(xl)-diff(xl)/zoomFactorXY/2 mean(xl)+diff(xl)/zoomFactorXY/2]);
    ylim([mean(yl)-diff(yl)/zoomFactorXY/2 mean(yl)+diff(yl)/zoomFactorXY/2]);
end

axis off;

if addScaleBarToFig
    hand.hSB=scalebar;
else
    hand.hSB=[];
end

if savePositionData
    disp('Saving data and printing figures...');
    if isempty(positionFileName)
        positionFileName='positionData';
        disp('Saving position data in current folder');
    end
    print(hand.hFigure,positionFileName,'-djpeg','-r1000');
    print(hand.hFigure,positionFileName,'-dtiff','-r1000');
    
    imgInfo=imfinfo([positionFileName '.tif']);
    [X,Y]=meshgrid(1:imgInfo.Width,1:imgInfo.Height);
    xLim=get(gca,'xLim');
    yLim=get(gca,'yLim');
    pos = plotboxpos(gca);
    
    pixSizeX_um=imgInfo.Width*pos(3)./diff(xLim);
    pixSizeY_um=imgInfo.Height*pos(4)./diff(yLim);
    
    pixPosRealSpaceX=X/pixSizeX_um+xLim(1) -pos(1)*imgInfo.Width/pixSizeX_um;
    pixPosRealSpaceY=flipud(Y/pixSizeY_um+yLim(1) -pos(2)*imgInfo.Height/pixSizeY_um);
    
    save(positionFileName,'pixPosRealSpaceX','pixPosRealSpaceY','treeInit');
end