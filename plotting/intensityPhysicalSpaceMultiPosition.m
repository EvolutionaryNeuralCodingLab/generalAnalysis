% [varOut,hOut]=intensityPhysicalSpaceMultiPosition(A,neuPosXY,catPosXY,varargin)
% Function purpose : Plots the raw voltage data on the real space of the electrode
%
% Function recives :    A - intensity matrix for different categories
%                       neuPosXY - neuron positions
%                       catPosXY - category position (can be in other scales than neuron positions)
%
% Function give back :  h - axis handle
%
% Last updated : 22/03/15
function [varOut,hOut]=intensityPhysicalSpaceMultiPosition(A,neuPosXY,catPosXY,varargin)
%A is [neurons x positions]

%default parameters
parseObj = inputParser;
addParameter(parseObj,'normalizeCategories',false,@isnumeric); %[ms]
addParameter(parseObj,'a',[1;1],@isnumeric);
addParameter(parseObj,'b',[0;0],@isnumeric);
addParameter(parseObj,'rot',0,@isnumeric); %rotation in degrees
addParameter(parseObj,'inputParams',0,@isnumeric); %if one plots input params
addParameter(parseObj,'additionPositionMarker',[],@isnumeric); %[2 x number of categories] - addes an additional spatial marker per category

%plotting parameters
addParameter(parseObj,'h',[],@ishandle); %rotation in degrees
addParameter(parseObj,'gap',0.1,@isnumeric); %rotation in degrees
addParameter(parseObj,'dotSize',[],@isnumeric); %in units of position, if empty calculates automatically
addParameter(parseObj,'plotNeuronNumbers',false,@isnumeric); %plot numbers of neurons/channels on points
addParameter(parseObj,'plotCatNumbers',false,@isnumeric); %plot numbers of categories on subplots
addParameter(parseObj,'colorMap','MyHotColorMap'); %color map either string or M x 3 matrix
addParameter(parseObj,'addGrid',false,@isnumeric); %a rectangle surrounding every catergory
addParameter(parseObj,'filledCircles',true,@isnumeric); %a rectangle surrounding every catergory
addParameter(parseObj,'transparency',0.6,@isnumeric); % alpha transparency of points
addParameter(parseObj,'removeZeroPoints',0,@isnumeric); %a rectangle surrounding every catergory


%Parse variables
parseObj.parse(varargin{:});

%evaluate all input parameters in workspace
if parseObj.Results.inputParams || nargin==0, disp(parseObj.Results), return, end

%make parameter structure
par=parseObj.Results;

%% Main code
if isempty(par.h)
    hOut.axes=axes;
else
    axes(par.h);
    hOut.axes=par.h;
end
hold(hOut.axes,'on');
axis(hOut.axes,'equal','off');

[nNeu,nPos]=size(A);

if nNeu~=size(neuPosXY,2)
    disp('The number of neurons/channels does not correspond to the number of raws in the activity input matrix');
    return;
end
if nPos~=size(catPosXY,2)
    disp('The number of conditions does not correspond to the number of columns in the activity input matrix');
    return;
end

%{
if rot~=0
    rotMat=[cos(par.rot/180*pi) -sin(par.rot/180*pi);sin(par.rot/180*pi) cos(par.rot/180*pi)];
    [XL, YL] = (rotMat*[XL;YL])*par.a+par.b;
end
%}

dX=max(neuPosXY,[],2)-min(neuPosXY,[],2); %width of neuron positions

%determine the minimal distance between any two categories
minDist=dist(catPosXY);
minDist(minDist==0)=Inf;
minDist=min(minDist(:));

minPositionInCategory=[min(neuPosXY(1,:)) min(neuPosXY(2,:))];% min positions for neurons

X=neuPosXY(1,:)'*ones(1,nPos)-minPositionInCategory(1); %X positions for neurons starting from zero
catX=catPosXY(1,:)/minDist*dX(1)*(1+par.gap);%X move category positions to neurons units
X=bsxfun(@plus,X,catX);%X add category position to neurons

Y=neuPosXY(2,:)'*ones(1,nPos)-minPositionInCategory(2);%Y positions for neurons starting from zero
catY=catPosXY(2,:)/minDist*dX(2)*(1+par.gap);%Y move category positions to neurons units
Y=bsxfun(@plus,Y,catY);%Y add category position to neurons

if isempty(par.dotSize)
    allDist=dist(neuPosXY);
    par.dotSize=sqrt(mean(allDist(:)));
end

% Plotting Activity points

if par.normalizeCategories
    A=normZeroOne(A);
end

%plot grid lines
if par.addGrid
    olX=[catX;catX+dX(1);catX+dX(1);catX;catX];
    olY=[catY;catY;catY+dX(2);catY+dX(2);catY];
    hOut.outline=line(olX,olY,'color',[0.7 0.7 0.7]);
end

if par.filledCircles
    hOut.scatter=scatter(X(:),Y(:),par.dotSize,A(:),'filled','MarkerFaceAlpha',par.transparency);
else
    hOut.scatter=scatter(X(:),Y(:),par.dotSize,A(:),'linewidth',1.5,'MarkerEdgeAlpha',par.transparency);
end

if isstr(par.colorMap)
    if strcmp(par.colorMap,'MyHotColorMap')
        cMap=load('MyHotColorMap');cMap=cMap.MyHotColorMap;
    else
        cMap=colormap(par.colorMap);
    end
elseif isnumeric(par.colorMap)
    cMap=par.colorMap;
end
colormap(hOut.axes,cMap);

if par.plotNeuronNumbers
    strCell=mat2cell((1:nNeu)'*ones(1,nPos),ones(1,nNeu),ones(1,nPos))
    hOut.txt=text(X(:),Y(:),strCell,'fontsize',8,'Parent',hOut.axes,'horizontalAlignment','center');
end
if par.plotCatNumbers
    hOut.txtCat=text(catX(:),catY(:),num2str((1:numel(catY))'),'fontsize',8,'Parent',hOut.axes,'horizontalAlignment','center');
end
if ~isempty(par.additionPositionMarker)
    amX=catX+par.additionPositionMarker(1,:)-minPositionInCategory(1);
    amY=catY+par.additionPositionMarker(2,:)-minPositionInCategory(2);
    hOut.additionalMarker=plot(amX,amY,'+k');
end
set(hOut.axes,'Box','off');

hOut.scalebar=addScaleBar(hOut.axes,'scaleBarAxes','x','unitsInScaleBar',0,'XUnitStr','\mum','transparentScale',0);
hOut.cBar=colorbar;
hOut.cBar.Position=[0.8356    0.2048    0.0126    0.1169];
annotation('textbox',hOut.cBar.Position,'FitBoxToText','off','FaceAlpha',par.transparency,'EdgeColor','none','BackgroundColor',[1 1 1]);
ylabel(hOut.cBar,'intensity');

hold(hOut.axes,'off');

if nargout>0
    varOut.dX=dX;
    varOut.X=X;
    varOut.Y=Y;
    varOut.catX=catX;
    varOut.catY=catY;
    varOut.minPositionInCategory=minPositionInCategory;
end

