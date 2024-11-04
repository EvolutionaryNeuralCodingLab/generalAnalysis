%  [hCbar,h,hNodes,hEdges]=ConnectivityPhysicalSpacePlot(Channels,Intensity ,Connectivity,ConnectivityIntensity,En,DrawWhat,Allignment,h,varargin)
% Function purpose : Plots an connectivity channel map on real space with connecting lines according to intensity
%
% Function recives :    Channels - the channels vector
%                       Intensity [1 X N] - the intensity of every channel point (if Intensity is empty on channel marker are added).
%                       Connectivity [2 x N] or [NxN symettric matrix] 
%                                                       if [2 x N] - a vector for N connecting lines between the ith channel and the jth channel
%                                                       vector should be given as a serial order in Channel (e.g. [1 4 5 3 1 3;2 3 4 2 5 5]);
%                                                       if [NxN symettric matrix] - the connectivity matrix should be symettrix. Edges that should not
%                                                       be displayed have to contain NaNs.
%                       ConnectivityIntensity [1 X N] - a vector of the intensity of the N connecting lines 
%                       En - Electrode setup matrix
%                       DrawWhat - 0-Grid and electrode numbers
%                                                1-Grid only
%                                                2-Electrode numbers only
%                                                3-nothing
%                       Allignment [2 X 3] - two allignment markers of electrode coordinates on the overlayed figure (if one is added)
%                                                      every marker includes [ElectrodeNumber X-coordinate Y-coordinate]
%                                                      if no allignment is needed put empty value []
%                       h - axis handle
%                       varargin - additional arguments: 'Transparancy' - [0 1] - line/dot transparacy
%                                                                                                      'CircleSizeFactor' - Factor for circle size [0 - Inf]
%                                                                                                      'LineSizeFactor' - Factor for line size [0 - Inf]
%                                                                                                      'cm' - color map - 64x3 matrix - [0 1]
%                                                                                                      'ElecNumberSize' - point size for electrode numbers [8 Inf]
%                                                                                                      'NormIntensity' - 'yes'/'no'/'posNeg' - If no, intensity should be [0 1], if 'posNeg'  should be [-1 1]
%                                                                                                      'NoiseFac'-the noise line correlation line places [0 1]
%
% Function give back :  h - axis handle
%                       hcbar - color bar handle
%                       hNodes - the node handles
%                       hEdges - the edge (links) handles
%
% Last updated : 02/10/11
function [hCbar,h,hNodes,hEdges]=ConnectivityPhysicalSpacePlot(Channels,Intensity ,Connectivity,ConnectivityIntensity,En,DrawWhat,Allignment,h,varargin)
%default values
h=[];hNodes=[];hEdges=[];hcbar=[];
NoiseFac=0; %NoiseFactor in the location of the lines (~0.05) 
Transparancy=0.3;
CircleSizeFactor=3;
LineSizeFactor=1;
ElecNumberSize=12;
cm=colormap('jet');
NormIntensity='yes';
plotColorBar=1;
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end

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

if ~exist('h','var') || isempty(h)
    h=axes;
else
    if ishandle(h) %checks that the axes object of the handle h exists (was not deleted)
        axes(h);
    else
       disp('The input axis handle is not valid');
       return;
    end
end

% In case connectivity input was given as a matrix
if size(Connectivity,1)==size(Connectivity,2) && isempty(ConnectivityIntensity)
    [Connectivity,ConnectivityIntensity]=connectivityMatrix2ConnectivityIntensity(Connectivity);
end
if isempty(ConnectivityIntensity)
    ConnectivityIntensity=ones(1,size(Connectivity,2));
end
if isempty(Intensity)
    Intensity=ones(1,length(Channels));
end

XL=xlim;
YL=ylim;
[max_grid_y,max_grid_x]=size(En);

if isempty(Allignment)
    aX=(XL(2)-XL(1))/max_grid_x;
    aY=(YL(2)-YL(1))/max_grid_y;
    bX=-0.5*aX;
    bY=-0.5*aY;
else
    [n1,m1]=find(En==Allignment(1,1));
    [n2,m2]=find(En==Allignment(2,1));
    %According to linear transformation: x'=(x2'-x1')/(x2-x1) * (x-x1) +x1' - where x,x' are two coordinate systems
    aX=(Allignment(2,2)-Allignment(1,2))/(m2-m1);
    aY=(Allignment(2,3)-Allignment(1,3))/(n2-n1);
    bX=-(Allignment(2,2)-Allignment(1,2))/(m2-m1)*m1+Allignment(1,2);
    bY=-(Allignment(2,3)-Allignment(1,3))/(n2-n1)*n1+Allignment(1,3);
end

if strcmp(NormIntensity,'yes')
    if length(unique(Intensity))>1
        Act_norm=(Intensity-min(Intensity))./max(Intensity-min(Intensity));
    else
        Act_norm=ones(1,length(Intensity));
    end
    if length(unique(ConnectivityIntensity))>1
        Act_norm_con=(ConnectivityIntensity-min(ConnectivityIntensity))./max(ConnectivityIntensity-min(ConnectivityIntensity));
    else
        Act_norm_con=0.5*ones(1,length(ConnectivityIntensity));
    end
    

elseif strcmp(NormIntensity,'posNeg')
    Act_norm=Intensity/2+0.5;
    Act_norm_con=ConnectivityIntensity/2+0.5;
elseif strcmp(NormIntensity,'no')
    Act_norm=Intensity;
    Act_norm_con=ConnectivityIntensity;
else
    error('*****Norm Intensity value not defined properly. Should be ''yes'',''no'',''posNeg''****');
end



%creating translation between electrodes and locations
translation=[];
%En2=fliplr(En);
Elecs=sort(En(~isnan(En)));
for i=1:length(Elecs)
    [n,m]=find(En==Elecs(i));
    translation=[translation;m n Elecs(i)];
end
W=0.01*(XL(2)-XL(1)+YL(2)-YL(1))/2; %Width

if DrawGrid
    for i=1.5:(max_grid_y-0.5)
        TransparentLine([XL(1) XL(2)],[i*aY+bY i*aY+bY],[.7 .7 .7],W,Transparancy);
        %line([XL(1) XL(2)],[i*aY+bY i*aY+bY],'LineWidth',2,'Color',[.1 .5 .7]);
    end
    for i=1.5:(max_grid_x-0.5)
        TransparentLine([i *aX+bX i *aX+bX],[YL(1) YL(2)],[.7 .7 .7],W,Transparancy);
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
for i=1:length(Channels)
    x(i)=aX*translation(find(translation(:,3)==Channels(i)),1)+bX;
    y(i)=aY*translation(find(translation(:,3)==Channels(i)),2)+bY;
    if ~isempty(Intensity) && ~isnan(Act_norm(i))
        hNodes(i)=TransparentCircle(x(i),y(i),CircleSizeFactor*W,cm(ceil(1+Act_norm(i)*63),:),Transparancy);
        %plot(x(i),y(i),'.','color',cm(ceil(1+Act_norm(i)*63),:),'markersize',40);
    end
end

if ~isempty(Connectivity)
    dx=aX*NoiseFac*rand(2,size(Connectivity,2));
    dy=aX*NoiseFac*rand(2,size(Connectivity,2));
    for i=1:size(Connectivity,2)
        %         line([x(Connectivity(1,i))+dx(1,i) x(Connectivity(2,i))+dx(2,i)],[y(Connectivity(1,i))+dy(1,i) y(Connectivity(2,i))+dy(2,i)],...
        %             'color',cm(ceil(1+Act_norm_con(i)*63),:),'LineWidth',3,'Marker','o','MarkerSize',2);
        hEdges(i)=TransparentCircledLine([x(Connectivity(1,i))+dx(1,i) x(Connectivity(2,i))+dx(2,i)],[y(Connectivity(1,i))+dy(1,i) y(Connectivity(2,i))+dy(2,i)],...
            cm(ceil(1+Act_norm_con(i)*63),:),LineSizeFactor*W,Transparancy);
    end
end

%Plotting numbers on top of the propagation points
if DrawElectrodeNumbers
    for i=1:length(Elecs)
        if Elecs(i)<10
            text(aX*(translation(i,1)-0.15)+bX,aY*(translation(i,2))+bY,num2str(Elecs(i)),'fontsize',ElecNumberSize);
        else
            text(aX*(translation(i,1)-0.29)+bX,aY*(translation(i,2))+bY,num2str(Elecs(i)),'fontsize',ElecNumberSize);
        end
    end
end

colormap(cm);
if plotColorBar
    hCbar=colorbar;
    if ~exist('Allignment','var') || isempty(Allignment)
        set(hCbar,'ytick',[eps 1-eps]);
    else
        set(hCbar,'ytick',[eps 63.5]);
    end
    if strcmp(NormIntensity,'no')
        set(hCbar,'yticklabel',[0 1]);
    elseif strcmp(NormIntensity,'PosNeg')
        set(hCbar,'yticklabel',[-1 1]);
    elseif ~isempty(ConnectivityIntensity)
        set(hCbar,'yticklabel',[min(ConnectivityIntensity) max(ConnectivityIntensity)]);
    else
        set(hCbar,'yticklabel',[min(Intensity) max(Intensity)]);
    end
    Str=get(hCbar,'YTickLabel');
    set(hCbar,'YTickLabel',Str(:,1:min(5,size(Str,2))));
end

hold off;

function [Connectivity,ConnectivityIntensity]=connectivityMatrix2ConnectivityIntensity(C)
%find off diagonal element from one side only
NElectrodes=size(C,1);
pDiag=find(tril(ones(NElectrodes),-1)==0);
C(pDiag)=Inf;
[XP,YP]=find(~isnan(C) & ~isinf(C));
Connectivity=[YP';XP'];
ConnectivityIntensity=C(sub2ind([NElectrodes NElectrodes],XP,YP))';

