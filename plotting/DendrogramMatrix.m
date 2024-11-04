% [DC,order,clusters,h,Z]=DendrogramMatrix(C,varargin);
%
% Function purpose : Order matrix according to the hierarchical dendrogram algorithm (Euclidean distance)
%
% Function recives :    C - Correlation matrix
%                       varargin - 
%                           toPlotBinaryTree - if true plots the binary tree (default==0)
%                       
% Function give back : DC - Ordered Correlation matrix
%                      order - the new ordering of rows in the ordered matrix, DC=C(order,order);
%                      clusters - the clusters associated with dendrogram devision (cluster numbers go from top to bottom on the tree)
%                      h - figure axes
%                      Z - linkage in dendrogram
%
% Recommended usage: [DC,order]=DendrogramMatrix(C);
% To show matrix use : imagesc(DC); or pcolor(DC); - to show grid
function [DC,order,clusters,h,Z]=DendrogramMatrix(C,varargin)
%default options
toPlotBinaryTree=0;
figureHandle=[];
linkMethod='ward';
maxClusters=6;
cMap=lines(64);
cLim=[];
treeDepth=2;
clusteringCriterion='inconsistent'; %or 'distance'
cutoff=[];
plotOrderLabels=true;
linkMetric='euclidean';
Orientation='left';
change2SequentialColorOrder=true;
plotRectanglesOnDendrogram=true;
hDendro=[]; %a 1x2 handle array, the first for the tree and the second for the matrix

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
if isempty(figureHandle)
    figureHandle=figure;
else
    if isempty(hDendro)
        figure(figureHandle);
    end
end

if size(C,1)~=size(C,2)
    disp('Analyzing non-square matrix!!!!!!!!!!!!');
    nonSquareMatrix=1;
else
    nonSquareMatrix=0;
    disp('Analyzing square matrix');
end

if any(isnan(C))
    warning('NaN values appear in the activity matrix. This may effect the clustering results!!!');
end
%calculate euclidean distance
Y=pdist(C,linkMetric);
Z=linkage(Y,linkMethod);
%C = cophenet(Z,Y);
%I=inconsistent(Z,4);

if strcmp(Orientation,'right')
    pO1=[1 3];
    pOP={3,1:2};
elseif strcmp(Orientation,'left')
    pO=[1 3];
    pOP={1,2:3};
elseif strcmp(Orientation,'top')
    pO2=[3 1];
    pOP={1,2:3};
elseif strcmp(Orientation,'bottom')
    pO=[3 1];
    pOP={3,1:2};
end


if isempty(hDendro)
    hDendro=subplot(pO(1),pO(2),pOP{1});
else
    hDendro(1).Visible='on';
    axes(hDendro(1));
end

if isempty(cutoff)
    clusters = cluster(Z, 'maxclust',maxClusters);
    [h0,T,order] = dendrogram2(Z,0,'Orientation',Orientation,'colorthreshold',mean(Z(end-maxClusters+1:end-maxClusters+2,3)),'cMap',cMap);
elseif strcmp('clusteringCriterion','distace')
    clusters = cluster(Z, 'cutoff',cutoff,'depth',treeDepth,'criterion','distance');
    [h0,T,order] = dendrogram2(Z,0,'Orientation',Orientation,'ColorThreshold',cutoff,'cMap',cMap);
else
    clusters = cluster(Z, 'cutoff',cutoff,'depth',treeDepth);
    [h0,T,order] = dendrogram2(Z,0,'Orientation',Orientation,'cMap',cMap);
end

if ~plotOrderLabels
    hDendro(1).YTickLabel=[];
end

if change2SequentialColorOrder
    hLines=findobj(hDendro,'Type','Line');
    lineColors=cell2mat(get(hLines,'color'));
    YVal=cell2mat(get(hLines,'YData'));
    inMap=unique(lineColors,'rows');
    inMap=inMap(2:end,:); %remove black color lines
    finalClusters=min(maxClusters,size(inMap,1));
    for i=1:finalClusters
        p{i}=find(lineColors(:,1)==inMap(i,1) & lineColors(:,2)==inMap(i,2) & lineColors(:,3)==inMap(i,3));
        maxVal(i)=max(YVal(p{i},1));
    end
    [~,colorOrder]=sort(maxVal);
    for i=1:finalClusters
        set(hLines(p{colorOrder(i)}),'color',cMap(i,:));
    end
    
    %change cluster number to be sequential
    [~,IA]=unique(clusters(order));
    [~,clusterOrder]=sort(IA);
    tmpClusters=clusters;
    for i=1:finalClusters
        clusters(tmpClusters==clusterOrder(i))=i;
    end
end



if nonSquareMatrix
    DC=C(order,:);
else
    DC=C(order,order);
end

hRect=[];
if ~toPlotBinaryTree
    close(figureHandle);
    hDendro=[];cb=[];
else
    if strcmp(Orientation,'left') || strcmp(Orientation,'bottom')
        set(hDendro(1),'YDir','Reverse');
    end
    
    if isempty(hDendro) | numel(hDendro)==1
        hDendro(2)=subplot(pO(1),pO(2),pOP{2});
    else
        axes(hDendro(2));
    end
    if isempty(cLim)
        imagesc(DC);
    else
        imagesc(DC,cLim);
    end
    cb=colorbar;
    if plotRectanglesOnDendrogram
        hold on;
        ordClu=clusters(order);
        pStart=[1 find(diff(ordClu)==1)'+1];
        pEnd=[pStart(2:end)-1 numel(clusters)];
        for i=1:finalClusters
            hRect(i)=rectangle('position',[pStart(i) pStart(i) pEnd(i)-pStart(i) pEnd(i)-pStart(i)],'edgeColor',cMap(i,:),'lineWidth',3);
        end
    end
    %linkaxes(hDendro,'y');
end
h=[hDendro cb hRect];