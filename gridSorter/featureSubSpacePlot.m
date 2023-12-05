function [hFigure]=featureSubSpacePlot(spikeFeatures,idx,varargin)
hFigure=[];
cmap=lines;

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

%create figure if no figure handle given in input
if isempty(hFigure)
    hFigure=figure('Position',[100 100 1200 800],'color','w');
end

%plot
[nSpikes,nFeatures]=size(spikeFeatures);
set(hFigure,'DefaultAxesColorOrder',cmap);
xlimFeatures=zeros(nFeatures,2);
for j=1:nFeatures
    subaxis(nFeatures,nFeatures,j,j,'S',0.005,'M',0.005);
    hist(spikeFeatures(:,j),20);
    axis tight;
    xlimFeatures(j,:)=xlim;
    set(gca,'XTickLabel',[],'YTickLabel',[]);
end
for k=1:nFeatures
    for j=(k+1):nFeatures
        h=subaxis(nFeatures,nFeatures,k,j,'S',0.005,'M',0.005);
        scatter(spikeFeatures(:,k),spikeFeatures(:,j),2,cmap(idx,:));%scatter function with cmap color input is slow but used here instead of color index to allign colors to previous plot
        xlim(xlimFeatures(k,:));
        ylim(xlimFeatures(j,:));
        set(gca,'XTickLabel',[],'YTickLabel',[]);
    end
end
