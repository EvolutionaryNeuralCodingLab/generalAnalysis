function [hFigure]=PCAfeaturesSpikeShapePlot(spikeFeatures,spikeShapes,Fs,idx1,idx2,En,ch,varargin)
removeOutliers=0;
hFigure=[];
cmap=lines;
rejectPCAStd=20; %the multiplicator in PCA space for determining a thershold for outlier rejection
maxNumberOutliers=50;

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

if isempty(hFigure)
    hFigure=figure('Position',[100 100 1200 800],'color','w');
end
set(hFigure,'DefaultAxesColorOrder',cmap);

hHBox = uix.HBox( 'Parent', hFigure, 'Padding', 5, 'Spacing', 1 );
hVBox = uix.VBox( 'Parent', hHBox, 'Padding', 5, 'Spacing', 1 );
hPanel = uipanel( 'Parent', hHBox, 'BorderType','none','BackgroundColor','w');
set( hHBox, 'Widths', [-1 -2] );
h1=axes('Parent',hVBox);
h2=axes('Parent',hVBox);
set( hVBox, 'Heights', [-1 -1] );

[PCAFeatureSimMat,PCAspikeFeatures] = pca(spikeFeatures); %run PCA for visualization purposes
sPC1=median(abs(  PCAspikeFeatures(:,1)-median(PCAspikeFeatures(:,1))  )) / 0.6745;
sPC2=median(abs(  PCAspikeFeatures(:,2)-median(PCAspikeFeatures(:,2))  )) / 0.6745;
pOutliers=find( PCAspikeFeatures(:,1)>rejectPCAStd*sPC1 | PCAspikeFeatures(:,1)<-rejectPCAStd*sPC1 | PCAspikeFeatures(:,2)>rejectPCAStd*sPC2 | PCAspikeFeatures(:,2)<-rejectPCAStd*sPC2);
nOutliers=numel(pOutliers);

if nOutliers<maxNumberOutliers && removeOutliers
    disp([num2str(nOutliers) ' outliers removed']);
    PCAspikeFeatures(pOutliers,:)=[];
    spikeFeatures(pOutliers,:)=[];
    nSpikes=nSpikes-nOutliers;
    [PCAFeatureSimMat,PCAspikeFeatures] = princomp(spikeFeatures); %run PCA for visualization purposes
end

axes(h1); %PCA space before matching
cent=[];
nClusters1=numel(unique(idx1));
for k=1:nClusters1
    cent(k,:)=mean(PCAspikeFeatures(idx1==k,1:2),1);
end
h1s=scatter(PCAspikeFeatures(:,1),PCAspikeFeatures(:,2),2,cmap(idx1,:));
text(cent(:,1),cent(:,2),mat2cell(1:nClusters1,1,ones(1,nClusters1)),'FontWeight','Bold');
xlabel('PCA 1');
ylabel('PCA 2');
freezeColors(h1s);

axes(h2); %PCA space after matching
%recalculate centers
cent=[NaN NaN];
nClusters2=numel(unique(idx2));
for k=1:nClusters2
    cent(k,:)=mean(PCAspikeFeatures(idx2==k,1:2));
end
h2s=scatter(PCAspikeFeatures(:,1),PCAspikeFeatures(:,2),2,cmap(idx2,:));
text(cent(:,1),cent(:,2),mat2cell(1:nClusters2,1,ones(1,nClusters2)),'FontWeight','Bold');
xlabel('PCA 1');
ylabel('PCA 2');
freezeColors(h2s);
%set(h2,'position',[0.0533    0.0575    0.3632    0.4238]);

[nSpikeSamples,nSpikes,nOverheadSurroundingCh]=size(spikeShapes);
avgSpikeWaveforms=zeros(nSpikeSamples,max(1,nClusters2),nOverheadSurroundingCh);
for j=1:nClusters2
    avgSpikeWaveforms(:,j,:)=median(spikeShapes(:,idx2==j,:),2);
end

[h,hParent]=spikeDensityPlotPhysicalSpace(spikeShapes,Fs,ch,En,...
    'hParent',hPanel,'avgSpikeWaveforms',avgSpikeWaveforms,'logDensity',true,'cmap',cmap);
