function [h,hParent,hScaleBar]=spikeDensityPlotPhysicalSpace(spikeShapes,Fs,ch,En,varargin)
% [h,hParent,hScaleBar]=spikeDensityPlotPhysicalSpace(spikeShapes,Fs,ch,En,varargin)
% [nSpikeSamples,nSpikes,nCh]
%default argumnets
avgSpikeWaveforms = [];
hParent = [];
lineWidth = 2;
lineColor=[0.85 0.325 0.098];
yl = [];
scaleBar = true;
plotChannelNumbers = true;
logColorScale=true;
cmap=lines(64);

spikeShapesIsAHist=false;
xBin=[];
yBin=[];

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

%Collects all options
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end

%calculate general variables
tBins=[];Vbins=[];
if ~isempty(spikeShapes)
    if ~spikeShapesIsAHist
        [nSpikeSamples,nSpikes,nCh]=size(spikeShapes);
        sampleTimes=((1:nSpikeSamples)/Fs*1000)';
    else
        [nCh,tBins,Vbins]=size(spikeShapes);
        spikeShapes=double(spikeShapes);
    end
else
    [nSpikeSamples,~,nCh]=size(avgSpikeWaveforms);
    sampleTimes=((1:nSpikeSamples)/Fs*1000)';
end

%Remove empty channels from layout (En) to make a nicer plot
p1=all(isnan(En));
p2=all(isnan(En'));
En(:,p1)=[];En(p2,:)=[];

%calculate the surrounding channel grid according to channel layout
En=flipud(En);
for i=1:numel(ch)
    [xCh(i),yCh(i)]=find(En==ch(i));
end
xCh=xCh-min(xCh)+1;
yCh=yCh-min(yCh)+1;
n=max(yCh);
m=max(xCh);

%check if first argument is a valid parent handle, else make one
if isempty(hParent)
    hParent=figure('position',[100 200 1000 800]);
end
%P = panel(hParent);
%P.pack(m,n);
%P.margin=0.005;

%calculate limits for plot
if isempty(yl)
    if ~isempty(avgSpikeWaveforms)
        yl=[min(avgSpikeWaveforms(:)-10) max(avgSpikeWaveforms(:)+10)];
        sampleTimesAvg=((1:size(avgSpikeWaveforms,1))/Fs*1000)';
    else
        yl=[min(spikeShapes(:)) max(spikeShapes(:))];
    end
end

mDens=mean(spikeShapes(:));
sDens=std(spikeShapes(:));
densRange=[max(0,mDens-3*sDens) min(max(spikeShapes(:)),mDens+3*sDens)+eps];
%plot channels
for i=1:nCh
    %h(i)=P(yCh(i),xCh(i)).select();
    h(i)=subaxis(hParent,m,n,yCh(i),xCh(i),'S',0.002,'M',0.002);
    if ~isempty(spikeShapes)
        if spikeShapesIsAHist
            colormap(h(i),flipud(gray(256)));
            if logColorScale
                imagesc(xBin,yBin,log10(1+squeeze(spikeShapes(i,:,:))'),log10(1+densRange));
            else
                imagesc(xBin,yBin,squeeze(spikeShapes(i,:,:))',densRange);
            end
        else
            hist2(sampleTimes*ones(1,nSpikes),squeeze(spikeShapes(:,:,i)),'c1',sampleTimes,'n2',100,'h',h(i),'logColorScale',logColorScale,'plotColorBar',0);
        end
        set(h(i),'YDir','normal');
    end
    hold on;
    if ~isempty(avgSpikeWaveforms)
        plot(sampleTimesAvg,squeeze(avgSpikeWaveforms(:,:,i)),'lineWidth',lineWidth);
    end
    if plotChannelNumbers
        text(sampleTimesAvg(1)+(sampleTimesAvg(end)-sampleTimesAvg(1))*0.1,yl(1)+(yl(2)-yl(1))*0.13,num2str(ch(i)),'fontWeight','bold','fontSize',10);
    end
    ylim(yl);
    xlim(sampleTimesAvg([1 end]));
    set(h(i),'YTickLabel',[]);
    set(h(i),'XTickLabel',[]);
    hold off;
end

%add scale bar to plot
if scaleBar
    [yCorder]=max(yCh(xCh==m));
    [~,pCh]=find(xCh==m & yCh==yCorder);
    if yCorder<n
        [hScaleBar]=addScaleBar(h(pCh),'scaleFac',8,'xShift',1);
    else
        [hScaleBar]=addScaleBar(h(pCh),'scaleFac',5);
    end
end