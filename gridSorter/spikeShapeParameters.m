function [par,h]=spikeShapeParameters(WF,preMs,Fs,varargin)
% [par,h]=spikeShapeParameters(WF,preMs,postMs,Fs,varargin)
% Function purpose : Calculate parameters of spike shapes
%
% Function recives :    WF - average spike waveforms over all electrodes in ch [Double [neurons x ch x samples]
%                       preSpikeMs - pre spike time in WF [ms]
%                       postSpikeMs - pre spike time in WF [ms]
%                       Fs - sampling frequency of the WFs
%                       varargin ('property name','property value')
%
% Function give back :  par - a structure of output parameters
%                       h - a structure of handles from generated plots
%
% Last updated : 12/12/14

%default parameters
smoothingDuration=0.5; %[ms]
plotWidth=false;
plotWidthSpkNumbers=false; %plot spike numbers for each point
plotWidthAddNoise=true; %add noise a noise with a max shift of the sampling frequency to eliminate quantization effects in plot
plotWidthColorBar=true; %add color bar
plotWidthLogColormap=true; %plot colorbar in log scale

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
    eval([varargin{i} '=' 'varargin{i+1};']);
end

h=[];
% Find neuron with reliable spiking
[nSpikes,nSpikeSamples]=size(WF);

% calculate pre and post samples and their position in the waveform
preSpkSamples=preMs*Fs/1000;
postSpkSamples=nSpikeSamples-preSpkSamples;
pPre=1:preSpkSamples;
pPost=(preSpkSamples+1):(preSpkSamples+postSpkSamples);


%smooth spike shapes with local linear regression
smoothingSamples=round(smoothingDuration*Fs/1000);
spikeShapeSmooth=zeros([nSpikes,nSpikeSamples]);
for i=1:nSpikes
    spikeShapeSmooth(i,:) = smooth(WF(i,:),smoothingSamples,'loess');
end

%extract spike maxima
par.extermalAmp=max(abs(spikeShapeSmooth),[],2);
par.maxAmp=max(spikeShapeSmooth,[],2);
par.minAmp=min(spikeShapeSmooth,[],2);

%extract width from crossing at half width
spkWidthSamples=zeros(nSpikes,1);
widthVolt=(par.maxAmp+par.minAmp)/2;
for i=1:nSpikes
    pDown=find(spikeShapeSmooth(i,pPre)>widthVolt(i),1,'last');
    pUp=preSpkSamples+find(spikeShapeSmooth(i,pPost)>widthVolt(i),1,'first');
    if ~(isempty(pDown) || isempty(pUp))
        spkWidthSamples(i)=(pUp-pDown);
    else
        spkWidthSamples(i)=NaN;
    end
end
par.width=spkWidthSamples/Fs*1000;

%calculate width from time difference between minimal and maximal derivative
spikeShapeDer=diff(spikeShapeSmooth,[],2);
[~,derMinP]=min(spikeShapeDer(:,pPre),[],2);
[~,derMaxP]=max(spikeShapeDer(:,pPost(1:end-1)),[],2);
par.widthDer=( (preSpkSamples+derMaxP)-derMinP)/Fs*1000;

if plotWidth
    if plotWidthAddNoise
        noiseAdd=(rand(nSpikes,2)-0.5)/Fs*1000;
    else
        noiseAdd=zeros(nSpikes,2);
    end
    
    f=figure;h.widthAxes=axes;
    if plotWidthLogColormap
        h.widthScatter=scatter(par.widthDer+noiseAdd(:,1),par.width+noiseAdd(:,2),25,log(par.extermalAmp));
    else
        h.widthScatter=scatter(par.widthDer+noiseAdd(:,1),par.width+noiseAdd(:,2),25,par.extermalAmp);
    end
    xlabel('Width (derivative)');
    ylabel('Width @ half maxima');
    
    if plotWidthSpkNumbers
        text(par.widthDer+noiseAdd(:,1),par.width+noiseAdd(:,2),num2str((1:nSpikes)'),'fontSize',8);
    end
    
    if plotWidthColorBar
        h.widthCb=colorbar;
        if plotWidthLogColormap
            cbTicks=get(h.widthCb,'YTick');
            set(h.widthCb,'YTickLabel',num2str(exp(cbTicks)',3));
        end
        set(gca,'position',[0.1300    0.1100    0.7164    0.8150]);
        set(h.widthCb,'position',[0.8607    0.1071    0.0179    0.8167]);
        ylabel(h.widthCb,'Extermal Amp [\muV]');
    end

end

