

function [fig mx mn] = PlotRawDataNP(obj,opts,params)

arguments (Input) 
    obj 
    opts.fig (1,1) matlab.ui.Figure
    opts.chan double 
    opts.startTimes double
    opts.window double
    opts.spikeTimes double
    params.freq = "AP"
    params.typeD = "line"
    params.multFactor = 1
    params.stdMult = 2
    params.plotSpikeLine = false

end

tic
NP = obj.dataObj;
raw_signal = squeeze(NP.getData(opts.chan,opts.startTimes,opts.window));
toc
%raw_signal = squeeze(aaa(:,1,:));

%aaa = raw_signal;

fc = 300;

if params.freq == "AP"
    [b, a] = butter(4, fc/(NP.samplingFrequency/2), 'high');
end

t = (0:length(raw_signal)-1) /NP.samplingFrequency;
offset = mean(params.stdMult*std(raw_signal));

if size(opts.spikeTimes,1) > size(opts.spikeTimes,2)
    opts.spikeTimes = opts.spikeTimes';
end

if size(raw_signal,1) > size(raw_signal,2)
    raw_signal = raw_signal';
end

fig = figure(opts.fig);

if  params.typeD == "line"

    mx = -inf;
    mn = inf;

for tr = 1:size(raw_signal,1)

    if params.freq == "AP"

        if size(raw_signal,2) == 1
            raw_signal =  raw_signal'; 
        end
        signal = params.multFactor.*filtfilt(b, a, raw_signal(tr,:));

    else
        if size(raw_signal,2) == 1
            raw_signal =  raw_signal'; 
        end
        signal = params.multFactor.*raw_signal(tr,:);
    end

    j = size(raw_signal,1)-tr;
    plot(t, signal+offset*(j), 'LineWidth',0.5,'Color','k');

    if size(opts.spikeTimes,2) > 0

        %Convert ms of raster to seconds
        spikeTimesIndex = find( opts.spikeTimes(tr,:)>0)/1000;
        spikeSamples = round(spikeTimesIndex * NP.samplingFrequency);

        hold on
        %plot(spikeTimesIndex,repmat(min(signal-offset*(j)),1,length(spikeTimesIndex)),'.','Color','b','MarkerSize',7);
        win = round(0.001* NP.samplingFrequency); %ms
        % Overlay red segments for each spike

        %if size(raw_signal,1) >1
        for i = 1:numel(spikeSamples)
            idx1 = max(1, spikeSamples(i)-win);
            idx2 = min(length(signal), spikeSamples(i)+win);
            plot(t(idx1:idx2), signal(idx1:idx2)+offset*(j), 'r', 'LineWidth', 1.5);
        end

        %else
        %%Improve with plot function for several trials
        if params.plotSpikeLine
            try
                xline(spikeTimesIndex,'LineWidth',1.5,'Color','b','Alpha',0.3) %Plot spikes.
            catch
                fprintf('Selected trial has no spikes')
            end
        end
        %end
        
    end
    hold on

 
    if max(signal) > mx
        mx = max(signal);
    end

    if min(signal) < mn
        mn = min(signal);
    end
    hold on
end

%yticks([])

xlim([0 length(raw_signal)/NP.samplingFrequency]);
lims = xlim;

if tr >1
    limsY = ylim;
    ylim([limsY(1)+std(raw_signal(:)) limsY(2)-std(raw_signal(:))]);
end

% ax = gca; % Get current axes
% ax.YAxis.FontSize = 7; % Change font size of y-axis tick labels

hold off

end
% cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')
% print(gcf, sprintf('%s-10TRRawData-MovBall-%s-U%d-W%d-%dW-speed-500.pdf',NP.recordingName,orderNames{k},u,window_size(1),window_size(2)), '-dpdf', '-r300', '-vector');

end