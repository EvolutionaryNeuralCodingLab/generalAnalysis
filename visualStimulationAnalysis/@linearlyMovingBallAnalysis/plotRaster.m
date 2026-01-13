function plotRaster(obj,params)

arguments (Input) 
    obj
    params.overwrite logical = false
    params.analysisTime = datetime('now')
    params.inputParams = false
    params.preBase = 200
    params.bin = 15
    params.exNeurons = 1
    params.AllSomaticNeurons = false
    params.AllResponsiveNeurons = false
    params.SelectedWindow = true
    params.speed string = "max" %or "1", "2", etc
    params.MergeNtrials =1
    params.oneTrial = false
    params.GaussianLength = 10
    params.MaxVal_1 =false
    params.useNormTrialWindow = false
    params.OneDirection string = "all"
    params.OneLuminosity string = "all"
    
end

NeuronResp = obj.ResponseWindow;
Stats = obj.ShufflingAnalysis;

if params.speed ~= "max"
    fieldName =  sprintf('Speed%d',  str2double(params.speed));
else
    fn = fieldnames(NeuronResp);

    % Extract numbers from field names
    nums = nan(numel(fn),1);
    for i = 1:numel(fn)
        tok = regexp(fn{i}, '\d+', 'match');
        if ~isempty(tok)
            nums(i) = str2double(tok{end}); % use last number if multiple
        end
    end

    % Find field with highest number
    [~, idx] = max(nums);
    fieldName = fn{idx};

end

goodU = NeuronResp.goodU;
p = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
phy_IDg = p.phy_ID(string(p.label') == 'good');
pvals = Stats.(fieldName).pvalsResponse;
stimDur = NeuronResp.(fieldName).stimDur;
stimInter = NeuronResp.stimInter;

C = NeuronResp.(fieldName).C;

if params.OneDirection ~= "all"
    switch params.OneDirection
        case "up"
            C = NeuronResp.(fieldName).C(round(C(:,2), 2)==0,:);
        case "left"
            C = NeuronResp.(fieldName).C(round(C(:,2), 2)==1.57,:);
        case "down"
            C = NeuronResp.(fieldName).C(round(C(:,2), 2)==3.14,:);
        case "right"
            C = NeuronResp.(fieldName).C(round(C(:,2), 2)==4.71,:);
        otherwise
            error("Unknown inputPa value: %s", params.OneDirection)
    end
end

if params.OneLuminosity ~= "all"
    switch params.OneLuminosity
        case "black"
            C = NeuronResp.(fieldName).C(round(C(:,6), 2)==1,:);
        case "white"
            C = NeuronResp.(fieldName).C(round(C(:,6), 2)==255,:);
        otherwise
            error("Unknown inputPa value: %s", params.OneLuminosity)
    end
end


[C indexS] = sortrows(C,[2 6 3 4 5]);
directimesSorted = C(:,1)';

%Unique parmeters of the different categories
uDir = unique(C(:,2));
uOffset = unique(C(:,3));
uSize = unique(C(:,4));
uSpeed = unique(C(:,5));
uLums= unique(C(:,6));


%Number of unique parameters per category
offsetN = length(uOffset);
direcN = length(uDir);
sizeN = length(uSize);
speedN = length(uSpeed);
nT = numel(C(:,1));
lumsN = length(uLums);
trialDivision = nT/(offsetN*direcN*speedN*sizeN*lumsN); %Number of trials per unique conditions

preBase = round(stimInter-stimInter/4);

if params.AllSomaticNeurons    
    eNeuron = 1:size(goodU,2);
    pvals = [eNeuron;pvals(eNeuron)];
elseif params.AllResponsiveNeurons
    eNeuron = find(pvals<0.05);
    pvals = [eNeuron;pvals(eNeuron)];% Select all good neurons if not specified  
    if isempty(eNeuron) 
        fprintf('No responsive neurons.\n')
        return
    end
else
    eNeuron = params.exNeurons;
    pvals = [eNeuron;pvals(eNeuron)];
end

[Mr] = BuildBurstMatrix(goodU(:,eNeuron),round(p.t/params.bin),round((directimesSorted-preBase)/params.bin),round((stimDur+preBase*2)/params.bin));

[Mr]=ConvBurstMatrix(Mr,fspecial('gaussian',[1 params.GaussianLength],3),'same');


[nT,nN,nB] = size(Mr);

Mr2 = [];

ur = 1;

for u = eNeuron


    fig =  figure;

    title(sprintf('U.%d-Unit-phy-%d-p-%d',u,phy_IDg(u),pvals(2,ur)));

    %sizeN=1;

    j=1;

    mergeTrials = params.MergeNtrials;

    Mr2 = zeros(size(Mr,1),size(Mr,3));

    if mergeTrials > 1 %Merge trials

        for i = 1:mergeTrials:nT

            meanb = mean(squeeze(Mr(i:min(i+mergeTrials-1, end),ur,:)),1);

            Mr2(i:i+mergeTrials-1,:) = repmat(meanb,[mergeTrials 1]);

            j = j+1;

        end
    else
        Mr2 = Mr(:,ur,:);
    end

    [nT,nN,nB] = size(Mr2);

    if sum(Mr2,'all') ==0
        close
        ur = ur+1;
        continue
    end


    subplot(20,1,[7 18]);
    imagesc(squeeze(Mr2).*(1000/params.bin));colormap(flipud(gray(64)));
    %Plot stim start:
    xline(preBase/params.bin,'k', LineWidth=1.5)
    %Plot stim end:
    xline(stimDur/params.bin+preBase/params.bin,'k',LineWidth=1.5)


    if params.MaxVal_1
        caxis([0 1])
    end
    dirStart = C(1,2);
    offStart = C(1,3);
    lumStart = C(1,6);
    sizeStart = C(1,4);
    for t = 1:nT
        if dirStart ~= C(t,2)
            yline(t-0.5,'k',LineWidth=2);
            dirStart = C(t,2);
        end
        if offStart ~= C(t,3)
            yline(t-0.5,'k',LineWidth=0.5);
            offStart = C(t,3);
        end
        if lumStart ~= C(t,6)
            yline(t-0.5,'--b',LineWidth=1);
            lumStart = C(t,6);
        end
        if sizeStart ~= C(t,4)
            yline(t-0.5,'--r',LineWidth=0.05);
            sizeStart = C(t,4);
        end

    end

    hold on


    xticklabels([])
    xlim([0 round(stimDur+preBase*2)/params.bin])
    xticks([0 preBase/params.bin:300/params.bin:(stimDur+preBase*2)/params.bin (round((stimDur+preBase*2)/100)*100)/params.bin])
    xticklabels([]);

    yt =[0];
    for d = 1:direcN
        yt = [yt [1:trialDivision*2*sizeN:(nT/direcN)-1+trialDivision*sizeN]+max(yt)+trialDivision-1];

    end

    yt = yt(2:end);
    yticks(yt)

    yticklabels(repmat([trialDivision:trialDivision*2*sizeN:(nT/direcN)-1+trialDivision*sizeN],1,direcN))

    ax = gca; % Get current axes
    ax.YAxis.FontSize = 7; % Change font size of y-axis tick labels


    if params.SelectedWindow  %%Select highest window stim type
        j =1;
        meanMr = zeros(1,nT/trialDivision);
        for i = 1:trialDivision:nT
            meanMr(j) = mean(Mr2(i:i+trialDivision-1,:),'all');
            j = j+1;
        end
        
        %Find max trial category
        [maxTrialCat,maxRespIn]= max(meanMr);

        X = Mr2(maxRespIn*trialDivision+1:maxRespIn*trialDivision + trialDivision,:);
        window = 500;
        % Moving mean across 2nd dimension
        mm = movmean(X, window, 2, 'Endpoints', 'discard');
        % Average across rows to get kernel score
        score = mean(mm, 1);
        % Find max kernel location
        [maxVal, idx] = max(score);

        % Kernel column range
        start = idx;

    else
        if params.useNormTrialWindow
            [maxResp,maxRespIn]= max(NeuronResp.(fieldName).NeuronVals(u,:,1));
        else
            [maxResp,maxRespIn]= max(NeuronResp.(fieldName).NeuronVals(u,:,4));
        end
        start = NeuronResp.(fieldName).NeuronVals(u,maxRespIn,3)*NeuronResp.params.binRaster-20;  
        window = 500;
    end

    maxRespIn = maxRespIn-1;
    trials = maxRespIn*trialDivision+1:maxRespIn*trialDivision + trialDivision;
    y1 = maxRespIn*trialDivision + trialDivision;
    y2 = maxRespIn*trialDivision;


  
    patch([(preBase+start)/params.bin (preBase+start+window)/params.bin (preBase+start+window)/params.bin (preBase+start)/params.bin],...
        [y2 y2 y1 y1],...
        'b','FaceAlpha',0.3,'EdgeColor','none')

    patch([1 (preBase*2+stimDur)/params.bin (preBase*2+stimDur)/params.bin 1],...
        [y2 y2 y1 y1],...
        'b','FaceAlpha',0.2,'EdgeColor','b')


    %%%%%% Plot PSTH
    %%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subplot(20,1,[19 20])

    MRhist = BuildBurstMatrix(goodU(:,u),round(p.t),round((directimesSorted-preBase)),round((stimDur+preBase*2)));

    MRhist = squeeze(MRhist(trials,:,:));

    [nT2,nB2]= size(MRhist);

    spikeTimes = repmat([1:nB2],nT2,1);

    spikeTimes = spikeTimes(logical(MRhist));
    % Define bin edges (adjust for resolution)
    binWidth = 125; % 10 ms bins
    if nB>300
         binWidth = 250;
    end
    edges = [1:binWidth:round((stimDur+preBase*2))]; % Adjust time window as needed

    % Compute histogram
    psthCounts = histcounts(spikeTimes, edges);

    % Convert to firing rate (normalize by bin width)
    psthRate = (psthCounts / (binWidth * nT2))*1000;

    b=bar(edges(1:end-1), psthRate,'histc');
    b.FaceColor = 'b';
    b.FaceAlpha = 0.3;
    b.MarkerEdgeColor = "none";
    ylabel('Spikes/sec','FontSize',10);
    xlabel('Seconds','FontSize',10);
    xlim([0 round((stimDur+preBase*2)/100)*100])

    try %zero spiking in selection
        ylim([0 max(psthRate)+std(psthRate)])
    catch
        ur = ur+1;
        close
        continue
    end

    xticks([0 preBase:300:(stimDur+preBase*2) round((stimDur+preBase*2)/100)*100])

    xline([preBase stimDur+preBase],'LineWidth',1.5)

    xticklabels([-(preBase) 0:300:round((stimDur/100))*100 round((stimDur/100))*100 + preBase]./1000)

    ax = gca; % Get current axes
    ax.YAxis.FontSize = 7; % Change font size of y-axis tick labels
    ax.XAxis.FontSize = 7;


    %%%%PLot raw data several trials one
    %%%%channel%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Mark selected trial
   
    bin3 = 2;
    trialM = BuildBurstMatrix(goodU(:,u),round(p.t/bin3),round((directimesSorted+start)/bin3),round((window)/bin3));
    TrialM = squeeze(trialM(trials,:,:))';
    
    [mxTrial TrialNumber] = max(sum(TrialM));

    RasterTrials = trials(TrialNumber);

    chan = goodU(1,u);

    subplot(20,1,[1 3])

    startTimes = directimesSorted(RasterTrials)+start;

    freq = "AP"; %or "LFP"

    typeData = "line"; %or heatmap

    spikes = squeeze(BuildBurstMatrix(goodU(:,u),round(p.t),round(startTimes),round((window))));


    [fig, mx, mn] = PlotRawDataNP(obj,fig = fig,chan = chan, startTimes = startTimes,...
        window = window,spikeTimes = spikes);

    xlabel(string(chan))
    xline(-start/1000,'LineWidth',1.5)
    xticklabels([])
    xlabel([]);xticks([])
    ylabel('uV')
    title(sprintf('U.%d-Unit-phy-%d-p-%d',u,phy_IDg(u),pvals(2,ur)));

    %%%%%%%%%%% Plot raster of selected trials
    %%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(20,1,[4 5])

    RasterTrials =  trials;

    bin3 = 2;
    trialM = BuildBurstMatrix(goodU(:,u),round(p.t/bin3),round((directimesSorted+start)/bin3),round((window)/bin3));

    if numel(RasterTrials) == 1
        TrialM = squeeze(trialM(RasterTrials,:,:))';
    else
        TrialM = squeeze(trialM(trials,:,:));
    end

    TrialM(TrialM~=0) = 0.3;
    spikes1 = TrialM(TrialNumber,:);
    spikeLoc = find(spikes1 >0);
    if isempty(spikeLoc)
        close
        ur = ur+1;
        continue
    end
    TrialM(TrialNumber,spikeLoc) = 1;

    %Select offset in which selected trial belongs (10 trials)
    imagesc(TrialM);colormap(flipud(gray(64)));
    %caxis([0 1])
    xline([preBase stimDur+preBase],'LineWidth',1.5)
    ylabel([sprintf('%d trials',numel(trials))])
    if params.fixedWindow
        xticks([1 abs(start)/bin3:20/bin3:window/bin3])
        xticklabels([start 0:20:window-abs(start)])
    else
        xticks([1:20/bin3:window/bin3])
        xticklabels([start:20:window+abs(start)])
    end

    xlabel('Milliseconds')
    set(gca,'FontSize',7)

    xline(-start/bin3,'LineWidth',1.5)

    xline(spikeLoc,'LineWidth',1,'Color','r','Alpha',0.3)

    yline(TrialNumber,'LineWidth',3,'Color','r','Alpha',0.3) %Mark trial

    fig.Position = [680     5   296   9734];

    if params.overwrite,obj.printFig(fig,sprintf('%s-%s-MovBall-SelectedTrials-eNeuron-%d',obj.dataObj.recordingName,fieldName,u)),end

    if ur ~= length(eNeuron)
        close
    end
    
    ur = ur+1;

end %end eNeuron for loop

end %end plotRaster