function plotRaster(obj,params)

arguments (Input) 
    obj
    params.overwrite logical = false
    params.analysisTime = datetime('now')
    params.inputParams = false
    params.preBase = 200
    params.bin = 30
    params.exNeurons = 1
    params.exNeuronsPhyID double = []   % alternative to exNeurons: specify neurons by phy cluster ID
    params.AllSomaticNeurons = false
    params.AllResponsiveNeurons = false
    params.SelectedWindow = true
    params.speed string = "max" %or "1", "2", etc
    params.MergeNtrials =1
    params.oneTrial = false
    params.GaussianLength = 10
    params.Gaussian logical = false
    params.MaxVal_1 =true
    params.useNormTrialWindow = false
    params.OneDirection string = "all"
    params.OneLuminosity string = "all"
    params.PaperFig logical = false
    params.statType string = "maxPermuteTest"
    params.sortingOrder string = ["direction","luminosity","offset","size","speed"]

    
end

NeuronResp = obj.ResponseWindow;

if params.statType == "BootstrapPerNeuron"
    Stats = obj.BootstrapPerNeuron;
else
    Stats = obj.StatisticsPerNeuron;
end



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

p = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
phy_IDg = p.phy_ID(string(p.label') == 'good');
pvals = Stats.(fieldName).pvalsResponse;
stimDur = NeuronResp.(fieldName).stimDur;
stimInter = NeuronResp.stimInter;
label = string(p.label');
goodU = p.ic(:,label == 'good'); %somatic neurons

% Convert phy IDs to unit indices if exNeuronsPhyID is provided.
% This overrides exNeurons if both are set — phy ID is more explicit.
if ~isempty(params.exNeuronsPhyID)
    [found, neuronIdx] = ismember(params.exNeuronsPhyID, phy_IDg);
    if any(~found)
        warning('The following phy IDs were not found in good units and will be skipped: %s', ...
            num2str(params.exNeuronsPhyID(~found)));
    end
    params.exNeurons = neuronIdx(found);  % convert to regular indices
    fprintf('  Converted phy IDs [%s] -> unit indices [%s]\n', ...
        num2str(params.exNeuronsPhyID(found)), num2str(params.exNeurons));
end

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
            C = C(round(C(:,6), 2)==1,:);
        case "white"
            C = C(round(C(:,6), 2)==255,:);
        otherwise
            error("Unknown inputPa value: %s", params.OneLuminosity)
    end
end


sortMap = struct( ...
    'direction', 2, ...
    'offset',      3, ...
    'size',       4,...
    'speed', 5, ...
    'luminosity',6);

sortCols = zeros(1,numel(params.sortingOrder));

for k = 1:numel(params.sortingOrder)
    sortCols(k) = sortMap.(params.sortingOrder(k));
end

[C, sortIdx] = sortrows(C, sortCols);

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

%Calculate size of ball in degrees:
%Standard measurements for last set of experiments:
eye_to_monitor_distance = 21.5000;
pixel_size = 33;
resolution = 1080;
pixel_size =  pixel_size/resolution;
monitor_resolution = [1920 1080];
[theta_x,theta_y] = pixels2eyeDegrees(eye_to_monitor_distance,pixel_size,monitor_resolution);

for i = 1:sizeN
    sizeBall(i) = round(abs(abs(theta_x(1,uSize(i)))-abs(theta_x(1,1))),2);
end

sizesString = strjoin(string(sizeBall), "_");

% ==========================================================
% OUTER GROUP INFO
% ==========================================================

outerGroup = params.sortingOrder(1);
outerCol   = sortMap.(outerGroup);

outerVals = C(:,outerCol);

uOuter = unique(outerVals,'stable');

groupStarts = zeros(numel(uOuter),1);
groupEnds   = zeros(numel(uOuter),1);

for i = 1:numel(uOuter)

    idx = find(outerVals == uOuter(i));

    groupStarts(i) = idx(1);
    groupEnds(i)   = idx(end);

end

groupCenters = (groupStarts + groupEnds)/2;

% ==========================================================
% BUILD GROUP LABELS
% ==========================================================

groupLabels = strings(size(uOuter));

switch outerGroup

    % ======================================================
    case "direction"

        for i = 1:numel(uOuter)

            val = round(uOuter(i),2);

            switch val

                % Direction group labels use actual arrow glyphs (Unicode arrow
                % code points, via char() for encoding-independent rendering)
                % instead of letter codes, so the group label reads visually as
                % the ball's motion direction rather than as text.
                case 0
                    groupLabels(i) = char(8593); % U+2191 up arrow

                case 1.57
                    groupLabels(i) = char(8592); % U+2190 left arrow

                case 3.14
                    groupLabels(i) = char(8595); % U+2193 down arrow

                case 4.71
                    groupLabels(i) = char(8594); % U+2192 right arrow

                otherwise
                    groupLabels(i) = string(val);

            end
        end

    % ======================================================
    case "size"

        groupLabels = strings(size(uOuter));

        for i = 1:numel(uOuter)

            idxSize = find(uSize == uOuter(i));

            if ~isempty(idxSize)
                groupLabels(i) = sprintf('%.1f°', sizeBall(idxSize));
            else
                groupLabels(i) = string(uOuter(i));
            end

        end

    % ======================================================
    case "luminosity"

        for i = 1:numel(uOuter)

            if uOuter(i) == 1
                groupLabels(i) = "B";
            elseif uOuter(i) == 255
                groupLabels(i) = "W";
            else
                groupLabels(i) = string(uOuter(i));
            end

        end

    % ======================================================
    case "offset"

        for i = 1:numel(uOuter)
            groupLabels(i) = sprintf('O%d',uOuter(i));
        end

    % ======================================================
    case "speed"

        for i = 1:numel(uOuter)
            groupLabels(i) = sprintf('S%d',uOuter(i));
        end

end




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

if params.Gaussian
    [Mr]=ConvBurstMatrix(Mr,fspecial('gaussian',[1 params.GaussianLength],3),'same');
end

channels = goodU(1,eNeuron);

[nT,nN,nB] = size(Mr);

Mr2 = [];

ur = 1;


for u = eNeuron


    fig =  figure;



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


    subplot(18,1,[6 16]);
    imagesc(squeeze(Mr2).*(1000/params.bin));colormap(flipud(gray(64)));
    set(gca,'Clipping','off')
    %Plot stim start:
    xline(preBase/params.bin,'k', LineWidth=1.5)
    %Plot stim end:
    xline(stimDur/params.bin+preBase/params.bin,'k',LineWidth=1.5)


    if params.MaxVal_1
        caxis([0 1])
    end
    % ==========================================================
    % DYNAMIC GROUP LINES
    % ==========================================================

    hold on

    % --- outer group (thick lines) ---
    for i = 2:numel(groupStarts)

        % LineWidth reduced from 2 to 1: at width 2 these category-boundary lines
        % were visually covering part of the first/last trial row on each side of
        % the boundary. Still solid black (vs. the dashed gray inner-group lines
        % below) so the outer/inner hierarchy stays visually distinct.
        yline(groupStarts(i)-0.5, ...
            'b', ...
            'LineWidth',1);

    end

    % ==========================================================
    % INNER GROUPS
    % ==========================================================

    innerOrders = params.sortingOrder(2:end);

    lineStyles = {'-','--',':'};
    lineWidths = [0.8 0.8 0.8];

    for s = 1:numel(innerOrders)

        thisGroup = innerOrders(s);

        thisCol = sortMap.(thisGroup);

        vals = C(:,thisCol);

        prevVal = vals(1);

        for t = 2:nT

            if vals(t) ~= prevVal

                yline(t-0.5, ...
                    'Color',[0.4 0.4 0.4], ...
                    'LineStyle',lineStyles{min(s,numel(lineStyles))}, ...
                    'LineWidth',lineWidths(min(s,numel(lineWidths))));

                prevVal = vals(t);

            end

        end
    end

    % ==========================================================
    % GROUP LABELS
    % ==========================================================

    xText = -8;

    % Direction group labels (now arrow glyphs, see the "direction" case above)
    % are placed to the right of the raster instead of the left; every other
    % group-label type (size/luminosity/offset/speed) keeps the original left
    % placement at xText.
    if outerGroup == "direction"
        xTextPos = round(stimDur+preBase*2)/params.bin + 3;  % just past the raster's right edge
        horzAlign = 'left';
        labelFontSize = 14;  % larger than other group labels so the arrow glyphs read as bold/thick
    else
        xTextPos = xText;
        horzAlign = 'right';
        labelFontSize = 8;
    end

    for i = 1:numel(groupCenters)

        text(xTextPos, ...
            groupCenters(i), ...
            groupLabels(i), ...
            'HorizontalAlignment',horzAlign, ...
            'VerticalAlignment','middle', ...
            'FontWeight','bold', ...
            'FontSize',labelFontSize, ...
            'FontName','helvetica', ...
            'Clipping','off');

    end
    hold on

    xticklabels([])
    xlim([0 round(stimDur+preBase*2)/params.bin])
    xticks([0 preBase/params.bin:600/params.bin:(stimDur+preBase*2)/params.bin (round((stimDur+preBase*2)/100)*100)/params.bin])
    xticklabels([]);

    yt =[0];
    for d = 1:direcN
        yt = [yt [1:trialDivision*2*sizeN:(nT/direcN)-1+trialDivision*sizeN]+max(yt)+trialDivision-1];

    end

    yt = yt(2:end-1);
    yticks(yt)

    yticklabels(repmat([trialDivision:trialDivision*2*sizeN:(nT/direcN)-1+trialDivision*sizeN],1,direcN))

    ax = gca; % Get current axes
    ax.YAxis.FontSize = 8; % Change font size of y-axis tick labels
    ax.YAxis.FontName = 'helvetica';

    ylabel('Trials','FontSize',10,'FontName','helvetica')

    if params.SelectedWindow  %%Select the trial set whose best in-stimulus 500 ms window is max
        window = 500; %in ms
        winBins = round(window/params.bin);   % number of bins spanning the fixed extraction window

        % Full merged/averaged raster as (trials x bins). squeeze() handles both
        % layouts: mergeTrials>1 leaves Mr2 already 2D, mergeTrials==1 leaves it
        % trials x 1 x bins. The search runs over ALL trials/categories at once so
        % the winning trial set is chosen by its single best 500 ms window (below),
        % not by overall epoch spiking.
        Xall = squeeze(Mr2);
        Xall(Xall>1) = 1;
        [nAllRows, n_cols] = size(Xall);

        % Restrict the sliding-window search so the selected window lies ENTIRELY
        % within the stimulus period — never in the pre/post baseline. Column c of
        % Xall covers ms offset [(c-1)*bin, c*bin) relative to the epoch start (which
        % is preBase ms before stim onset), so a window starting at column c spans ms
        % [(c-1)*bin, (c-1)*bin + window]. Requiring that interval to sit inside
        % [preBase, preBase+stimDur] gives these column bounds (also clamped to the
        % matrix so a window is never read past the available bins).
        firstStimCol = max(1, ceil(preBase/params.bin) + 1);                                     % earliest start col with window start >= stim onset
        lastStimCol  = min(n_cols - winBins + 1, floor((preBase+stimDur-window)/params.bin) + 1); % latest start col with window end <= stim offset

        % If no window of this fixed width fits inside the stimulus (e.g. stimDur <
        % window, or the epoch is too short), there is no valid in-stimulus
        % selection. This is a recoverable, experiment-level parameter mismatch (not
        % an alignment invariant), so warn explicitly and skip this neuron rather
        % than erroring the whole batch run.
        if lastStimCol < firstStimCol
            warning('plotRaster:WindowExceedsStim: neuron phyID=%d (unit idx %d), experiment %s: no %d ms window (%d bins) fits inside the stimulus (stimDur=%d ms, bin=%d ms). Skipping this neuron.', ...
                phy_IDg(u), u, obj.dataObj.recordingName, window, winBins, stimDur, params.bin);
            close   % discard this neuron's partially-built figure
            ur = ur+1;
            continue   % move to the next neuron in the eNeuron loop
        end

        % Mean over every candidate in-stimulus window, for EVERY trial in EVERY
        % category (the full raster), so the peak window is found jointly across
        % trial sets rather than within a pre-chosen one.
        candCols = firstStimCol:lastStimCol;                 % valid window-start columns (each keeps the whole window inside the stimulus)
        window_means = zeros(nAllRows, numel(candCols));
        for k = 1:numel(candCols)
            col = candCols(k);                               % actual column for this candidate window
            window_means(:, k) = mean(Xall(:, col:col+winBins-1), 2);
        end

        % The global maximum picks the single highest 500 ms window across all trials
        % and positions. The winning trial set (category) is simply whichever block
        % of trialDivision trials that winning row falls in — this replaces the old
        % "highest overall epoch spiking" category selection with "highest peak
        % windowed response".
        [~, linear_idx] = max(window_means(:));
        [best_row_all, best_col_local] = ind2sub(size(window_means), linear_idx);
        best_col  = candCols(best_col_local);                        % window-start column on the real time axis
        maxRespIn = floor((best_row_all-1)/trialDivision);           % 0-indexed winning trial set (category)
        best_row  = best_row_all - maxRespIn*trialDivision;          % winning row within that trial set (1..trialDivision)

        % When mergeTrials>1, Mr2 held trial-averaged data at search time, so best_row
        % reflects an arbitrary "first row of the winning merge-group", not the real
        % trial with the most spikes (rows within a merge-group are numerically
        % identical after averaging, so max() ties break to the first one). Re-derive
        % best_row from the REAL unmerged Mr for this neuron, restricted to the same
        % best_col time window, so the displayed trial is the actual highest-firing
        % one. When mergeTrials==1, Mr2 already equals real per-trial data (see the
        % 'else' branch above), so this re-derivation is skipped.
        if mergeTrials > 1
            trialsForBestCol = maxRespIn*trialDivision+1 : maxRespIn*trialDivision+trialDivision;  % real trial rows in the winning condition-block
            XrealWindow = squeeze(Mr(trialsForBestCol, ur, best_col:best_col+winBins-1));  % real (non-averaged) spike counts, this neuron, in the identified window
            realRowSums = sum(XrealWindow, 2);   % total real spikes per trial row in this window
            [~, best_row] = max(realRowSums);    % pick the real trial with the most spikes
        end

        % Kernel column range.
        % best_col is a 1-indexed BuildBurstMatrix column; its bin covers ms offset
        % [(best_col-1)*bin, best_col*bin) (confirmed in BuildBurstMatrix.c), so the
        % window START is (best_col-1)*bin — using best_col*bin here shifted the
        % extracted window and the red patch one full bin later than the window
        % actually selected by the search above.
        start = (best_col-1)*params.bin;


        % % --- Plot ---
        % figure;
        % imagesc(X);
        % colorbar;
        % axis tight;
        % hold on;
        % 
        % % Highlight the full best row (horizontal span)
        % rectangle('Position', [0.5, best_row - 0.5, n_cols, 1], ...
        %     'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
        % 
        % % Highlight the selected window (column span within best row)
        % rectangle('Position', [best_col - 0.5, best_row - 0.5, round(window/params.bin), 1], ...
        %     'EdgeColor', 'y', 'LineWidth', 2.5);

    else
        if params.useNormTrialWindow
            [maxResp,maxRespIn]= max(NeuronResp.(fieldName).NeuronVals(u,:,1));
        else
            [maxResp,maxRespIn]= max(NeuronResp.(fieldName).NeuronVals(u,:,4));
        end
        start = NeuronResp.(fieldName).NeuronVals(u,maxRespIn,3)*NeuronResp.params.binRaster-20;  
        window = 500;
        maxRespIn = maxRespIn-1;
    end

    trials = maxRespIn*trialDivision+1:maxRespIn*trialDivision + trialDivision;
    y1 = maxRespIn*trialDivision + trialDivision+0.5;
    y2 = maxRespIn*trialDivision+0.5;

  
    %patch([(preBase+start)/params.bin (preBase+start+window)/params.bin (preBase+start+window)/params.bin (preBase+start)/params.bin],...
        % [y2 y2 y1 y1],...
        % 'r','FaceAlpha',0.2,'EdgeColor','none')

    patch([0 (preBase*2+stimDur)/params.bin (preBase*2+stimDur)/params.bin 0],...
        [y2 y2 y1 y1],...
        'k','FaceAlpha',0.1,'EdgeColor','none')


    % TrialM = squeeze(Mr2(trials,round((preBase+start)/params.bin):round((preBase+start+window)/params.bin)))';
    % 
    % [mxTrial TrialNumber] = max(sum(TrialM));

    RasterTrials = trials(best_row);

    % patch([(preBase+start)/params.bin (preBase+start+window)/params.bin (preBase+start+window)/params.bin (preBase+start)/params.bin],...
    %     [RasterTrials-0.5 RasterTrials-0.5 RasterTrials+0.5 RasterTrials+0.5],...
    %     'r','FaceAlpha',0.3,'EdgeColor','none')

    % A red fill confined to exactly one trial row (tried previously at FaceAlpha
    % 0.3 and 0.6) is essentially unreadable in a dense raster, since one row is
    % only a few screen pixels tall no matter the opacity or z-order — a fixed
    % geometry problem, not a colour problem. Two elements identify the
    % extraction instead: (1) a translucent red vertical band spanning the whole
    % winning trial-category block (y1:y2, already computed above) marks the
    % extracted TIME WINDOW at a height that is actually visible; (2) the exact
    % TRIAL ROW is identified below by recoloring its own real spikes in red.
    hBand = patch([(start)/params.bin (start+window)/params.bin (start+window)/params.bin (start)/params.bin],...
        [y2 y2 y1 y1],...
        'r','FaceAlpha',0.25,'EdgeColor','none');
    uistack(hBand,'top');

    % A plot()-based triangle/arrow marker placed past the raster's right edge
    % did not render reliably. The direction-arrow group labels use text() at
    % that same past-the-edge position and do render correctly (confirmed
    % above), so the same text()-based approach is reused here for a red
    % left-pointing triangle glyph at the exact RasterTrials row — a fixed-size
    % marker that does not depend on how many spikes that trial happens to have.
    % text(round(stimDur+preBase*2)/params.bin + 3, RasterTrials, char(9668), ...   % U+25C4 black left-pointing triangle
    %     'Color','r', 'FontSize',16, 'FontWeight','bold', ...
    %     'HorizontalAlignment','left', 'VerticalAlignment','middle', 'Clipping','off');

    % Recolor in red the raster's OWN displayed bins for the selected row, inside
    % the extraction window. The red MUST come from Mr2 (the matrix actually drawn
    % by imagesc above), not from the real single trial Mr: with MergeNtrials>1 the
    % displayed row is the merge-averaged pattern, so Mr's real-trial spikes sit at
    % different time columns than the black pixels on screen and the red looked
    % time-shifted relative to the raster. Using Mr2 makes the red overlay the
    % exact bins the raster shows, within the window.
    % Columns are derived from start/window (the same quantities passed to the raw
    % extraction), so the red span matches the red band above in both the
    % SelectedWindow and the NeuronVals branches.
    % Bin count comes from size(Mr,3): Mr is always 3D (trials x neurons x bins),
    % whereas Mr2 is 2D only when merged (it is trials x 1 x bins when unmerged, so
    % size(Mr2,2) would be 1 there). The 2-subscript Mr2(row,col) read below is
    % correct in both layouts because Mr2's middle dim is 1 when 3D.
    winColStart = round(start/params.bin) + 1;                                       % first raster column of the extraction window (col c covers [(c-1)*bin, c*bin) ms)
    winCols = winColStart : min(winColStart + round(window/params.bin) - 1, size(Mr,3));  % window columns, clamped to the available bins
    dispRowSpikes = Mr2(RasterTrials, winCols);                         % displayed (raster) values for this row, inside the window
    spikeColsFull = winCols(dispRowSpikes > 0);                         % window bin columns the raster shows as active
    if ~isempty(spikeColsFull)
        % Each spike is drawn as a vertical line rather than a one-row-tall
        % patch. A patch sized to exactly one trial row is invisible here: with
        % many trials in this short subplot a single row is sub-pixel, so it is
        % rasterized away completely (verified: it rendered 0 red pixels at this
        % figure's real geometry). A line's LineWidth is in points (device
        % units), so it always survives rasterization, and the vertical extent is
        % scaled to the trial count to stay visible without swamping neighbours.
        halfH = max(0.5, size(Mr,1)/120);   % half-height in trial rows, scaled so the mark reads at any trial count
        plot([spikeColsFull; spikeColsFull], ...                                             % x: one column per active bin, same bin index top and bottom
            [repmat(RasterTrials-halfH,1,numel(spikeColsFull)); repmat(RasterTrials+halfH,1,numel(spikeColsFull))], ...  % y: spans halfH rows either side of the selected trial
            'r-', 'LineWidth', 1);   % red vertical tick per displayed spike bin, drawn on top
    end




    %%%%%% Plot PSTH
    %%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subplot(18,1,[17 18])

    MRhist = BuildBurstMatrix(goodU(:,u),round(p.t),round((directimesSorted-preBase)),round((stimDur+preBase*2)));

    MRhist = squeeze(MRhist(trials,:,:));

    [nT2,nB2]= size(MRhist);

    spikeTimes = repmat([1:nB2],nT2,1);

    spikeTimes = spikeTimes(logical(MRhist));
    % Define bin edges (adjust for resolution)
    binWidth = 125; % 10 ms bins
    % Bin count must come from size(Mr,3), not nB: the loop above does
    % [nT,nN,nB]=size(Mr2), and when mergeTrials>1 Mr2 is 2D (trials x bins), so
    % that call leaves nB=1 and this condition never fired — the PSTH silently
    % used 125 ms bins for long epochs whenever trials were merged, making the
    % bin width depend on MergeNtrials rather than on the epoch length.
    if size(Mr,3)>300
         binWidth = 250;
    end
    edges = [1:binWidth:round((stimDur+preBase*2))]; % Adjust time window as needed

    % Compute histogram
    psthCounts = histcounts(spikeTimes, edges);

    % Convert to firing rate (normalize by bin width)
    psthRate = (psthCounts / (binWidth * nT2))*1000;

    b=bar(edges(1:end-1), psthRate,'histc');
    b.FaceColor = 'k';
    b.FaceAlpha = 0.3;
    b.MarkerEdgeColor = "none";
    xlim([0 round((stimDur+preBase*2)/100)*100])

    try %zero spiking in selection
        ylim([0 max(psthRate)+std(psthRate)])
    catch
        ur = ur+1;
        close
        continue
    end

    xticks([0 preBase:600:(stimDur+preBase*2) round((stimDur+preBase*2)/100)*100])

    xline([preBase stimDur+preBase],'LineWidth',1.5)

    xticklabels([-(preBase) 0:600:round((stimDur/100))*100 round((stimDur/100))*100 + 2*preBase]./1000)

    ax = gca; % Get current axes
    ax.XAxis.FontSize = 8;
    ax.XAxis.FontName = 'helvetica';

    ax.YAxis.FontSize = 8;
    ax.YAxis.FontName = 'helvetica';

    ylabel('[spk/s]','FontSize',10,'FontName','helvetica');
    xlabel('Time [s]','FontSize',10,'FontName','helvetica');

    ylims = ylim;
    yticks([round(ylims(2)/10)*5 ceil(ylims(2)/10)*10])


    %%%%PLot raw data several trials one
    %%%%channel%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Mark selected trial (RasterTrials was already picked above by the real-spike search)
    RasterTrials = trials(best_row);

    chan = goodU(1,u);

    subplot(18,1,[1 3])

    startTimes = directimesSorted(RasterTrials)+start-preBase;

    % A negative startTimes means the requested raw-data window would begin before
    % the recording/trial could exist — a directly-detectable error condition. Warn
    % explicitly and skip this neuron's raw-data panel rather than passing a negative
    % time into NP.getData/BuildBurstMatrix.
    if startTimes < 0
        warning('plotRaster:NegativeStartTime: neuron phyID=%d (unit idx %d), experiment %s: computed startTimes=%.2f ms is negative (RasterTrials=%d, start=%.2f ms, preBase=%d ms). Skipping raw-data extraction for this neuron.', ...
            phy_IDg(u), u, obj.dataObj.recordingName, startTimes, RasterTrials, start, preBase);
        close
        ur = ur+1;
        continue
    end

    freq = "AP"; %or "LFP"

    typeData = "line"; %or heatmap

    spikes = squeeze(BuildBurstMatrix(goodU(:,u),round(p.t),round(startTimes),round((window))));


    [fig, mx, mn] = PlotRawDataNP(obj,fig = fig,chan = chan, startTimes = startTimes,...
        window = window,spikeTimes = spikes);
    
    ax = gca; 
    ax.YAxis.FontSize = 8;
    ax.YAxis.FontName = 'helvetica';
   
    xlims= xlim;
    xticks([0:(xlims(2)/5):xlims(2)])
    xticklabels([0:100:window])
    ax.XAxis.FontSize = 8;
    ax.XAxis.FontName = 'helvetica';

    xlabel(string(chan))
    xline(-start/1000,'LineWidth',1.5)
    xlabel('Time [ms]','FontName','helvetica','FontSize',10) 

    ax.XRuler.TickDirection = 'out';   % ticks only outward (bottom)
    ax.XAxisLocation = 'bottom';
    ylabel('[\muV]','FontSize',10,'FontName','helvetica')
    title({sprintf('U.%d-Chan-%d-U.phy-%d-p=%.4f',u,channels(ur),phy_IDg(u),pvals(2,ur)),...
        sprintf('Ball-sizes-deg-%s',sizesString)});

    %%%%%%%%%%% Plot raster of selected trials
    %%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % subplot(20,1,[5 7])
    % 
    % RasterTrials =  trials;
    % 
    % bin3 = 4;
    % trialM = BuildBurstMatrix(goodU(:,u),round(p.t/bin3),round((directimesSorted+start)/bin3),round((window)/bin3));
    % 
    % if numel(RasterTrials) == 1
    %     TrialM = squeeze(trialM(RasterTrials,:,:))';
    % else
    %     TrialM = squeeze(trialM(trials,:,:));
    % end
    % 
    % %TrialM(TrialM~=0) = 0.3;
    % spikes1 = TrialM(TrialNumber,:);
    % spikeLoc = find(spikes1 >0);
    % if isempty(spikeLoc)
    %     close
    %     ur = ur+1;
    %     continue
    % end
    % %TrialM(TrialNumber,spikeLoc) = 1;
    % 
    % %Select offset in which selected trial belongs (10 trials)
    % imagesc(TrialM);colormap(flipud(gray(64)));
    % caxis([0 1])
    % xline([preBase stimDur+preBase],'LineWidth',1.5)
    % ylabel([sprintf('%d trials',numel(trials))])
    % 
    % xticks([1:50/bin3:window/bin3+1])
    % xticklabels([0:50:window])
    % 
    % ax = gca;
    % ax.XRuler.TickDirection = 'out';   % ticks only outward (bottom)
    % ax.XAxisLocation = 'bottom';
    % 
    % xlabel('Milliseconds')
    % set(gca,'FontSize',7)
    % 
    % xline(-start/bin3,'LineWidth',1.5)
    % 
    % xline(spikeLoc,'LineWidth',1,'Color','r','Alpha',0.3)
    % 
    % yline(TrialNumber,'LineWidth',3,'Color','r','Alpha',0.3) %Mark trial

    %fig.Position = [1     1   366   379];

    set(fig, 'Units', 'centimeters');
    set(fig, 'Position', [20 20 9 12]);

    if params.PaperFig
        obj.printFig(fig,sprintf('%s-%s-MovBall-SelectedTrials-eNeuron-%d',obj.dataObj.recordingName,fieldName,u),PaperFig = params.PaperFig)
    elseif params.overwrite
        obj.printFig(fig,sprintf('%s-%s-MovBall-SelectedTrials-eNeuron-%d',obj.dataObj.recordingName,fieldName,u))
    end

    if ur ~= length(eNeuron)
        close
    end
    
    ur = ur+1;

end %end eNeuron for loop

end %end plotRaster