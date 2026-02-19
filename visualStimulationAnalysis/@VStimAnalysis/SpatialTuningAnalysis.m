%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate tuning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% indexes

function fig = SpatialTuningAnalysis(expList, Stims2Comp,params)

arguments
    expList  (1,:) double  %%Number of experiment from excel list
    Stims2Comp cell %% Comparison order {'MB','RG','MBR'} would select neurons responsive to moving ball and 
    % compare this neurons responses to other stimuli. 
    params.threshold = 0.05;
    params.overwrite = false;
    params.overwriteResponse = false;
    params.percentile = 90; %Tope percentile selected to calculate tuning
    params.nShuffle = 1000;
end

%1. Bin the receptive field into offsetxoffset grid and
%take the mean within each grid.


rowT = cell(1,numel(respU));
colT = cell(1,numel(respU));

SpatialTuningIndex = zeros(1,numel(respU));

figure;tiledlayout(9,10,"TileSpacing","tight")

MrC = ConvBurstMatrix(Mr,fspecial('gaussian',[1 10],3),'same');

shuffC = ConvBurstMatrix(squeeze(mean(shuffledData,4)),fspecial('gaussian',[1 10],3),'same');

for u =1:length(respU)

    % Example 54x54 matrix
    rfi = squeeze(RFuST(:,:,u));


    % Define block size
    blockSize = redCoorY/offsetN;

    % Reshape and compute the mean for each 9x9 block
    meanGrid = blockproc(rfi, [blockSize blockSize], @(x) mean(x.data(:)));

    % Find mean of top blocks and rest, and the location of
    % the top blocks
    TOPper = prctile(meanGrid(:), params.percentile);

    Rtop = mean(meanGrid(meanGrid >= TOPper));

    [rowT{u}, colT{u}] = find(meanGrid >= TOPper);

    Rrest = mean(meanGrid(meanGrid < TOPper));

    % Calculate first term

    realTerm = (Rtop - Rrest)/mean(meanGrid,'all');

    %Do the same for the shuffled

    shuffTerm =zeros(1,numel(respU));

    for s = 1:params.nShuffle 

        rfiShuff = squeeze(RFuShuffST(:,:,u,s));

        meanGridShuff = blockproc(rfiShuff, [blockSize blockSize], @(x) mean(x.data(:)));
        %
        %                         % Use the same positions of the top blocks in the
        %                         % shuffled data
        %
        %                         RtopS = mean(meanGridShuff(rowT{u},colT{u}),'all');
        %
        %                         Mask = true(size(meanGridShuff));
        %
        %                         Mask(rowT{u},colT{u}) = false;
        %
        %                         RrestS = mean(meanGridShuff(Mask));

        % Find mean of top blocks and rest, and the location of
        % the top blocks
        TOPper = prctile(meanGridShuff(:), params.percentile);

        RtopS = mean(meanGridShuff(meanGridShuff >= TOPper));

        RrestS = mean(meanGridShuff(meanGridShuff < TOPper));

        % Calculate first term

        shuffTerm(s) = (RtopS - RrestS)/mean(meanGridShuff,'all');

    end

    SpatialTuningIndex(u) = realTerm-median(shuffTerm);

    nexttile
    imagesc(rfi);title(string(SpatialTuningIndex(u)))
    axis off;

    %imagesc(squeeze(shuffC(:,u,:)));title(string(shuffTerm(u)))

    %                     imagesc(squeeze(MrC(:,respU(u),:)));title(string(SpatialTuningIndex(u)))
    %                     yline([10:10:size(Mr,1)],'w','LineWidth',1);
    %                     yline([90:90:size(Mr,1)],'w','LineWidth',3);

end
