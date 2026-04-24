function [crossingFrame, dwellFrames, validGridPerDirection, nTrialsPerCellDir] = ...
    computeBallGridCrossings(obj, speedIndx, params)
% computeBallGridCrossings - Detect when ball centre crosses each grid cell.
%
% For each trial and each grid cell, finds the frame at which the ball centre
% first enters the spatial bin corresponding to that grid cell. Grid cells are
% defined on the original screen coordinates (not the reduced coordinates used
% elsewhere for receptive field plotting).
%
% Inputs:
%   obj       - experiment object with VST metadata and stimulus category matrix C
%   speedIndx - speed index into VST trajectory arrays (1 or 2)
%   params    - parameter struct containing GridSize (e.g. 9)
%
% Outputs:
%   crossingFrame         : [nTrials × nGridCells] frame at which ball centre
%                           first enters grid cell. NaN if ball never enters.
%   dwellFrames           : [nTrials × nGridCells] number of frames ball centre
%                           remains in cell. 0 if never entered.
%   validGridPerDirection : [nGridCells × nDirections] logical. True if at least
%                           one trial in that direction crosses the cell.
%   nTrialsPerCellDir     : [nGridCells × nDirections] trial count per cell×dir.

    % -------------------------------------------------------------------------
    % Load ball trajectory data
    % ChangePosX, ChangePosY: [nTrials × nFrames] ball centre pixel coordinates
    % -------------------------------------------------------------------------
    Xpos = obj.VST.ballTrajectoriesX;  % original 4D: [nSpeeds × nOffsets × nDirs × nFrames]
    Ypos = obj.VST.ballTrajectoriesY;


    if size(Xpos,1) > 1
        Xpos = Xpos(speedIndx,:,:,1:unique(obj.VST.nFrames(speedIndx,:,:)));  % select trajectories for this speed
        Ypos = Ypos(speedIndx,:,:,1:unique(obj.VST.nFrames(speedIndx,:,:)));
    end

    % Flatten trajectories to [nTrials × nFrames] matching trial ordering in C
    sizeX         = size(Xpos);          % [nSpeeds × nOffsets × nDirs × nFrames]
    nSizes        = length(unique(obj.VST.ballSizes));
    C             = obj.ResponseWindow.(sprintf('Speed%d', speedIndx)).C;  % category matrix
    trialDivVid   = size(C,1) / numel(unique(C(:,2))) / numel(unique(C(:,3))) ...
                  / numel(unique(C(:,4))) / numel(unique(C(:,5)));  % trials per unique video

    nFrames  = sizeX(end);
    nTrials  = size(C,1);

    % Build [nTrials × nFrames] position arrays matching C order
    % Loop structure MUST match the order used to build C (dir × offset × speed)
    ChangePosX = zeros(nTrials, nFrames);
    ChangePosY = zeros(nTrials, nFrames);
    j = 1;
    for d  = 1:sizeX(3)           % directions
        for of = 1:sizeX(2)       % offsets
            for sp = 1:sizeX(1)   % speeds (one after speedIndx selection)
                % Replicate trajectory across size categories and trial divisions
                traj = squeeze(Xpos(sp, of, d, :))';   % [1 × nFrames]
                ChangePosX(j:j+nSizes*trialDivVid-1, :) = repmat(traj, nSizes*trialDivVid, 1);
                trajY = squeeze(Ypos(sp, of, d, :))';
                ChangePosY(j:j+nSizes*trialDivVid-1, :) = repmat(trajY, nSizes*trialDivVid, 1);
                j = j + nSizes * trialDivVid;
            end
        end
    end
changePosDir1x = Xpos(:,:,1,:);
changePosDir1y = Ypos(:,:,1,:);
    % -------------------------------------------------------------------------
    % Define grid cells on original screen coordinates
    % Screen is [obj.VST.rect(3) × obj.VST.rect(4)] pixels
    % Each grid cell is cellW × cellH pixels
    % Cell (gx, gy) spans: x ∈ [(gx-1)*cellW, gx*cellW], y ∈ [(gy-1)*cellH, gy*cellH]
    % -------------------------------------------------------------------------
    screenW = obj.VST.rect(3);                  % full screen width in pixels
    screenH = obj.VST.rect(4);                  % full screen height in pixels
    nGrid   = params.GridSize;                  % e.g. 9 → 9×9 = 81 cells
    cellW   = screenH / nGrid;                  % width of each cell in pixels
    cellH   = screenH / nGrid;                  % height of each cell in pixels
    nCells  = nGrid * nGrid;                    % total number of grid cells

    cropOffsetX = (screenW - screenH)/2;


    % -------------------------------------------------------------------------
    % For each trial and cell: find first frame of entry and dwell duration
    % -------------------------------------------------------------------------
    crossingFrame = nan(nTrials, nCells);       % initialise as NaN (never entered)
    dwellFrames   = zeros(nTrials, nCells);     % initialise dwell time as 0

    for t = 1:nTrials
        % For each frame, determine which grid cell the ball centre is in
        gxPerFrame = floor((ChangePosX(t,:)  - cropOffsetX )/ cellW) + 1;  % [1 × nFrames] grid x index
        gyPerFrame = floor(ChangePosY(t,:) / cellH) + 1;  % [1 × nFrames] grid y index

        % % Clamp to valid grid range (ball may be off-screen during entry/exit)
        % gxPerFrame = max(1, min(nGrid, gxPerFrame));
        % gyPerFrame = max(1, min(nGrid, gyPerFrame));
        %
        % % Flatten (gx, gy) to linear cell index: cellIdx = (gy-1)*nGrid + gx
        % cellIdxPerFrame = (gyPerFrame - 1) * nGrid + gxPerFrame;  % [1 × nFrames]

        valid = gxPerFrame >= 1 & gxPerFrame <= nGrid & ...
            gyPerFrame >= 1 & gyPerFrame <= nGrid;

        cellIdxPerFrame = nan(size(gxPerFrame));
        cellIdxPerFrame(valid) = (gyPerFrame(valid)-1)*nGrid + gxPerFrame(valid);

        visitedCells = unique(cellIdxPerFrame(~isnan(cellIdxPerFrame)));

        for c = visitedCells
            inCell = cellIdxPerFrame == c;
            frames = find(inCell);

            if isempty(frames)
                continue
            end

            % --- cell center ---
            gx = mod(c-1, nGrid) + 1;
            gy = floor((c-1)/nGrid) + 1;

            cx = (gx - 0.5) * cellW;
            cy = (gy - 0.5) * cellH;

            % --- distance to center ---
            dx = ChangePosX(t, frames) - cx;
            dy = ChangePosY(t, frames) - cy;
            dist = sqrt(dx.^2 + dy.^2);

            % --- center crossing ---
            [~, minIdx] = min(dist);
            centerFrame = frames(minIdx);

            % --- exit ---
            exitFrame = frames(end);

            crossingFrame(t, c) = centerFrame;
            %dwellFrames(t, c)   = exitFrame - centerFrame + 1;
            dwellFrames(t,c) = numel(frames);
        end
    end

    figure;imagesc(reshape(mean(dwellFrames,1),[9,9]))
    % -------------------------------------------------------------------------
    % Identify valid grid cells per direction
    % Cell is valid for direction d if at least one trial in that direction crosses it
    % -------------------------------------------------------------------------
    directions   = C(:,2);              % direction label per trial
    uDirs        = unique(directions);  % unique direction values
    nDirs        = numel(uDirs);

    validGridPerDirection = false(nCells, nDirs);
    nTrialsPerCellDir     = zeros(nCells, nDirs);

    for d = 1:nDirs
        trialsThisDir = directions == uDirs(d);
        % Cell is valid if any trial in this direction has a non-NaN crossing
        validGridPerDirection(:,d) = any(~isnan(crossingFrame(trialsThisDir,:)), 1)';
        % Trial count per cell for this direction
        nTrialsPerCellDir(:,d)     = sum(~isnan(crossingFrame(trialsThisDir,:)), 1)';
    end
end