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
% -------------------------------------------------------------------------
% Reconstruct trajectories from stimulus geometry rather than raw data
% Raw obj.VST.ballTrajectoriesX/Y has sampling glitches — using min/max
% offsets from obj.VST.parallelsOffset and screen centre gives clean
% constant-velocity trajectories independent of sampling artefacts.
% -------------------------------------------------------------------------

% Frame count per (offset, direction) for this speed condition
nFramesFull = obj.VST.nFrames;  % [nSpeeds × nOffsets × nDirections]
if ndims(nFramesFull) == 3
    nFramesPerOffsetDir = squeeze(nFramesFull(speedIndx, :, :));  % [nOffsets × nDirs]
else
    nFramesPerOffsetDir = nFramesFull;
end

useOriginalCorrs = true;

if useOriginalCorrs

    Xpos = obj.VST.ballTrajectoriesX;
    Ypos = obj.VST.ballTrajectoriesY;

    if size(Xpos,1) > 1
        Xpos = Xpos(speedIndx,:,:,1:unique(obj.VST.nFrames(speedIndx,:,:)));
        Ypos = Ypos(speedIndx,:,:,1:unique(obj.VST.nFrames(speedIndx,:,:)));
    end

else

    % Reconstruct from stimulus design parameters
    [Xpos, Ypos] = reconstructBallTrajectoriesFromGeometry(obj, speedIndx, nFramesPerOffsetDir);

end
% Xpos, Ypos are [1 × nOffsets × nDirs × nFramesMax]

% Flatten trajectories to [nTrials × nFrames] matching trial ordering
%  in C
sizeX         = size(Xpos);          % [nSpeeds × nOffsets × nDirs × nFrames]
nSizes        = length(unique(obj.VST.ballSizes));
C             = obj.ResponseWindow.(sprintf('Speed%d', speedIndx)).C;  % category matrix
trialDivVid   = size(C,1) / numel(unique(C(:,2))) / numel(unique(C(:,3))) ...
    / numel(unique(C(:,4))) / numel(unique(C(:,5)));  % trials per unique video

nFrames  = sizeX(end);
nTrials  = size(C,1);

% Build [nTrials × nFrames] position arrays matching C order
% Loop structure MUST match the order used to build C (dir × offset ×
% speed)member you can have more than one speed
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
dwellFrames   = nan(nTrials, nCells);     % initialise dwell time as 0

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

        cx = cropOffsetX + (gx - 0.5) * cellW;   % CORRECT
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

%figure;imagesc(reshape(mean(dwellFrames,1),[9,9]))
% 
% figure;imagesc(reshape(mean(dwellFrames, 1, 'omitnan'), [nGrid nGrid]));
% colorbar; title('Original trajectories');

% test = Xpos(1,:,1,:);
% figure;hist(test(:))

% -------------------------------------------------------------------------
% Identify valid grid cells per direction
% Cell is valid for direction d if at least one trial in that direction crosses it
% -------------------------------------------------------------------------
directions   = C(:,2);              % direction label per trial
uDirs        = unique(directions);  % unique direction values
nDirs        = numel(uDirs);
nOffsets     = numel(obj.VST.parallelsOffset);
centerX     = obj.VST.centerX;
centerY     = obj.VST.centerY;

validGridPerDirection = false(nCells, nDirs);
nTrialsPerCellDir     = zeros(nCells, nDirs);

for d = 1:nDirs
    trialsThisDir = directions == uDirs(d);
    % Cell is valid if any trial in this direction has a non-NaN crossing
    validGridPerDirection(:,d) = any(~isnan(crossingFrame(trialsThisDir,:)), 1)';
    % Trial count per cell for this direction
    nTrialsPerCellDir(:,d)     = sum(~isnan(crossingFrame(trialsThisDir,:)), 1)';
end

% figure;
% for d = 1:nDirs
%     subplot(1, nDirs, d);
%     hold on;
%     for o = 1:nOffsets
%         x = squeeze(Xpos(1, o, d, ~isnan(Xpos(1,o,d,:))));
%         y = squeeze(Ypos(1, o, d, ~isnan(Ypos(1,o,d,:))));
%         plot(x, y, 'b-');
%     end
%     plot(centerX, centerY, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
%     rectangle('Position', [0 0 screenW screenH], 'EdgeColor', 'k', 'LineWidth', 1.5);
%     title(sprintf('Direction %d', d));
%     axis equal;
%     xlim([-screenW screenW*2]);
%     ylim([-screenH screenH*2]);
% end
end

function [dx_dir, dy_dir] = inferDirectionVector(XposRaw, YposRaw, dirIdx)
% Infer motion direction by comparing trajectory endpoints
% Returns unit vector (dx_dir, dy_dir) for the motion direction

% Average across offsets to get robust endpoints (immune to glitches in single trajectories)
xStart = mean(XposRaw(:, dirIdx, 1),  'omitnan');     % first frame across offsets
xEnd   = mean(XposRaw(:, dirIdx, end), 'omitnan');    % last frame across offsets
yStart = mean(YposRaw(:, dirIdx, 1),  'omitnan');
yEnd   = mean(YposRaw(:, dirIdx, end), 'omitnan');

% Motion vector
dx = xEnd - xStart;
dy = yEnd - yStart;

% Normalise to unit vector
mag    = sqrt(dx^2 + dy^2);
dx_dir = dx / mag;
dy_dir = dy / mag;
end

function [Xpos, Ypos] = reconstructBallTrajectoriesFromGeometry(obj, speedIndx, nFramesPerOffsetDir)
% reconstructBallTrajectoriesFromGeometry - Reconstruct ball trajectories from
% stimulus design parameters rather than (potentially glitched) raw trajectory data.
%
% Direction of motion is inferred from raw trajectory endpoints (first vs last
% frame across offsets), avoiding any assumption about angle conventions.
% Offset is applied perpendicular to motion direction.

% Stimulus geometry
centerX     = obj.VST.centerX;
centerY     = obj.VST.centerY;
screenW     = obj.VST.rect(3);
screenH     = obj.VST.rect(4);
offsets     = obj.VST.parallelsOffset;
directions  = unique(obj.VST.directions);

nOffsets   = numel(offsets);
nDirs      = numel(directions);
nFramesMax = max(nFramesPerOffsetDir(:));

Xpos = nan(1, nOffsets, nDirs, nFramesMax);
Ypos = nan(1, nOffsets, nDirs, nFramesMax);

travelDist = sqrt(screenW^2 + screenH^2);

% Load raw trajectories for direction inference
% Squeeze speed dim so we have [nOffsets × nDirs × nFrames]
XposRawFull = obj.VST.ballTrajectoriesX;
YposRawFull = obj.VST.ballTrajectoriesY;
if size(XposRawFull,1) > 1
    XposRaw = squeeze(XposRawFull(speedIndx, :, :, :));
    YposRaw = squeeze(YposRawFull(speedIndx, :, :, :));
else
    XposRaw = squeeze(XposRawFull);
    YposRaw = squeeze(YposRawFull);
end

for d = 1:nDirs
    % Infer direction vector from raw data — robust to angle convention
    [dx_dir, dy_dir] = inferDirectionVector(XposRaw, YposRaw, d);

    % Perpendicular vector (rotate 90° counterclockwise in screen coordinates)
    dx_perp = -dy_dir;
    dy_perp =  dx_dir;

    for o = 1:nOffsets
        offsetVal = offsets(o);
        nFr       = nFramesPerOffsetDir(o, d);

        % Trajectory midpoint = screen centre + offset perpendicular to motion
        midX = centerX + offsetVal * dx_perp;
        midY = centerY + offsetVal * dy_perp;

        % Start/end points along motion direction
        xStart = midX - (travelDist/2) * dx_dir;
        yStart = midY - (travelDist/2) * dy_dir;
        xEnd   = midX + (travelDist/2) * dx_dir;
        yEnd   = midY + (travelDist/2) * dy_dir;

        % Linear interpolation across frames — constant velocity
        Xpos(1, o, d, 1:nFr) = linspace(xStart, xEnd, nFr);
        Ypos(1, o, d, 1:nFr) = linspace(yStart, yEnd, nFr);
    end
end
end