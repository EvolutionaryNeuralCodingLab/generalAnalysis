function [colorbarLims] = PlotReceptiveFields(obj, params)
% PlotReceptiveFields  Plots the spatial receptive fields of neurons recorded
%                      during a moving-ball stimulus, optionally filtered by
%                      direction, size, and luminance.
%
% OUTPUT
%   colorbarLims  – [cMin cMax] limits of the last colorbar drawn.

% ── Input argument block ─────────────────────────────────────────────────────
arguments (Input)
    obj                                   % Parent analysis object
    params.overwrite           logical = false   % Overwrite existing figure files
    params.analysisTime                = datetime('now')  % Timestamp for provenance
    params.inputParams                 = false           % If true, print params and exit
    params.exNeurons                   = nan;            % Explicit neuron indices to plot
    params.AllSomaticNeurons           = false           % Plot every somatic unit
    params.AllResponsiveNeurons        = true            % Plot only statistically responsive units
    params.fixedWindow                 = false           % Use a fixed time window
    params.speed                       = 1;              % Stimulus speed index (1 = slow, 2 = fast)
    params.noEyeMoves                  = false           % Use no-eye-movement recording mode
    params.reduceFactor                = 20             % Spatial downsampling factor for the RF map
    params.allCombined                 = false           % Also plot the direction-summed RF
    params.eye_to_monitor_distance     = 21.5            % Eye-to-monitor distance in cm
    params.pixel_size                  = 33              % Physical size of one pixel in µm
    params.resolution                  = 1080            % Monitor vertical resolution in pixels
    params.meanAllNeurons              = false           % Average RF across all neurons before plotting
    params.PaperFig           logical  = false           % Apply paper-quality rendering settings
    params.OneDirection        string  = "all"           % Restrict plot to one direction ("up","left","down","right","all")
    params.OneLuminosity       string  = "all"           % Restrict plot to one luminance ("black","white","all")
    params.OneSize             string  = "all"           % Restrict plot to one size ("small","middle","big","all")
    params.colorbarLims                = []              % Optional manual colorbar limits [cMin cMax]
    params.tickNum                     = 3               % Number of ticks per axis   
end

% ── Debug helper: print parameter struct and exit ────────────────────────────
if params.inputParams, disp(params), return, end

% ── Load pre-computed statistics and receptive fields ───────────────────────
Stats = obj.ShufflingAnalysis;                              % Shuffling-based significance statistics
RFs   = obj.CalculateReceptiveFields('speed', params.speed); % Receptive-field maps at the chosen speed

% Build the field name that indexes speed-specific substructs (e.g., "Speed1")
fieldName = sprintf('Speed%d', params.speed);

% Extract per-unit response p-values for the chosen speed
pvals = Stats.(fieldName).pvalsResponse;

% Load the response window data to recover the stimulus condition table
responses = obj.ResponseWindow;

% Extract the unique values of each stimulus dimension from the condition matrix
% Columns of C are assumed to encode: [?, direction, ?, size, ?, luminance]
uDir  = unique(responses.(fieldName).C(:, 2));  % Unique motion directions (radians)
uSize = unique(responses.(fieldName).C(:, 4));  % Unique stimulus sizes
uLum  = unique(responses.(fieldName).C(:, 6));  % Unique luminance levels

% ── Resolve which direction(s) to include ────────────────────────────────────
if params.OneDirection ~= "all"
    switch params.OneDirection
        case "up"
            dirIDX = find(uDir == 0);       % 0 rad = upward motion
        case "left"
            dirIDX = find(uDir == 1.57);    % π/2 rad ≈ leftward motion
        case "down"
            dirIDX = find(uDir == 3.14);    % π rad ≈ downward motion
        case "right"
            % BUG FIX: was find(uDir==1.57) which is identical to "left".
            % Right corresponds to 3π/2 ≈ 4.71 rad.
            dirIDX = find(uDir == 4.71);
        otherwise
            error("Unknown OneDirection value: %s", params.OneDirection)
    end
    DirectionSelected = uDir(dirIDX);   % Scalar: the single selected direction value
else
    dirIDX         = 1:numel(uDir);     % All direction indices
    DirectionSelected = uDir;           % All direction values
end

% ── Resolve which luminance(s) to include ────────────────────────────────────
if params.OneLuminosity ~= "all"
    switch params.OneLuminosity
        case "black"
            lumIDX = find(uLum == 1);    % Luminance value 1 = black stimulus
        case "white"
            lumIDX = find(uLum == 255);  % Luminance value 255 = white stimulus
        otherwise
            error("Unknown OneLuminosity value: %s", params.OneLuminosity)
    end
    LuminositySelected = uLum(lumIDX);  % Scalar: the single selected luminance value
else
    lumIDX            = 1:numel(uLum);  % All luminance indices
    LuminositySelected = uLum;          % All luminance values
end

% ── Resolve which size(s) to include ─────────────────────────────────────────
if params.OneSize ~= "all"
    switch params.OneSize
        case "small"
            sizeIDX = 1;                % First (smallest) size
        case "middle"
            sizeIDX = 2;                % Middle size
        case "big"
            sizeIDX = 3;                % Last (largest) size
        otherwise
            % BUG FIX: error message incorrectly referenced params.OneLuminosity
            error("Unknown OneSize value: %s", params.OneSize)
    end
    % BUG FIX: variable was named SizeSelected (no 's') but referenced later
    % as SizesSelected, causing "Undefined variable" at runtime.
    SizesSelected = uSize(sizeIDX);     % Scalar: the single selected size value
else
    sizeIDX       = 1:numel(uSize);     % All size indices
    SizesSelected = uSize;              % All size values
end

% ── Select which neurons to process ──────────────────────────────────────────
if isnan(params.exNeurons)
    if params.AllSomaticNeurons
        % Use every unit regardless of responsiveness
        eNeuron = 1:numel(pvals);
        pvals   = [eNeuron; pvals(eNeuron)]; % Row 1: indices, Row 2: p-values
    elseif params.AllResponsiveNeurons
        % Keep only units whose response p-value is below α = 0.05
        eNeuron = find(pvals < 0.05);
        pvals   = [eNeuron; pvals(eNeuron)];
        if isempty(eNeuron)
            fprintf('No responsive neurons.\n')
            return
        end
    end
else
    % Use the explicitly provided unit indices
    eNeuron = params.exNeurons;
    pvals   = [eNeuron; pvals(eNeuron)];
end

% ── Build the reduced-resolution coordinate grid ─────────────────────────────
coorRect    = obj.VST.rect';                                    % Screen rectangle [x y w h] (pixels), transposed to column
% SUGGESTION: consider using obj.VST.rect directly with named fields to
% avoid silent breakage if the rectangle layout changes.

% Downsampling factor must not exceed the smallest ball radius
reduceFactor = min([params.reduceFactor, min(obj.VST.ballSizes)]);

% Reduced-resolution grid dimensions
redCoorX = round(coorRect(3) / reduceFactor);   % Width in reduced pixels
redCoorY = round(coorRect(4) / reduceFactor);   % Height in reduced pixels

% Convert pixel size to cm at the reduced resolution
pixel_size_cm = params.pixel_size / (params.resolution / reduceFactor);

% Build the visual-angle (degrees) coordinate arrays for the reduced grid
monitor_resolution = [redCoorX, redCoorY];                       % [width height] in reduced pixels
[theta_x, theta_y] = pixels2eyeDegrees(params.eye_to_monitor_distance, ...
                                        pixel_size_cm, monitor_resolution);

% Crop theta_x horizontally so it is square (same spatial extent as theta_y)
% The crop removes the equal-width margins on left and right
theta_x = theta_x(:, 1 + (redCoorX - redCoorY)/2 : (redCoorX - redCoorY)/2 + redCoorY);

% ── Pre-compute symmetric, degree-based tick positions (shared across tiles) ─
% SUGGESTION: Computing ticks once here and reusing inside the loop avoids
% repeated interp1 calls and guarantees all subplots share identical ticks.

% X-axis: find the largest absolute visual angle present in the cropped grid
maxDeg_x    = max(abs(theta_x(1, :)));
% Define 5 tick positions in degrees, symmetric about 0, stopping 1° inside the edge
tickDeg_x   = linspace(-(maxDeg_x - 5), maxDeg_x - 5, params.tickNum);
% Map degree values back to reduced-pixel column indices via linear interpolation
tickPix_x   = interp1(theta_x(1, :), 1:size(theta_x, 2), tickDeg_x, 'linear', 'extrap');

% Y-axis: same procedure on the first column of theta_y
maxDeg_y    = max(abs(theta_y(:, 1)));
tickDeg_y   = linspace(-(maxDeg_y - 5), maxDeg_y - 5,  params.tickNum);
tickPix_y   = interp1(theta_y(:, 1), 1:size(theta_y, 1), tickDeg_y, 'linear', 'extrap');

% Pre-round degree labels so they are computed only once
tickLbl_x   = round(tickDeg_x);   % Rounded x-axis degree labels for display
tickLbl_y   = round(tickDeg_y);   % Rounded y-axis degree labels for display

% ── Load the appropriate RF array ────────────────────────────────────────────
if params.noEyeMoves
    % Load the no-eye-movement RF (previously saved per-recording)
    RFu = squeeze(load(sprintf('NEM-RFuST-Q1-Div-X-%s', NP.recordingName)).RFuST);
else
    RFu            = RFs.RFuST;              % Direction-summed RF map [y × x × unit]
    RFuDirSizeLum  = RFs.RFuDirSizeLumFilt; % RF array split by [dir × size × lum × y × x × unit]
end

% Build a 2-D Gaussian smoothing kernel scaled to the RF map size
% Kernel size is proportional to the number of spatial offsets used
offsetN      = numel(unique(obj.VST.offsets));                           % Number of unique stimulus offsets
TwoDGaussian = fspecial('gaussian', floor(size(RFu, 2) / (offsetN / 2)), ... % Kernel window (px)
                         redCoorY / offsetN);                            % Kernel sigma (px)
% SUGGESTION: fspecial is from the Image Processing Toolbox.
% imgaussfilt or a manually constructed kernel can be used as a fallback.

% ═══════════════════════════════════════════════════════════════════════════════
% Main loop: one iteration per selected neuron
% ═══════════════════════════════════════════════════════════════════════════════
for u = eNeuron

    %ru = find(eNeuron == u);   % Index of neuron u within the eNeuron vector
    ru =u;
    % ── Optional: plot the direction-summed (combined) RF ──────────────────
    if params.allCombined

        figRF = figure;   % Open a new figure window
        % Display the Gaussian-smoothed combined RF as a colour image
        imagesc(squeeze(conv2(RFu(:, :, ru), TwoDGaussian, 'same')));

        c = colorbar;                   % Attach a colourbar
        title(c, 'spk/s')              % Label colourbar units

        colormap('turbo')              % Apply perceptually-uniform colour map
        title(sprintf('u-%d', u))      % Title with unit number

        % Apply symmetric x ticks (degrees)
        xticks(tickPix_x);
        xticklabels(tickLbl_x);

        % Apply symmetric y ticks (degrees)
        yticks(tickPix_y);
        yticklabels(tickLbl_y);

        axis image    % Equal aspect ratio, no white space

        figRF.Position = [680 577 156 139];   % Set figure size (pixels)

        % Save figure to disk if overwrite flag is set
        if params.noEyeMoves
            if params.overwrite
                if params.PaperFig
                    obj.printFig(figRF, sprintf('%s-NEM-MovBall-ReceptiveField-eNeuron-%d.pdf', ...
                                  obj.dataObj.recordingName, u), PaperFig = params.PaperFig);
                else
                    obj.printFig(figRF, sprintf('%s-NEM-MovBall-ReceptiveField-eNeuron-%d.pdf', ...
                                  obj.dataObj.recordingName, u));
                end
            end
        else
            if params.overwrite
                if params.PaperFig
                    obj.printFig(figRF, sprintf('%s-MovBall-ReceptiveField-eNeuron-%d', ...
                                  obj.dataObj.recordingName, u), PaperFig = params.PaperFig);
                else
                    obj.printFig(figRF, sprintf('%s-MovBall-ReceptiveField-eNeuron-%d', ...
                                  obj.dataObj.recordingName, u));
                end
            end
        end

    end % allCombined

    % ── Extract this neuron's RF slice and apply any dimension filters ──────

    if params.meanAllNeurons
        % Average the RF across all neurons (dimension 6) and collapse it
        RFuRed = reshape(mean(RFuDirSizeLum, 6), ...
                         [size(RFuDirSizeLum, 1), size(RFuDirSizeLum, 2), ...
                          size(RFuDirSizeLum, 3), size(RFuDirSizeLum, 4), ...
                          size(RFuDirSizeLum, 5)]);
        % BUG: hasNotString is never defined; the intent seems to be to
        % additionally average over the non-compared dimensions (e.g., if
        % only direction is compared, average over size and luminance).
        % Define hasNotString before this block, e.g.:
        %   hasNotString = [];
        %   if params.OneSize ~= "all",      hasNotString(end+1) = 2; end
        %   if params.OneLuminosity ~= "all", hasNotString(end+1) = 3; end
        %   if params.OneDirection ~= "all",  hasNotString(end+1) = 1; end
        % Then the loop below is correct:
        for i = 1:numel(hasNotString)   % Average over dimensions not being compared
            RFuRed = mean(RFuRed, hasNotString(i));
        end
    else
        % Extract the RF for neuron ru; keep all 5 condition dimensions explicit
        RFuRed = reshape(RFuDirSizeLum(:, :, :, :, :, ru), ...
                         [size(RFuDirSizeLum, 1), size(RFuDirSizeLum, 2), ...
                          size(RFuDirSizeLum, 3), size(RFuDirSizeLum, 4), ...
                          size(RFuDirSizeLum, 5)]);   % [dir × size × lum × y × x]

        % Apply size filter if requested (select single size slice)
        if params.OneSize ~= "all"
            RFuRed = RFuRed(:, sizeIDX, :, :, :);
        end

        % Apply luminance filter if requested
        if params.OneLuminosity ~= "all"
            RFuRed = RFuRed(:, :, lumIDX, :, :);
        end

        % Apply direction filter if requested
        if params.OneDirection ~= "all"
            RFuRed = RFuRed(dirIDX, :, :, :, :);
        end
    end

    % ── Determine colour-axis limits from this neuron's data ────────────────
    cMax = max(RFuRed, [], 'all');   % Global maximum firing rate across all conditions
    cMin = min(RFuRed, [], 'all');   % Global minimum firing rate across all conditions

    % ── Compute tile-grid layout for the tiled figure ───────────────────────
    % BUG FIX: was prod(size(RFuRed,[1 2 3])) which always produces a scalar,
    % making the numel==3 branch unreachable. Corrected to size().
    tilesSize = size(RFuRed, [1 2 3]);   % [nDir, nSize, nLum] – counts of condition tiles

    % Flatten to a [rows × cols] layout: rows = directions, cols = size × lum
    tilesSize = [tilesSize(1), tilesSize(2) * tilesSize(3)];

    % ── Create the tiled figure ──────────────────────────────────────────────
    figRF = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Full-screen window
    NeuronLayout = tiledlayout(tilesSize(1), tilesSize(2), ...
                               "TileSpacing", "tight", "Padding", "compact");

    j = 0;   % Running tile counter (used to attach colorbar to the last tile)

    % ── Inner loops: one tile per (direction × size × luminance) combination ─
    for d = 1:size(RFuRed, 1)     % Iterate over direction slices
        for s = 1:size(RFuRed, 2) % Iterate over size slices
            for l = 1:size(RFuRed, 3) % Iterate over luminance slices

                ax = nexttile;   % Advance to the next tile in the layout

                % Display the RF heat map for condition (d, s, l)
                imagesc(squeeze(RFuRed(d, s, l, :, :)));

                % Lock colour axis to the per-neuron global range for comparability
                % SUGGESTION: clim() is the modern replacement for the deprecated caxis()
                clim([cMin, cMax]);

                % Style y-axis tick labels
                axi = gca;
                axi.YAxis.FontSize = 8;
                axi.YAxis.FontName = 'helvetica';

                % Axis labels (degrees of visual angle)
                xlabel('Degrees', 'FontSize', 10, 'FontName', 'helvetica')
                ylabel('Degrees', 'FontSize', 10, 'FontName', 'helvetica')

                % Style x-axis tick labels
                axi.XAxis.FontSize = 8;
                axi.XAxis.FontName = 'helvetica';

                colormap('turbo')   % Perceptually-uniform colour map

                % BUG FIX: was uDir(d), uSize(s), uLum(l) – these index the
                % FULL unfiltered arrays, so the label is wrong whenever a
                % filter is active.  Use the Selected vectors instead.
                title(sprintf('Dir-%.2f-Size-%.0f-Lum-%.0f', ...
                              DirectionSelected(d), SizesSelected(s), LuminositySelected(l)), ...
                      'FontSize', 4)

                % ── Apply symmetric degree-based ticks ──────────────────────
                xticks(tickPix_x);         % Set x tick positions (reduced-pixel units)
                xticklabels(tickLbl_x);    % Label with rounded degree values
                yticks(tickPix_y);         % Set y tick positions (reduced-pixel units)
                yticklabels(tickLbl_y);    % Label with rounded degree values

                j = j + 1;   % Increment tile counter

                % Attach colourbar only to the final tile so it doesn't clutter the layout
                if j == size(RFuRed, 1) * size(RFuRed, 2) * size(RFuRed, 3)
                    c = colorbar;
                    title(c, 'spk/s', 'FontSize', 8, 'FontName', 'helvetica')

                    % Override colour limits if the caller supplied explicit bounds
                    if ~isempty(params.colorbarLims)
                        clim(params.colorbarLims);
                    end

                    colorbarLims = c.Limits;   % Return the final colourbar limits
                end

                axis(ax, 'image');         % Equal aspect ratio so the RF map is not distorted
                % SUGGESTION: axis equal already enforces equal scaling;
                % pbaspect([1 1 1]) is therefore redundant here and can be removed.

            end % luminance
        end % size
    end % direction

    % ── Build filename strings from the selected condition labels ────────────
    Sdir  = strjoin(string(DirectionSelected),  "-");   % e.g. "0-1.57-3.14-4.71"
    Ssize = strjoin(string(SizesSelected),      "-");   % e.g. "5-10-20"
    Slum  = strjoin(string(LuminositySelected), "-");   % e.g. "1-255"

    % ── Handle the mean-all-neurons special case ─────────────────────────────
    if params.meanAllNeurons
        title(NeuronLayout, 'MeanAllUnits');   % Label the whole layout
        if params.overwrite
            obj.printFig(figRF, sprintf('%s-%s-MovBall-RF-sep-%s-Mean', ...
                          obj.dataObj.recordingName, fieldName, ...
                          sprintf('Dir-%s-Size-%s-Lum-%s', Sdir, Ssize, Slum)));
        end
        return   % Only one figure is produced; no per-neuron loop needed
    end

    % ── Resize figure for single-row layouts ────────────────────────────────
    if tilesSize(1) == 1
        set(figRF, 'Units', 'centimeters');
        set(figRF, 'Position', [2 2 4 4]);   % Compact format for one-row grids
    end

    % ── Save figure ──────────────────────────────────────────────────────────
    if ~params.noEyeMoves   % Standard (eye-movement) recording mode
        if params.overwrite
            if params.PaperFig
                obj.printFig(figRF, sprintf('%s-%s-MovBall-RF-sep-%s-eNeuron-%d', ...
                              obj.dataObj.recordingName, fieldName, ...
                              sprintf('Dir-%s-Size-%s-Lum-%s', Sdir, Ssize, Slum), u), ...
                              PaperFig = params.PaperFig);
            else
                obj.printFig(figRF, sprintf('%s-%s-MovBall-RF-sep-%s-eNeuron-%d', ...
                              obj.dataObj.recordingName, fieldName, ...
                              sprintf('Dir-%s-Size-%s-Lum-%s', Sdir, Ssize, Slum), u));
            end
        end
    end

    % Close the figure unless this is the last neuron (keep last open for inspection)
    if u ~= eNeuron(end)
        close
    end

end % neuron loop

end % function