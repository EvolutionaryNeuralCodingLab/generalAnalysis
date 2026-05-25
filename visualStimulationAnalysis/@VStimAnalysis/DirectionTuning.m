function results = DirectionTuning(vsObj, params)
% DirectionTuning  Compute orientation & direction selectivity indices for
%   a single experiment from the ResponseWindow of any stimulus object that
%   contains a 'direction' category in its condition matrix.
%
%   OSI is computed via the circular-variance method (1 − CV at 2θ):
%       OSI = |Σ R(θ) · e^(i·2θ)| / Σ R(θ)
%   equivalent to equation (4) in Ringach et al. (2002) J Neurosci 22:5639
%   and equation (2) in Mazurek et al. (2014) Front Neural Circuits 8:130.
%
%   DSI is computed analogously at 1θ:
%       DSI = |Σ R(θ) · e^(i·θ)| / Σ R(θ)
%   which is the standard direction selectivity metric (Mazurek et al. 2014).
%
%   Preferred direction is the angle of the resultant vector at 1θ.
%   Preferred orientation is the angle of the resultant vector at 2θ,
%   divided by 2 to map back to orientation space.
%
%   INPUTS
%     vsObj  – any stimulus analysis object (e.g. linearlyMovingBallAnalysis,
%              StaticDriftingGratingAnalysis) that has a ResponseWindow method.
%     params – name-value pairs (see arguments block).
%
%   OUTPUTS
%     results – struct with fields:
%       .tuningCurve      [nN × nDir]  mean spike rate per direction
%       .tuningCurveSEM   [nN × nDir]  SEM across trials per direction
%       .OSI              [nN × 1]     orientation selectivity index
%       .DSI              [nN × 1]     direction selectivity index
%       .prefDirRad       [nN × 1]     preferred direction (radians)
%       .prefDirDeg       [nN × 1]     preferred direction (degrees)
%       .prefOriRad       [nN × 1]     preferred orientation (radians)
%       .prefOriDeg       [nN × 1]     preferred orientation (degrees)
%       .uDirRad          [1 × nDir]   unique directions in radians (sorted)
%       .uDirDeg          [1 × nDir]   unique directions in degrees (sorted)
%       .nTrialsPerDir    [nN × nDir]  trial count per neuron × direction
%       .goodU            [2 × nGood]  ic columns of good somatic units
%       .phyID            [nGood × 1]  phy cluster IDs of good units
%       .respMask         [nN × 1]     logical: neuron is responsive (p < threshold)
%       .pvals            [nN × 1]     p-values from StatisticsPerNeuron
%       .params           struct       copy of params used for this run
%
%   EXAMPLE
%     NP  = loadNPclassFromTable(49);
%     obj = StaticDriftingGratingAnalysis(NP);
%     res = DirectionTuning(obj, 'fieldName', 'Moving');
%
%   REFERENCES
%     Ringach DL et al. (2002) J Neurosci 22:5639-5651
%     Mazurek M et al. (2014) Front Neural Circuits 8:130
%
%   See also: AllExpDirectionTuning, hierBoot, plotSwarmBootstrapWithComparisons

% =========================================================================
%                         ARGUMENTS BLOCK
% =========================================================================
arguments
    vsObj                                                   % stimulus analysis object (handle class)

    % --- Which subfield of ResponseWindow to use ---
    params.fieldName        string  = "auto"                % e.g. 'Speed1', 'Moving'. "auto" attempts detection.

    % --- Which column in NeuronVals dim3 holds direction (radians) ---
    params.dirColumn        double  = NaN                   % NeuronVals dim3 index for direction.
                                                            % NaN = auto-detect from colNames.

    % --- Spike-rate column in NeuronVals dim3 ---
    params.rateColumn       double  = 1                     % NeuronVals dim3 index for mean spike rate.

    % --- Tuning curve construction ---
    params.aggMethod        string  = "mean"                % "mean" (standard) or "max" across trials
                                                            % per direction.  "mean" is recommended for
                                                            % OSI/DSI (Mazurek et al. 2014).

    % --- Responsiveness filter ---
    params.threshold        double  = 0.05                  % p-value threshold for responsive neurons

    % --- Statistics ---
    params.overwriteRW      logical = false                 % force recomputation of ResponseWindow
    params.overwriteStats   logical = false                 % force recomputation of StatisticsPerNeuron

    % --- Saving ---
    params.save             logical = true                  % save results to disk
    params.overwrite        logical = false                 % overwrite existing saved results
end

% =========================================================================
%  1.  LOAD RESPONSE WINDOW AND SPIKE SORTING
% =========================================================================

% Compute or load the ResponseWindow for this stimulus object
vsObj.ResponseWindow('overwrite', params.overwriteRW);      % ensures cached result exists
rw = vsObj.ResponseWindow;                                  % retrieve the cached struct

% Compute or load per-neuron statistics (p-values for responsiveness)
vsObj.StatisticsPerNeuron('overwrite', params.overwriteStats);
Stats = vsObj.StatisticsPerNeuron;                          % retrieve cached statistics

% =========================================================================
%  2.  RESOLVE FIELD NAME (Speed1 / Speed2 / Moving / Static / flat)
% =========================================================================

% Determine which subfield of the ResponseWindow struct contains NeuronVals
fieldName = resolveFieldName(rw, vsObj, params.fieldName);

% =========================================================================
%  3.  EXTRACT NEURONVALS AND IDENTIFY DIRECTION COLUMN
% =========================================================================

% NeuronVals dimensions: [nGoodUnits, nConditions, nFeatures]
%   dim3 = 1          : mean spike rate
%   dim3 = 2..4       : internal features (duration, baseline, etc.)
%   dim3 = 5+         : stimulus condition columns matching rw.colNames{1}(5:end)
NeuronVals = rw.(fieldName).NeuronVals;                     % [nN × nCond × nFeat]
nN   = size(NeuronVals, 1);                                 % number of good somatic units
nFeat = size(NeuronVals, 3);                                % total features in dim3

% Auto-detect which dim3 column holds direction (radians)
if isnan(params.dirColumn)
    dirCol = findDirectionColumn(rw, nFeat);                % returns dim3 index
else
    dirCol = params.dirColumn;                              % user override
end

% Validate that the detected column is within bounds
assert(dirCol >= 1 && dirCol <= nFeat, ...
    'DirectionTuning:badDirCol', ...
    'Direction column index %d is out of range [1, %d].', dirCol, nFeat);

fprintf('Using NeuronVals dim3=%d for direction labels, dim3=%d for spike rate.\n', ...
    dirCol, params.rateColumn);

% =========================================================================
%  4.  EXTRACT DIRECTION LABELS AND SPIKE RATES
% =========================================================================

% Direction labels per condition — same for all neurons, read from row 1.
% Each column of dim 2 corresponds to one unique stimulus condition (collapsed
% across trials by ResponseWindow).
dirLabels  = squeeze(NeuronVals(1, :, dirCol));             % [1 × nCond] direction in radians
spikeRates = NeuronVals(:, :, params.rateColumn);           % [nN × nCond] spike rate per condition

% Sorted unique directions (radians)
uDirRad = unique(dirLabels);                                % [1 × nDir]
nDir    = numel(uDirRad);                                   % number of unique directions

% Sanity check: at least 3 directions needed for meaningful OSI/DSI
assert(nDir >= 3, ...
    'DirectionTuning:tooFewDirs', ...
    'Only %d unique directions found — need >= 3 for OSI/DSI.', nDir);

fprintf('Found %d unique directions: [%s] (deg)\n', ...
    nDir, num2str(rad2deg(uDirRad'), '%.1f '));

% =========================================================================
%  5.  BUILD TUNING CURVE — MEAN (OR MAX) SPIKE RATE PER DIRECTION
% =========================================================================

% Pre-allocate tuning curve and SEM arrays
tuningCurve    = zeros(nN, nDir);                           % [nN × nDir] mean rate per direction
tuningCurveSEM = zeros(nN, nDir);                           % [nN × nDir] SEM per direction
nTrialsPerDir  = zeros(nN, nDir);                           % [nN × nDir] trial counts for diagnostics

for d = 1:nDir
    % Logical mask: which conditions share this direction value.
    % Using a tolerance of 1e-6 rad (~0.00006 deg) to handle floating-point
    % imprecision in stored radian values.
    dirMask = abs(dirLabels - uDirRad(d)) < 1e-6;          % [1 × nCond] logical

    % Number of conditions (not trials) matching this direction.
    % NOTE: If the stimulus has additional factors (size, lum, speed), each
    %       combination is a separate condition sharing the same direction.
    %       The tuning curve aggregates across ALL such conditions, which is
    %       correct for a marginal direction tuning curve.
    nCond = sum(dirMask);                                   % scalar

    % Extract spike rates for all neurons at this direction
    ratesAtDir = spikeRates(:, dirMask);                    % [nN × nCond]

    % Store trial count per neuron (same across neurons since conditions
    % are shared, but stored per-neuron for generality)
    nTrialsPerDir(:, d) = nCond;                            % replicate across neurons

    if strcmpi(params.aggMethod, "mean")
        % RECOMMENDED: Mean across conditions sharing this direction.
        % Standard practice for OSI/DSI computation (Mazurek et al. 2014).
        tuningCurve(:, d) = mean(ratesAtDir, 2);           % [nN × 1]
    elseif strcmpi(params.aggMethod, "max")
        % Alternative: max across conditions — biased toward outliers,
        % inflates tuning curve peaks. Use only if justified.
        tuningCurve(:, d) = max(ratesAtDir, [], 2);        % [nN × 1]
    else
        error('DirectionTuning:badAggMethod', ...
            'Unknown aggMethod "%s". Use "mean" or "max".', params.aggMethod);
    end

    % SEM across conditions for error bars.
    % If only one condition matches, SEM = 0.
    if nCond > 1
        tuningCurveSEM(:, d) = std(ratesAtDir, 0, 2) / sqrt(nCond);  % [nN × 1]
    else
        tuningCurveSEM(:, d) = 0;                          % undefined for n=1
    end
end

% =========================================================================
%  6.  RECTIFY NEGATIVE RATES (GUARD FOR BASELINE-SUBTRACTED DATA)
% =========================================================================

% If rates were baseline-subtracted, some may be negative, which would make
% OSI/DSI uninterpretable (denominator could be zero or negative).
% Clamp to zero with a warning.
hasNegative = any(tuningCurve(:) < 0);                      % check for any negative rates
if hasNegative
    warning('DirectionTuning:negativeRates', ...
        'Tuning curve contains negative rates (likely baseline-subtracted). ' + ...
        'Clamping to zero. Consider using raw (non-subtracted) spike rates.');
    tuningCurve = max(tuningCurve, 0);                      % clamp negatives to zero
end

% =========================================================================
%  7.  COMPUTE OSI (CIRCULAR VARIANCE AT 2θ)
% =========================================================================

% OSI = |Σ_d R(θ_d) · e^(i·2·θ_d)| / Σ_d R(θ_d)
%
% This is 1 minus the circular variance of the response distribution at
% double the angle (so that 0° and 180° map to the same orientation).
% Range: [0, 1] where 0 = untuned, 1 = perfectly orientation-selective.

% Sum of responses across all directions per neuron (denominator)
sumR = sum(tuningCurve, 2);                                 % [nN × 1]

% Numerator components: project responses onto sin(2θ) and cos(2θ)
sinComponent2 = sum(tuningCurve .* sin(2 * uDirRad), 2);   % [nN × 1]
cosComponent2 = sum(tuningCurve .* cos(2 * uDirRad), 2);   % [nN × 1]

% OSI = magnitude of resultant vector / sum of responses
% Guard against division by zero: if sumR == 0, neuron is unresponsive → NaN
OSI = sqrt(sinComponent2.^2 + cosComponent2.^2) ./ sumR;   % [nN × 1]

% =========================================================================
%  8.  COMPUTE DSI (CIRCULAR VARIANCE AT 1θ)
% =========================================================================

% DSI = |Σ_d R(θ_d) · e^(i·θ_d)| / Σ_d R(θ_d)
%
% Same formula as OSI but without angle-doubling, so opposite directions
% do NOT collapse. Range: [0, 1] where 0 = no direction preference,
% 1 = responds to only one direction.

% Numerator components: project responses onto sin(θ) and cos(θ)
sinComponent1 = sum(tuningCurve .* sin(uDirRad), 2);       % [nN × 1]
cosComponent1 = sum(tuningCurve .* cos(uDirRad), 2);       % [nN × 1]

% DSI = magnitude of resultant vector at 1θ / sum of responses
DSI = sqrt(sinComponent1.^2 + cosComponent1.^2) ./ sumR;   % [nN × 1]

% =========================================================================
%  9.  COMPUTE PREFERRED DIRECTION AND ORIENTATION
% =========================================================================

% Preferred direction: angle of the resultant vector at 1θ.
% atan2 returns values in [−π, +π]; wrap to [0, 2π) for convenience.
prefDirRad = atan2(sinComponent1, cosComponent1);           % [nN × 1] in (−π, π]
prefDirRad = mod(prefDirRad, 2*pi);                         % [nN × 1] in [0, 2π)
prefDirDeg = rad2deg(prefDirRad);                           % [nN × 1] in [0, 360)

% Preferred orientation: angle of the resultant vector at 2θ, halved.
% Since the angle was doubled, dividing by 2 maps back to orientation space.
% Result is in [0°, 180°).
prefOriRad = atan2(sinComponent2, cosComponent2);           % [nN × 1] resultant at 2θ
prefOriRad = mod(prefOriRad, 2*pi) / 2;                    % [nN × 1] in [0, π)
prefOriDeg = rad2deg(prefOriRad);                           % [nN × 1] in [0, 180)

% =========================================================================
%  10. EXTRACT P-VALUES AND RESPONSIVENESS MASK
% =========================================================================

% Get p-values for each neuron from StatisticsPerNeuron.
% The subfield structure depends on stimulus type (e.g., Stats.Moving,
% Stats.Speed1, or flat Stats.pvalsResponse).
pvals = extractPvals(Stats, fieldName);                     % [nN × 1]

% Logical mask: which neurons are responsive at the given threshold
respMask = pvals < params.threshold;                        % [nN × 1] logical

fprintf('Responsive neurons: %d / %d (p < %.3f)\n', ...
    sum(respMask), nN, params.threshold);

% =========================================================================
%  11. EXTRACT GOOD-UNIT IDENTIFIERS
% =========================================================================

% Get spike sorting info: ic columns and phy cluster IDs for good units
p_sort = vsObj.dataObj.convertPhySorting2tIc(vsObj.spikeSortingFolder, 0, 1, 1);
label  = string(p_sort.label');                             % label for each unit
goodU  = p_sort.ic(:, label == 'good');                     % [2 × nGood] spike train indices
phyID  = p_sort.phy_ID(label == 'good');                    % [nGood × 1] phy cluster IDs

% =========================================================================
%  12. ASSEMBLE OUTPUT STRUCT
% =========================================================================

results.tuningCurve    = tuningCurve;                       % [nN × nDir]
results.tuningCurveSEM = tuningCurveSEM;                    % [nN × nDir]
results.OSI            = OSI;                               % [nN × 1]
results.DSI            = DSI;                               % [nN × 1]
results.prefDirRad     = prefDirRad;                        % [nN × 1]
results.prefDirDeg     = prefDirDeg;                        % [nN × 1]
results.prefOriRad     = prefOriRad;                        % [nN × 1]
results.prefOriDeg     = prefOriDeg;                        % [nN × 1]
results.uDirRad        = uDirRad;                           % [1 × nDir]
results.uDirDeg        = rad2deg(uDirRad);                  % [1 × nDir]
results.nTrialsPerDir  = nTrialsPerDir;                     % [nN × nDir]
results.goodU          = goodU;                             % [2 × nGood]
results.phyID          = phyID;                             % [nGood × 1]
results.respMask       = respMask;                          % [nN × 1] logical
results.pvals          = pvals;                             % [nN × 1]
results.params         = params;                            % copy of params used

% =========================================================================
%  13. SAVE TO DISK
% =========================================================================

if params.save
    saveDir  = fileparts(vsObj.getAnalysisFileName);        % experiment analysis folder
    saveName = sprintf('DirectionTuning_%s.mat', fieldName);% e.g. DirectionTuning_Moving.mat
    savePath = fullfile(saveDir, saveName);                  % full path

    if ~exist(savePath, 'file') || params.overwrite
        save(savePath, '-struct', 'results');                % save all fields as top-level vars
        fprintf('Saved: %s\n', savePath);
    else
        fprintf('File exists (use overwrite=true to replace): %s\n', savePath);
    end
end

end   % ===== END OF MAIN FUNCTION =====


% =========================================================================
%                      LOCAL HELPER FUNCTIONS
% =========================================================================


function fieldName = resolveFieldName(rw, vsObj, requested)
% resolveFieldName  Determine the correct subfield of ResponseWindow.
%
%   "auto" logic:
%     - StaticDriftingGratingAnalysis → 'Moving' (drifting phase)
%     - linearlyMovingBall/Bar       → highest Speed field (e.g. Speed2)
%     - everything else              → '' (flat struct, no subfield)
%
%   If the user explicitly provides a field name, use it directly.

    if requested ~= "auto"
        % User provided an explicit field name — trust it
        fieldName = char(requested);                        % convert string to char
        assert(isfield(rw, fieldName), ...
            'DirectionTuning:badField', ...
            'Field "%s" not found in ResponseWindow struct.', fieldName);
        return
    end

    % Auto-detect based on stimulus class name
    className = class(vsObj);                               % e.g. 'StaticDriftingGratingAnalysis'

    if contains(className, 'StaticDrifting', 'IgnoreCase', true)
        % For SDG, use the Moving (drifting) phase for direction tuning
        fieldName = 'Moving';

    elseif contains(className, {'linearlyMovingBall', 'linearlyMovingBar'}, 'IgnoreCase', true)
        % For moving ball/bar, use the highest speed field available
        fn = fieldnames(rw);                                % all field names
        speedFields = fn(startsWith(fn, 'Speed'));          % e.g. {'Speed1','Speed2'}
        if isempty(speedFields)
            error('DirectionTuning:noSpeed', ...
                'No Speed fields found in ResponseWindow for %s.', className);
        end
        fieldName = speedFields{end};                       % highest speed (sorted alphabetically)

    else
        % Generic stimulus — assume flat struct with NeuronVals at top level.
        % Check common field names.
        if isfield(rw, 'NeuronVals')
            fieldName = '';                                 % will be handled downstream
        else
            % Try first available subfield that contains NeuronVals
            fn = fieldnames(rw);
            found = false;
            for i = 1:numel(fn)
                if isstruct(rw.(fn{i})) && isfield(rw.(fn{i}), 'NeuronVals')
                    fieldName = fn{i};
                    found = true;
                    break
                end
            end
            if ~found
                error('DirectionTuning:noNeuronVals', ...
                    'Cannot locate NeuronVals in ResponseWindow for %s.', className);
            end
        end
    end

    fprintf('Auto-detected fieldName: "%s" (class: %s)\n', fieldName, className);
end


function dirCol = findDirectionColumn(rw, nFeat)
% findDirectionColumn  Auto-detect the NeuronVals dim3 index for direction.
%
%   Strategy: look up 'direction' in rw.colNames{1}(5:end) (the condition
%   parameter names, stripped of the first 4 internal features: nSpks,
%   duration, spikeRate, baseline).  The position in the stripped list plus
%   4 gives the dim3 index in NeuronVals.
%
%   Falls back to searching all colNames if the first 4 are not internal.

    % Retrieve full column names from the ResponseWindow
    if isfield(rw, 'colNames') && ~isempty(rw.colNames)
        fullColNames = rw.colNames{1};                      % cell array of all column names
    else
        error('DirectionTuning:noColNames', ...
            'rw.colNames is missing or empty. Cannot auto-detect direction column. ' + ...
            'Set params.dirColumn manually.');
    end

    % The first 4 entries of colNames are internal features stored in
    % NeuronVals dim3 indices 1–4 (spike count, duration, rate, baseline).
    % Condition parameter names start at index 5.
    condNames = fullColNames(5:end);                        % condition parameter names only

    % Search for 'direction' (case-insensitive) in condition names
    dirIdx = find(strcmpi(condNames, 'direction'), 1);      % index into condNames

    if ~isempty(dirIdx)
        dirCol = dirIdx + 4;                                % offset by the 4 internal features
        fprintf('Auto-detected direction column: colNames{%d} = "%s" → NeuronVals dim3=%d\n', ...
            dirIdx + 4, condNames{dirIdx}, dirCol);
    else
        % Fallback: search the full colNames array (in case the 4-offset
        % assumption is wrong for this stimulus type)
        dirIdxFull = find(strcmpi(fullColNames, 'direction'), 1);
        if ~isempty(dirIdxFull)
            dirCol = dirIdxFull;
            warning('DirectionTuning:offsetMismatch', ...
                'Found "direction" at colNames{%d} without the expected 4-column offset. Verify mapping.', ...
                dirIdxFull);
        else
            error('DirectionTuning:noDirCol', ...
                'Could not find "direction" in colNames: {%s}.\nSet params.dirColumn manually.', ...
                strjoin(fullColNames, ', '));
        end
    end

    % Validate range
    assert(dirCol <= nFeat, ...
        'DirectionTuning:dirColOOB', ...
        'Detected direction column %d exceeds NeuronVals dim3 size (%d).', dirCol, nFeat);
end


function pvals = extractPvals(Stats, fieldName)
% extractPvals  Pull per-neuron p-values from the Statistics struct.
%
%   Handles subfield structures:
%     Stats.Moving.pvalsResponse  (SDG)
%     Stats.Speed1.pvalsResponse  (MB)
%     Stats.pvalsResponse         (flat, e.g. RG/NI)

    if ~isempty(fieldName) && isfield(Stats, fieldName) && ...
            isfield(Stats.(fieldName), 'pvalsResponse')
        pvals = Stats.(fieldName).pvalsResponse;            % subfield path
    elseif isfield(Stats, 'pvalsResponse')
        pvals = Stats.pvalsResponse;                        % flat path
    else
        warning('DirectionTuning:noPvals', ...
            'Cannot find pvalsResponse in Stats. Setting all p = 1 (no filter).');
        pvals = ones(size(Stats.ZScoreU));                  % conservative: mark all as non-responsive
    end

    pvals = pvals(:);                                       % ensure column vector
end