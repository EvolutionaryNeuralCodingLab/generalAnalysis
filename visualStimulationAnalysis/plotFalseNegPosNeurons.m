

function plotFalseNegPosNeurons(D, params)
% plotFalseNegPosNeurons  Plot and save the single-neuron raster for every
%   false-negative / false-positive candidate from findFalseNegAndPos.
%
%   For each row of D, the stimulus it was flagged for (FN_for and/or FP_for)
%   is drawn with that analysis object's own plotRaster, restricted to the
%   candidate's Phy cluster ID, and the figure is saved to params.outDir.
%
%   Usage:
%     D = findFalseNegAndPos("...Ex_1-97_n44_Raster_MB-RG_union.mat");
%     plotFalseNegPosNeurons(D);                 % all candidates
%     plotFalseNegPosNeurons(D(1:5,:));          % small test batch first

arguments
    D table                                                              % candidate table from findFalseNegAndPos (already candidate-only)
    params.outDir (1,1) string = "W:\Large_scale_mapping_NP\lizards\Combined_lizard_analysis\FPandFNplots" % save location
    params.format (1,1) string = "png"                                   % "png" (raster) or "pdf" (vector)
    params.closeFigs logical   = true                                    % close each figure after saving (avoid 100s open at once)
end

if ~exist(params.outDir, 'dir'); mkdir(params.outDir); end               % ensure the output folder exists

uex = unique(D.Exp(:))';                                                 % experiments that contain at least one candidate
fprintf('Plotting FN/FP rasters for %d candidate rows across %d experiments.\n', height(D), numel(uex));

for ex = uex                                                             % one experiment at a time -> load NP only once
    try
        NP = loadNPclassFromTable(ex);                                   % Neuropixels object for this experiment
    catch ME
        warning('Exp %d: could not load NP (%s) — skipping its candidates.', ex, ME.message);
        continue                                                         % skip every candidate from this experiment
    end

    objCache = containers.Map('KeyType','char','ValueType','any');       % build each stimulus object at most once per experiment
    rows = find(D.Exp == ex)';                                           % candidate rows belonging to this experiment

    for r = rows                                                         % each candidate neuron in this experiment
        phy   = D.PhyID(r);                                              % Phy cluster ID to plot
        tags  = ["FN","FP"];                                             % the two possible flag directions
        stims = [D.FN_for(r), D.FP_for(r)];                             % stimulus named by each flag ("" if not flagged that way)

        for ti = 1:2                                                     % FN first, then FP
            stim = stims(ti);                                            % stimulus for this flag direction
            if stim == ""; continue; end                                % candidate is not flagged in this direction

            % --- get (or build once) the analysis object for this stimulus ---
            if ~isKey(objCache, char(stim))                             % build only on first use this experiment
                try
                    objCache(char(stim)) = buildStimObj(NP, stim);      % stimulus-specific analysis object
                catch ME
                    warning('Exp %d %s: could not build object (%s) — skipping.', ex, stim, ME.message);
                    continue
                end
            end
            obj = objCache(char(stim));                                 % reuse cached object

            % --- plot the single neuron, capturing whatever figure(s) it opens ---
            before = findall(0,'Type','figure');                       % figures open before the call
            try
                plotOneRaster(obj, stim, phy);                         % framework plotRaster for this neuron (see helper)
            catch ME
                warning('Exp %d %s phy %g: plotRaster failed (%s) — skipping.', ex, stim, phy, ME.message);
                continue
            end
            newFigs = setdiff(findall(0,'Type','figure'), before);     % figure(s) created by the call
            if isempty(newFigs); newFigs = gcf; end                    % fallback: the method reused an existing figure

            % --- save each new figure with a descriptive, sortable name ---
            for fi = 1:numel(newFigs)                                   % usually one; RG may open several panels
                suffix = "";                                           % no suffix when a single figure
                if numel(newFigs) > 1; suffix = "_p" + fi; end         % disambiguate multiple panels
                fname = makeName(D, r, stim, tags(ti), phy, ex, suffix, params.format); % build the filename
                exportgraphics(newFigs(fi), fullfile(params.outDir, fname), 'Resolution', 300); % write the image
                if params.closeFigs; close(newFigs(fi)); end           % free the figure
            end
        end
    end
end

fprintf('Done. Images saved to:\n  %s\n', params.outDir);
end

% =========================================================================
% LOCAL HELPERS
% =========================================================================
function obj = buildStimObj(NP, stim)
% buildStimObj  Analysis object for a stimulus abbreviation.
    switch stim
        case "RG",   obj = rectGridAnalysis(NP);                        % rectified grid
        case "MB",   obj = linearlyMovingBallAnalysis(NP);              % moving ball
        case "MBR",  obj = linearlyMovingBarAnalysis(NP);              % moving bar
        case {"SDGs","SDGm"}, obj = StaticDriftingGratingAnalysis(NP);  % grating
        case "NI",   obj = imageAnalysis(NP);                          % natural image
        case "NV",   obj = movieAnalysis(NP);                          % natural movie
        case "FFF",  obj = fullFieldFlashAnalysis(NP);                 % full-field flash
        otherwise,   error('Unknown stimulus "%s".', stim);
    end
end

function plotOneRaster(obj, stim, phy)
% plotOneRaster  Call the framework's single-neuron raster for one Phy ID.
%   RG follows your existing call. MB is written by analogy — CONFIRM it
%   matches linearlyMovingBallAnalysis.plotRaster (speed arg? names?).
    switch stim
        case "RG"
            obj.plotRaster(exNeuronsPhyID=phy, AllResponsiveNeurons=false, ...
                MergeNtrials=1, selectedLum=[], PaperFig=false);        % all luminances, this neuron only
        case {"MB","MBR"}
            obj.plotRaster(exNeuronsPhyID=phy, AllResponsiveNeurons=false, PaperFig=false); % <-- CONFIRM signature
        otherwise
            obj.plotRaster(exNeuronsPhyID=phy, AllResponsiveNeurons=false, PaperFig=false); % generic fallback
    end
end

function fname = makeName(D, r, stim, tag, phy, ex, suffix, fmt)
% makeName  Windows-safe, sortable filename for one candidate.
    animal    = string(D.Animal(r));                                    % animal ID (may be "<unresolved>")
    insertion = D.Insertion(r);                                         % insertion number (may be NaN)
    verdict   = "NA";                                                   % default when Stage 2 (confirm) was off
    if ismember("Verdict", string(D.Properties.VariableNames))         % include the verdict if present
        verdict = string(D.Verdict(r));
    end
    raw = sprintf('%s_%s_ins%g_%s_%s_phy%g_ex%g%s', ...                 % verdict first so files group by it
        sanitize(verdict), sanitize(animal), insertion, stim, tag, phy, ex, suffix);
    fname = raw + "." + fmt;                                            % append the extension
end

function s = sanitize(str)
% sanitize  Replace characters Windows forbids in filenames.
    s = regexprep(string(str), '[<>:"/\\|?*; ]+', '-');                 % forbidden / spacing chars -> '-'
    if strlength(s) == 0; s = "NA"; end                                % never return an empty token
end