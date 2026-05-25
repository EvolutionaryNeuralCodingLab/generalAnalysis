function fig = AllExpAnalysisV2(expList, Stims2Comp, params)
% PlotZScoreComparison  Pool z-scores and spike rates across recordings,
%   run hierarchical bootstrapping, and produce swarm + scatter figures.
%
% KEY CHANGES FROM PREVIOUS VERSION
%   SUGG-1  The 7×3 analysis dispatch and 7×4 extraction blocks are replaced
%           by loops over a stimulus registry plus two local helpers:
%             runStimAnalysis(vsObj,presentKey,stimKey,params,Stims2Comp)
%             extractStimVals(stats,rw,stimKey,StatMethod)
%   SUGG-2  Absent-stimulus data is filled with NaN instead of -Inf, so
%           standard MATLAB functions (mean, max, isnan) work without extra
%           guard code.
%   SUGG-3  Optional BH-FDR correction via params.useFDR / params.FDRmethod.
%           Applied per-recording to the raw p-value vectors before thresholding.
%   SUGG-4  Spike-rate scatter axes can be log-scaled via params.logScaleSpkR.
%   SUGG-5  Scatter plots show ALL neurons (no randiColors subsampling).
%           Alpha is reduced to 0.25 to handle overplotting.
%   SUGG-6  ComparePairs table building now uses a stimLookup struct instead
%           of `who` + `eval`, making variable access explicit and debuggable.
%   SUGG-7  All accumulator cell arrays are pre-allocated before the loop.
%
% BUG FIXES (retained from documented review)
%   BUG-1   Guard against empty TableStimComp  (no responsive units crash).
%   BUG-2   fprintf moved after NP = loadNPclassFromTable(ex).
%   BUG-3   AnimalChanged flag evaluated before AnimalI is updated, so the
%           insertion counter uses the correct pre-update animal state.
%   BUG-4   sumNeurSDGm and sumNeurSDGst reset to 0 when SDG is absent.
%   BUG-5   2+2 debug stub replaced with warning().
%   BUG-6   S.groupStats.P_ZscoreCompare consistently nested (was top-level).

% =========================================================================
% ARGUMENT BLOCK
% =========================================================================
arguments
    expList    (1,:) double   % Experiment IDs from the master Excel table
    Stims2Comp cell           % Stimulus comparison order, e.g. {'MB','RG','MBR'}.
                              %   First element is the anchor for neuron selection.
    params.threshold         = 0.05    % p-value cut-off for responsiveness
    params.diffResp          = false   % Use spike-diff (resp-base) instead of rate
    params.overwrite         = false   % Force recompute of combined .mat file
    params.StimsPresent      = {'MB','RG'}  % Stimuli present in all sessions
    params.StimsNotPresent   = {}
    params.StimsToCompare    = {}      % Override scatter pair (default: 1st & 2nd)
    params.overwriteResponse = false   % Force ResponseWindow recompute
    params.overwriteStats    = false   % Force per-neuron stats recompute
    params.overwriteGroupStats = false % Force group bootstrap recompute
    params.RespDurationWin   = 100     % Response-window duration (ms)
    params.shuffles          = 2000    % Bootstrap / shuffle iterations (per-neuron)
    params.StatMethod        = 'ObsWindow'  % 'ObsWindow' | 'bootsrapRespBase' | 'maxPermuteTest'
    params.ignoreNonSignif   = false   % Zero z-scores of non-sig secondary stimuli
    params.EachStimSignif    = false   % Use per-stimulus responsive sets (not anchor)
    params.ComparePairs      = {}      % Pair(s) for focused pairwise comparison
    params.PaperFig logical  = false   % Save figures via vs.printFig
    % [SUGG-3] FDR options
    params.useFDR    logical = false   % Apply FDR correction to p-values per recording
    params.FDRmethod char    = 'BH'    % FDR method: 'BH' (Benjamini-Hochberg)
    % [SUGG-4] Spike-rate axis scaling
    params.logScaleSpkR logical = false % Log-scale spike-rate scatter axes
end

% =========================================================================
% SECTION 1 – INITIALISE BOOKKEEPING  [SUGG-7: full pre-allocation]
% =========================================================================
n = numel(expList);

animalVector      = cell(1, n);   % per-neuron animal index (anchor-filtered)
insertionVector   = cell(1, n);   % per-neuron insertion index
totalU            = cell(1, n);   % total unit count per recording
responsiveNeurons = cell(1, n);   % union of neuron indices responsive to any stim

animal     = 0;    % unique-animal counter
insertion  = 0;    % unique-insertion counter
AnimalI    = "";   % animal ID from previous iteration (for change detection)
InsertionI = 0;    % insertion number from previous iteration
j = 1;             % experiment counter (1-based index into pre-allocated cell arrays)

% [SUGG-1] Single organised store replaces 60+ individual cell arrays.
%   Anchor-filtered fields:  zScores, spKr, diffSpk, pvals
%   Self-responsive fields:  zScoresg, spKrg, diffSpkGCells, sumNeur
stimNames_all = {'MB','RG','MBR','FFF','SDGm','SDGs','NI','NV'};
for sni = stimNames_all
    sn = sni{1};
    zScoresCells.(sn)   = cell(1, n);  spKrCells.(sn)    = cell(1, n);
    diffSpkCells.(sn)   = cell(1, n);  pvalsCells.(sn)   = cell(1, n);
    zScoresgCells.(sn)  = cell(1, n);  spKrGCells.(sn)   = cell(1, n);
    diffSpkGCells.(sn)  = cell(1, n);  sumNeurCells.(sn) = cell(1, n);
end

% =========================================================================
% SECTION 2 – OUTPUT PATH AND SAVE-FILE CHECK
% =========================================================================
NP = loadNPclassFromTable(expList(1));   % load first recording for path extraction
vs = linearlyMovingBallAnalysis(NP);     % need analysis object for getAnalysisFileName

nameOfFile = sprintf('\\Ex_%d-%d_Combined_Neural_responses_%s_filtered.mat', ...
    expList(1), expList(end), Stims2Comp{1});

p = [extractBefore(vs.getAnalysisFileName, 'lizards'), 'lizards'];

if ~exist([p '\Combined_lizard_analysis'], 'dir')
    cd(p); mkdir Combined_lizard_analysis
end
saveDir = [p '\Combined_lizard_analysis'];

% Decide whether the for-loop is needed
if exist([saveDir nameOfFile], 'file') == 2 && ~params.overwrite
    S = load([saveDir nameOfFile]);
    forloop = ~isequal(S.expList, expList);   % reprocess if experiment list changed
else
    forloop = true;
end

% =========================================================================
% SECTION 3 – INITIALISE LONG-FORMAT TABLES
% =========================================================================

% Per-neuron table for pairwise z-score / spike-rate comparison
longTablePairComp = table( ...
    categorical.empty(0,1), categorical.empty(0,1), ...
    categorical.empty(0,1), categorical.empty(0,1), ...
    double.empty(0,1),      double.empty(0,1), ...
    'VariableNames', {'animal','insertion','stimulus','NeurID','Z-score','SpkR'});

% Per-insertion table for fraction-responsive analysis
longTable = table( ...
    categorical.empty(0,1), categorical.empty(0,1), categorical.empty(0,1), ...
    double.empty(0,1),      double.empty(0,1), ...
    'VariableNames', {'animal','insertion','stimulus','respNeur','totalSomaticN'});

% =========================================================================
% SECTION 4 – PER-EXPERIMENT FOR-LOOP
% =========================================================================
if forloop
    for ex = expList

        % [BUG-2 fix] Load recording BEFORE printing its name
        NP  = loadNPclassFromTable(ex);
        fprintf('Processing recording: %s .\n', NP.recordingName)

        % Create analysis objects for the two always-present stimuli
        vs  = linearlyMovingBallAnalysis(NP);   % Moving Ball  (MB)
        vsR = rectGridAnalysis(NP);             % Rect Grid    (RG)

        % Extract animal ID (try 'PV##' then 'SA##' as fallback)
        Animal = string(regexp(vs.getAnalysisFileName, 'PV\d+', 'match', 'once'));
        if isequal(Animal, "")
            Animal = string(regexp(vs.getAnalysisFileName, 'SA\d+', 'match', 'once'));
        end

        % Add rows to fraction-responsive table for always-present stimuli
        longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("MB"), 0, 0};
        longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("RG"), 0, 0};

        % ------------------------------------------------------------------
        % 4a – Load optional stimuli; fall back to a dummy when absent.
        %      Dummy objects have the same unit count so NaN replacement
        %      [SUGG-2] works correctly later.
        % ------------------------------------------------------------------
        try
            vsBr = linearlyMovingBarAnalysis(NP);
            params.StimsPresent{3} = 'MBR';
            if isempty(vsBr.VST), error('absent'); end
            longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("MBR"), 0, 0};
        catch
            params.StimsPresent{3} = '';
            fprintf('Moving Bar stimulus not found.\n')
            vsBr = linearlyMovingBallAnalysis(NP);   % dummy with correct N
        end

        try
            vsG = StaticDriftingGratingAnalysis(NP);
            params.StimsPresent{4} = 'SDG';
            if isempty(vsG.VST), error('absent'); end
            longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("SDGm"), 0, 0};
            longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("SDGs"), 0, 0};
        catch
            params.StimsPresent{4} = '';
            fprintf('Gratings stimulus not found.\n')
            vsG = rectGridAnalysis(NP);
        end

        try
            vsNI = imageAnalysis(NP);
            params.StimsPresent{5} = 'NI';
            if isempty(vsNI.VST), error('absent'); end
            longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("NI"), 0, 0};
        catch
            params.StimsPresent{5} = '';
            fprintf('Natural images not found.\n')
            vsNI = rectGridAnalysis(NP);
        end

        try
            vsNV = movieAnalysis(NP);
            params.StimsPresent{6} = 'NV';
            if isempty(vsNV.VST), error('absent'); end
            longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("NV"), 0, 0};
        catch
            params.StimsPresent{6} = '';
            fprintf('Natural video not found.\n')
            vsNV = rectGridAnalysis(NP);
        end

        try
            vsFFF = fullFieldFlashAnalysis(NP);
            params.StimsPresent{7} = 'FFF';
            if isempty(vsFFF.VST), error('absent'); end
            longTable(end+1,:) = {categorical(Animal), categorical(j), categorical("FFF"), 0, 0};
        catch
            params.StimsPresent{7} = '';
            fprintf('Full-field flash not found.\n')
            vsFFF = rectGridAnalysis(NP);
        end

        % ------------------------------------------------------------------
        % 4b – Run analyses.  [SUGG-1] replaces 7 × ~8-line conditional
        %      blocks with a single loop and one local helper.
        % ------------------------------------------------------------------
        vsObjs   = { vs,    vsR,   vsBr,   vsG,    vsNI,   vsNV,   vsFFF  };
        presKeys = params.StimsPresent(1:7);
        % 'SDG' key covers vsG for both SDGm and SDGs extraction below
        stimKeys = {'MB',  'RG',  'MBR',  'SDG',  'NI',   'NV',   'FFF'  };

        statsAll = cell(1, 7);
        rwAll    = cell(1, 7);
        for k = 1:7
            % runStimAnalysis handles ResponseWindow + statistics dispatch,
            % computes or loads from cache based on presence and overwrite flags.
            [statsAll{k}, rwAll{k}] = runStimAnalysis( ...
                vsObjs{k}, presKeys{k}, stimKeys{k}, params, Stims2Comp);
        end

        % ------------------------------------------------------------------
        % 4c – Extract z-scores, p-values, spike rate, spike diff.
        %      [SUGG-1] extractStimVals handles Speed/Moving/Static subfields.
        %      All outputs are column vectors of length N (unit count).
        % ------------------------------------------------------------------
        [zS.MB,   pV.MB,   spkR.MB,   spkDiff.MB  ] = extractStimVals(statsAll{1}, rwAll{1}, 'MB',   params.StatMethod);
        [zS.RG,   pV.RG,   spkR.RG,   spkDiff.RG  ] = extractStimVals(statsAll{2}, rwAll{2}, 'RG',   params.StatMethod);
        [zS.MBR,  pV.MBR,  spkR.MBR,  spkDiff.MBR ] = extractStimVals(statsAll{3}, rwAll{3}, 'MBR',  params.StatMethod);
        [zS.SDGm, pV.SDGm, spkR.SDGm, spkDiff.SDGm] = extractStimVals(statsAll{4}, rwAll{4}, 'SDGm', params.StatMethod);
        [zS.SDGs, pV.SDGs, spkR.SDGs, spkDiff.SDGs] = extractStimVals(statsAll{4}, rwAll{4}, 'SDGs', params.StatMethod);
        [zS.NI,   pV.NI,   spkR.NI,   spkDiff.NI  ] = extractStimVals(statsAll{5}, rwAll{5}, 'NI',   params.StatMethod);
        [zS.NV,   pV.NV,   spkR.NV,   spkDiff.NV  ] = extractStimVals(statsAll{6}, rwAll{6}, 'NV',   params.StatMethod);
        [zS.FFF,  pV.FFF,  spkR.FFF,  spkDiff.FFF ] = extractStimVals(statsAll{7}, rwAll{7}, 'FFF',  params.StatMethod);

        % Total units in this recording (before any filtering)
        totalU{j} = numel(zS.MB);

        % ------------------------------------------------------------------
        % 4d – [SUGG-3] Optional FDR correction on raw p-values.
        %      Applied per-recording before the significance threshold is used.
        %      Corrects for the number of neurons tested simultaneously.
        % ------------------------------------------------------------------
        if params.useFDR
            for sni = stimNames_all
                pV.(sni{1}) = bhFDR(pV.(sni{1}));
            end
        end

        % ------------------------------------------------------------------
        % 4e – [SUGG-2] Set absent-stimulus data to NaN.
        %      NaN propagates cleanly through isnan/nanmean/nanmax and
        %      integrates with the significance masks below without extra
        %      -Inf guard code throughout.
        % ------------------------------------------------------------------
        %   Mapping: {StimsPresent index, stimulus field name(s)}
        absentMap = { 2, {'RG'};   3, {'MBR'};  4, {'SDGm','SDGs'}; ...
                      5, {'NI'};   6, {'NV'};    7, {'FFF'} };
        N = numel(zS.MB);   % total units (all same length by construction)
        for ai = 1:size(absentMap, 1)
            if isequal(params.StimsPresent{absentMap{ai,1}}, '')
                for sni = absentMap{ai,2}
                    sn = sni{1};
                    zS.(sn)      = nan(N, 1);
                    pV.(sn)      = nan(N, 1);
                    spkR.(sn)    = nan(N, 1);
                    spkDiff.(sn) = nan(N, 1);
                end
            end
        end

        % ------------------------------------------------------------------
        % 4f – Optional: suppress z-scores of non-significant secondary neurons
        %      (only meaningful when comparing across stimuli with a shared anchor)
        % ------------------------------------------------------------------
        if params.ignoreNonSignif
            for sni = stimNames_all
                zS.(sni{1})(pV.(sni{1}) > params.threshold) = NaN;
            end
        end

        % ------------------------------------------------------------------
        % 4g – Determine the anchor p-value vector.
        %      [SUGG-6] Direct struct field access replaces the eval-based
        %      variable-name lookup used in the original code.
        % ------------------------------------------------------------------
        anchorField = Stims2Comp{1};
        if strcmp(anchorField, 'SDG'), anchorField = 'SDGm'; end  % SDG → moving
        pvalsAnchor = pV.(anchorField);
        anchorMask  = ~isnan(pvalsAnchor) & (pvalsAnchor <= params.threshold);

        % ------------------------------------------------------------------
        % 4h – Populate pairwise comparison table.
        %      [SUGG-6] Uses stimLookup struct; no eval/who required.
        % ------------------------------------------------------------------
        if ~isempty(params.ComparePairs)
            cp1 = params.ComparePairs{1};
            cp2 = params.ComparePairs{2};

            % A neuron enters the table if it is significant for EITHER stimulus
            sigMask = (~isnan(pV.(cp1)) & pV.(cp1) < params.threshold) | ...
                      (~isnan(pV.(cp2)) & pV.(cp2) < params.threshold);
            unitIDs = find(sigMask);

            if ~isempty(unitIDs)
                nu  = numel(unitIDs);
                zC1 = zS.(cp1)(sigMask);      rC1 = spkR.(cp1)(sigMask);
                zC2 = zS.(cp2)(sigMask);      rC2 = spkR.(cp2)(sigMask);
                repAnimal = categorical(cellstr(repmat(Animal, nu, 1)));
                repInser  = categorical(repmat(j, nu, 1));

                T1 = table(repAnimal, repInser, ...
                    categorical(cellstr(repmat(cp1, nu, 1))), categorical(unitIDs)', ...
                    zC1', rC1', ...
                    'VariableNames', {'animal','insertion','stimulus','NeurID','Z-score','SpkR'});
                T2 = table(repAnimal, repInser, ...
                    categorical(cellstr(repmat(cp2, nu, 1))), categorical(unitIDs)', ...
                    zC2', rC2', ...
                    'VariableNames', {'animal','insertion','stimulus','NeurID','Z-score','SpkR'});

                longTablePairComp = [longTablePairComp; T1; T2];
            end
        end

        % ------------------------------------------------------------------
        % 4i – Filter data and count responsive neurons per stimulus.
        %      Two subsets per stimulus:
        %        'anchor-filtered' : neurons sig. for the anchor stimulus
        %        'self-responsive' : neurons sig. for this stimulus itself
        % ------------------------------------------------------------------
        respIndexes = [];   % accumulates union of responsive indices

        for sni = stimNames_all
            sn = sni{1};

            % Anchor-filtered subset (all stimuli indexed by the same mask)
            zScoresCells.(sn){j}  = zS.(sn)(anchorMask);
            spKrCells.(sn){j}     = spkR.(sn)(anchorMask)';    % row vector (matches original packing)
            diffSpkCells.(sn){j}  = spkDiff.(sn)(anchorMask);
            pvalsCells.(sn){j}    = pV.(sn)(anchorMask);

            % Self-responsive subset
            selfMask = ~isnan(pV.(sn)) & (pV.(sn) <= params.threshold);
            zScoresgCells.(sn){j}  = zS.(sn)(selfMask);
            spKrGCells.(sn){j}     = spkR.(sn)(selfMask);      % column vector
            diffSpkGCells.(sn){j}  = spkDiff.(sn)(selfMask);
            sumNeurCells.(sn){j}   = sum(selfMask);
            respIndexes = [respIndexes, find(selfMask)];

            % Update fraction-responsive table if this insertion/stimulus row exists
            try
                idx = (longTable.insertion == categorical(j)) & ...
                      (longTable.stimulus  == categorical(sn));
                if any(idx)
                    longTable.respNeur(idx)      = sum(selfMask);
                    longTable.totalSomaticN(idx) = N;   % same denominator for all rows
                end
            end
        end

        responsiveNeurons{j} = unique(respIndexes);

        % [BUG-5 fix] Sanity check: NI and NV filtered lengths must match
        nNI = sum(~isnan(pvalsCells.NI{j}));
        nNV = sum(~isnan(pvalsCells.NV{j}));
        if nNI ~= nNV
            warning('PlotZScoreComparison:sizeMismatch', ...
                'NI and NV anchor-filtered lengths differ (%d vs %d) in experiment %d.', ...
                nNI, nNV, ex);
        end

        % ------------------------------------------------------------------
        % 4j – Animal / insertion change detection.
        %      [BUG-3 fix] AnimalChanged is evaluated BEFORE AnimalI is
        %      updated, so the insertion counter can use the correct
        %      pre-update state.
        % ------------------------------------------------------------------
        Insertion = str2double(regexp( ...
            regexp(vs.getAnalysisFileName, 'Insertion\d+', 'match', 'once'), '\d+', 'match'));
        if isempty(Insertion), Insertion = 0; end   % safety for missing tag

        AnimalChanged = (Animal ~= AnimalI);   % evaluate BEFORE modifying AnimalI
        if AnimalChanged
            animal = animal + 1;
            AnimalNames{animal} = Animal;      %#ok<AGROW>
            AnimalI = Animal;
        end
        if Insertion ~= InsertionI || AnimalChanged   % [BUG-3 fix] use pre-evaluated flag
            InsertionI = Insertion;
            insertion  = insertion + 1;
        end

        % Replicate animal / insertion IDs once per anchor-filtered neuron
        nAnchor          = sum(anchorMask);
        animalVector{j}    = repmat(animal,    [1, nAnchor]);
        insertionVector{j} = repmat(insertion, [1, nAnchor]);

        j = j + 1;
        fprintf('Finished recording: %s .\n', NP.recordingName)
    end

    % ======================================================================
    % SECTION 4-END: PACK STRUCT S AND SAVE
    % Maintain original field-name conventions (prefixes: zScores, spKr,
    % diffSpk, spkR, spkDiff, sumNeur) for backward compatibility with any
    % code that loads the saved .mat file.
    % ======================================================================
    for sni = stimNames_all
        sn = sni{1};
        % Anchor-filtered (one cell per recording)
        S.stimValsSignif2oneStim.(['zScores' sn])  = zScoresCells.(sn);
        S.stimValsSignif2oneStim.(['spKr'    sn])  = spKrCells.(sn);
        S.stimValsSignif2oneStim.(['diffSpk' sn])  = diffSpkCells.(sn);
        S.pvals.(['pvals' sn])                     = pvalsCells.(sn);
        % Self-responsive (note different capitalisation to match original patterns)
        S.stimValsSignif.(['zScores' sn 'g'])      = zScoresgCells.(sn);
        S.stimValsSignif.(['spkR'    sn 'g'])      = spKrGCells.(sn);
        S.stimValsSignif.(['spkDiff' sn 'g'])      = diffSpkGCells.(sn);
        S.stimValsSignif.(['sumNeur' sn])          = sumNeurCells.(sn);
    end

    S.expList           = expList;
    S.animalVector      = animalVector;
    S.insertionVector   = insertionVector;
    S.totalUnits        = totalU;
    S.params            = params;
    S.responsiveNeurons = responsiveNeurons;
    S.TableRespNeurs    = longTable;
    S.TableStimComp     = longTablePairComp;

    save([saveDir nameOfFile], '-struct', 'S');
end   % if forloop

% =========================================================================
% SECTION 5 – PAIRWISE COMPARISON MODE  (ComparePairs is non-empty)
% =========================================================================
if ~isempty(params.ComparePairs)

    pairs = params.ComparePairs;   % cell of stimulus names; rows = pairs

    % [BUG-1 fix] Guard: splitapply crashes on an empty grouping vector
    if isempty(S.TableStimComp) || height(S.TableStimComp) == 0
        warning('PlotZScoreComparison:noUnits', ...
            'No significant units found for %s vs %s. Returning empty figure.', ...
            pairs{1}, pairs{2});
        fig = figure; return
    end

    % Treat residual NaN values as zero (conservative: no response)
    S.TableStimComp.('Z-score')(isnan(S.TableStimComp.('Z-score'))) = 0;
    S.TableStimComp.SpkR(isnan(S.TableStimComp.SpkR)) = 0;

    % Keep only insertions that contain both stimuli in every pair
    [G, ~]  = findgroups(S.TableStimComp.insertion);
    hasAll  = splitapply( ...
        @(s) all(ismember(unique(categorical(pairs)), s)), ...
        S.TableStimComp.stimulus, G);

    nBoot = 10000;   % bootstrap iterations for hierarchical resampling

    % Reload vs for figure saving (path may be stale if forloop did not run)
    NP = loadNPclassFromTable(expList(1));
    vs = linearlyMovingBallAnalysis(NP);

    % Build display labels once (RG→SB, SDGs→SG, SDGm→MG)
    s = replace(replace(replace(string(pairs), "RG","SB"), "SDGs","SG"), "SDGm","MG");

    % ------------------------------------------------------------------
    % 5a – Z-score: hierarchical bootstrap + swarm
    % ------------------------------------------------------------------
    j = 1;
    ps = zeros(1, size(pairs, 1));

    for i = 1:size(pairs, 1)
        [diffs, insers, animals] = collectPairDiffs( ...
            S.TableStimComp, pairs, j, 'Z-score');
        bootDiff = hierBoot(diffs, nBoot, insers, animals);
        ps(j)    = mean(bootDiff <= 0);
        j        = j + 1;
    end

    % Swarm; discard randiColors — scatter will use all points [SUGG-5]
    [fig, ~] = plotSwarmBootstrapWithComparisons(S.TableStimComp, pairs, ps, {'Z-score'}, ...
        yLegend='Z-score', yMaxVis=40, diff=true, plotMeanSem=false, Alpha=0.7);
    applyPaperAxes(gca);
    set(fig, 'Units','centimeters', 'Position',[20 20 4 6]);
    if params.PaperFig
        vs.printFig(fig, sprintf('Zscore-comparison-Swarm-%s-%s', pairs{1}, pairs{2}), ...
            PaperFig=params.PaperFig);
    end

    % ------------------------------------------------------------------
    % 5b – Z-score scatter  [SUGG-5: all points; alpha=0.25 for overlap]
    % ------------------------------------------------------------------
    fig = figure;
    [p1z, p2z, cA] = getPairScatterData(S.TableStimComp, pairs, 'Z-score');
    scatter(p1z, p2z, 7, cA, 'filled', 'MarkerFaceAlpha', 0.25)
    hold on; axis equal
    lims = [-5 40];   ylim(lims); xlim(lims)
    plot(lims, lims, 'k--', 'LineWidth', 1.5)
    xlabel(s{1}); ylabel(s{2})
    colormap(lines(numel(categories(S.TableStimComp.animal))))
    applyPaperAxes(gca);  title('Z-score')
    set(fig, 'Units','centimeters', 'Position',[20 20 5 5]);
    if params.PaperFig
        vs.printFig(fig, sprintf('Zscore-comparison-Scatter-%s-%s', pairs{1}, pairs{2}), ...
            PaperFig=params.PaperFig);
    end

    % ------------------------------------------------------------------
    % 5c – Spike rate: hierarchical bootstrap + swarm
    % ------------------------------------------------------------------
    j = 1;
    ps = zeros(1, size(pairs, 1));

    for i = 1:size(pairs, 1)
        [diffs, insers, animals] = collectPairDiffs( ...
            S.TableStimComp, pairs, j, 'SpkR');
        bootDiff = hierBoot(diffs, nBoot, insers, animals);
        ps(j)    = mean(bootDiff <= 0);
        j        = j + 1;
    end

    V1max = max(abs(diffs));   % set y-ceiling from actual data range
    [fig, ~] = plotSwarmBootstrapWithComparisons(S.TableStimComp, pairs, ps, {'SpkR'}, ...
        yLegend='SpkR', yMaxVis=V1max, diff=true, plotMeanSem=false, Alpha=0.7);
    applyPaperAxes(gca);
    set(fig, 'Units','centimeters', 'Position',[20 20 4 6]);
    if params.PaperFig
        vs.printFig(fig, sprintf('spkRate-comparison-Swarm-%s-%s', pairs{1}, pairs{2}), ...
            PaperFig=params.PaperFig);
    end

    % ------------------------------------------------------------------
    % 5d – Spike-rate scatter  [SUGG-4: optional log scale; SUGG-5: all pts]
    % ------------------------------------------------------------------
    fig = figure;
    [p1r, p2r, cA] = getPairScatterData(S.TableStimComp, pairs, 'SpkR');
    scatter(p1r, p2r, 7, cA, 'filled', 'MarkerFaceAlpha', 0.25)   % all points
    hold on; axis equal
    posVals = S.TableStimComp.SpkR(S.TableStimComp.SpkR > 0);
    if ~isempty(posVals)
        lims = [min(posVals), max(posVals)];
    else
        lims = [0, 1];
    end
    if params.logScaleSpkR && all(lims > 0)   % [SUGG-4]
        set(gca, 'XScale','log', 'YScale','log')
    end
    plot(lims, lims, 'k--', 'LineWidth', 1.5)
    ylim(lims); xlim(lims)
    xlabel(s{1}); ylabel(s{2})
    colormap(lines(numel(categories(S.TableStimComp.animal))))
    applyPaperAxes(gca);  title('Spk. rate')
    set(fig, 'Units','centimeters', 'Position',[20 20 5 5]);
    if params.PaperFig
        vs.printFig(fig, sprintf('spkRate-comparison-Scatter-%s-%s', pairs{1}, pairs{2}), ...
            PaperFig=params.PaperFig);
    end

else
    % ======================================================================
    % SECTION 5-ALT – MULTI-STIMULUS OVERVIEW  (no ComparePairs)
    %   Four-panel figure: Z-score swarm, Z-score scatter,
    %                      spike-rate swarm, spike-rate scatter.
    % ======================================================================
    fig = figure;
    tiledlayout(2, 2, 'TileSpacing', 'compact');

    % Choose field set based on filtering mode
    if ~params.EachStimSignif
        fn = fieldnames(S.stimValsSignif2oneStim);
    else
        fn = fieldnames(S.stimValsSignif);
    end
    fnp = fieldnames(S.pvals);

    % Expand 'SDG' into moving + static sub-conditions
    Stims2Comp2 = {};
    for i = 1:numel(Stims2Comp)
        if strcmp(Stims2Comp{i}, 'SDG')
            Stims2Comp2 = [Stims2Comp2, {'SDGs','SDGm'}]; %#ok<AGROW>
        else
            Stims2Comp2 = [Stims2Comp2, Stims2Comp(i)];   %#ok<AGROW>
        end
    end

    % Field suffix: '' for anchor-filtered, 'g' for self-responsive
    ending2 = repmat({'','g'}, 1, 2);
    ending2 = ending2{1 + params.EachStimSignif};

    % Assemble concatenated data arrays for swarm and scatter
    StimZS    = cell(numel(Stims2Comp2), 1);
    stimRSP   = cell(numel(Stims2Comp2), 1);
    stimPvals = cell(numel(Stims2Comp2), 1);
    x         = [];   % stimulus index label per neuron

    for i = 1:numel(Stims2Comp2)
        ending = Stims2Comp2{i};

        % Z-scores
        pattern = ['^zS.*' ending ending2 '$'];
        matches = fn(~cellfun('isempty', regexp(fn, pattern)));
        if ~params.EachStimSignif
            StimZS{i} = cell2mat(S.stimValsSignif2oneStim.(matches{1}))';
        else
            StimZS{i} = cell2mat(S.stimValsSignif.(matches{1}))';
        end

        % Spike rates (or diff-response if diffResp flag is set)
        if ~params.diffResp
            pattern = ['^spKr.*' ending ending2 '$'];
        else
            pattern = ['^diffSpk.*' ending ending2 '$'];
        end
        matches = fn(~cellfun('isempty', regexp(fn, pattern)));
        if params.EachStimSignif
            matches   = fn(~cellfun('isempty', regexp(fn, pattern, 'ignorecase')));
            C         = S.stimValsSignif.(matches{1});
            C         = cellfun(@(x) x', C, 'UniformOutput', false);
            stimRSP{i} = cell2mat(C');
        else
            try
                stimRSP{i} = cell2mat(S.stimValsSignif2oneStim.(matches{1})');
            catch
                try
                    stimRSP{i} = cell2mat(S.stimValsSignif2oneStim.(matches{1}));
                catch
                    Ccol       = cellfun(@(x) x(:), S.stimValsSignif2oneStim.(matches{1}), ...
                                         'UniformOutput', false);
                    stimRSP{i} = vertcat(Ccol{:})';
                end
            end
        end

        % p-values (for completeness; not plotted directly here)
        pattern      = ['^pvals.*' ending '$'];
        matches      = fnp(~cellfun('isempty', regexp(fnp, pattern)));
        stimPvals{i} = cell2mat(S.pvals.(matches{1}))';

        x = [x; ones(size(StimZS{i})) * i]; %#ok<AGROW>
    end

    % Per-neuron animal / insertion labels for colouring
    AnIndex      = cell2mat(S.animalVector)';
    InsIndex     = cell2mat(S.insertionVector)';
    colormapUsed = parula(max(AnIndex)) .* 0.6;
    allColorIdx  = repmat(AnIndex, numel(Stims2Comp2), 1);

    % ------------------------------------------------------------------
    % Panel 1: Z-score swarm
    % ------------------------------------------------------------------
    y = cell2mat(StimZS);
    nexttile
    if ~params.EachStimSignif
        swarmchart(x, y, 5, colormapUsed(allColorIdx,:), 'filled', 'MarkerFaceAlpha', 0.7);
    else
        swarmchart(x, y, 5, 'filled', 'MarkerFaceAlpha', 0.7);
    end
    xticks(1:8); xticklabels(Stims2Comp2); ylabel('Z-score');
    set(fig, 'Color','w'); yline(0, 'LineWidth', 2); ylim([-5 40])

    % Z-score group-level hierarchical bootstrap
    if params.overwriteGroupStats || ~isfield(S, 'groupStats')
        [probs, ps_boot] = runGroupBoot(y, x, InsIndex, AnIndex, numel(Stims2Comp2));
        S.groupStats.Bayes_ZscoreCompare = probs;
        S.groupStats.P_ZscoreCompare     = ps_boot;   % [BUG-6 fix: nested correctly]
        save([saveDir nameOfFile], '-struct', 'S');
    end

    % ------------------------------------------------------------------
    % Panel 2: Z-score scatter  [SUGG-5: all points, alpha=0.25]
    % ------------------------------------------------------------------
    nexttile
    if isempty(params.StimsToCompare)
        ind1 = 1; ind2 = 2;
    else
        ind1 = find(strcmp(Stims2Comp2, params.StimsToCompare{1}));
        ind2 = find(strcmp(Stims2Comp2, params.StimsToCompare{2}));
    end
    VTC = {StimZS{ind1}, StimZS{ind2}};
    if numel(VTC{1}) == numel(VTC{2})
        scatter(VTC{1}, VTC{2}, 10, AnIndex, 'filled', 'MarkerFaceAlpha', 0.25)
        colormap(colormapUsed); hold on; axis equal
        validY = y(~isnan(y) & ~isinf(y));
        lims   = [min(validY), max(validY)];
        plot(lims, lims, 'k--', 'LineWidth', 1.5)
        lims = [-5 40]; ylim(lims); xlim(lims)
        xlabel(Stims2Comp(ind1)); ylabel(Stims2Comp(ind2))
    end

    % ------------------------------------------------------------------
    % Panel 3: Spike-rate swarm
    % ------------------------------------------------------------------
    y = cell2mat(stimRSP);
    nexttile
    if ~params.EachStimSignif
        swarmchart(x, y, 5, colormapUsed(allColorIdx,:), 'filled', 'MarkerFaceAlpha', 0.7);
    else
        swarmchart(x, y, 5, 'filled', 'MarkerFaceAlpha', 0.7);
    end
    xticks(1:8); xticklabels(Stims2Comp2); ylabel('Spike Rate'); set(fig,'Color','w')

    % Spike-rate group-level hierarchical bootstrap
    if params.overwriteGroupStats || ~isfield(S, 'groupStats')
        [probs, ps_boot] = runGroupBoot(y, x, InsIndex, AnIndex, numel(Stims2Comp2));
        S.groupStats.Bayes_SpikeRateCompare = probs;
        S.groupStats.P_SpikeRateCompare     = ps_boot;
    end

    % ------------------------------------------------------------------
    % Panel 4: Spike-rate scatter  [SUGG-4: log scale; SUGG-5: all pts]
    % ------------------------------------------------------------------
    nexttile
    VTC = {stimRSP{ind1}, stimRSP{ind2}};
    if numel(VTC{1}) == numel(VTC{2})
        scatter(VTC{1}, VTC{2}, 10, AnIndex, 'filled', 'MarkerFaceAlpha', 0.25)
        colormap(colormapUsed); hold on; axis equal
        posY   = y(y > 0 & ~isinf(y) & ~isnan(y));
        lims   = [min(posY), max(posY)];
        if params.logScaleSpkR && all(lims > 0)    % [SUGG-4]
            set(gca, 'XScale','log', 'YScale','log')
        end
        plot(lims, lims, 'k--', 'LineWidth', 1.5)
        ylim(lims); xlim(lims)
        xlabel(Stims2Comp(ind1)); ylabel(Stims2Comp(ind2))
    end

end   % if/else ComparePairs

% =========================================================================
% SECTION 6 – FRACTION-RESPONSIVE ANALYSIS
% =========================================================================
if isempty(params.ComparePairs)
    pairs = {Stims2Comp{1}, Stims2Comp{2}};
else
    pairs = params.ComparePairs;
end

% Insertions with both stimuli present
[G, ~]    = findgroups(S.TableRespNeurs.insertion);
hasAll    = splitapply( ...
    @(s) all(ismember(unique(categorical(pairs)), s)), ...
    S.TableRespNeurs.stimulus, G);
tempTable = S.TableRespNeurs( ...
    hasAll(G) & ismember(S.TableRespNeurs.stimulus, unique(categorical(pairs))), :);

nBoot = 10000;
j     = 1;
ps    = zeros(1, size(pairs, 1));

for i = 1:size(pairs, 1)
    diffs = [];
    for ins = unique(S.TableRespNeurs.insertion)'
        idx1 = S.TableRespNeurs.insertion==categorical(ins) & S.TableRespNeurs.stimulus==pairs{j,1};
        idx2 = S.TableRespNeurs.insertion==categorical(ins) & S.TableRespNeurs.stimulus==pairs{j,2};
        if any(idx1) && any(idx2)
            tot  = S.TableRespNeurs.totalSomaticN(idx1);   % shared denominator
            diffs(end+1, 1) = S.TableRespNeurs.respNeur(idx1)/tot - ...
                              S.TableRespNeurs.respNeur(idx2)/tot; %#ok<AGROW>
        end
    end
    bootDiff = bootstrp(nBoot, @mean, diffs);
    ps(j)    = mean(bootDiff <= 0);
    j        = j + 1;
end

% Add total-responsive column for the plotting helper
[G, ~]                    = findgroups(tempTable.insertion);
totals                    = splitapply(@sum, tempTable.respNeur, G);
tempTable.TotalRespNeur   = totals(G);

fig = plotSwarmBootstrapWithComparisons(tempTable, pairs, ps, ...
    {'respNeur','totalSomaticN'}, fraction=true, ...
    yLegend='Responsive/total units', diff=false, filled=false, ...
    Xjitter='none', Alpha=0.6);

applyPaperAxes(gca);
set(fig, 'Units','centimeters', 'Position',[20 20 5 6]);

if params.PaperFig && ~isempty(params.ComparePairs)
    vs.printFig(fig, sprintf('ResponsiveUnits-comparison-%s-%s', ...
        params.ComparePairs{1}, params.ComparePairs{2}), PaperFig=params.PaperFig);
end

end   % ===== END OF MAIN FUNCTION =====


% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

% -------------------------------------------------------------------------
function [stats, rw] = runStimAnalysis(vsObj, presentKey, stimKey, params, Stims2Comp)
% runStimAnalysis  Compute or load statistics for one stimulus.
%
%   If presentKey is non-empty AND the stimulus appears in Stims2Comp the
%   ResponseWindow and the selected statistics method are (re-)computed.
%   Otherwise only the cached result is loaded.  Both branches return the
%   same output types, so callers need no conditional logic.
%
%   INPUTS
%     vsObj      – analysis object (e.g. linearlyMovingBallAnalysis)
%     presentKey – string key in params.StimsPresent, or '' if absent
%     stimKey    – canonical stimulus name ('MB','RG',…) for display only
%     params     – parameter struct from the main function arguments block
%     Stims2Comp – cell of stimulus names that should be analysed
%
%   OUTPUTS
%     stats – statistics struct (ShufflingAnalysis / Bootstrap / Statistics)
%     rw    – ResponseWindow struct

    shouldCompute = ~isequal(presentKey, '') && ismember(presentKey, Stims2Comp);

    if shouldCompute
        vsObj.ResponseWindow('overwrite', params.overwriteResponse, ...
                             'durationWindow', params.RespDurationWin);
        switch params.StatMethod
            case 'ObsWindow'
                vsObj.ShufflingAnalysis('overwrite', params.overwriteStats, ...
                                        'N_bootstrap', params.shuffles);
            case 'bootsrapRespBase'
                vsObj.BootstrapPerNeuron('overwrite', params.overwriteStats);
            case 'maxPermuteTest'
                vsObj.StatisticsPerNeuron('overwrite', params.overwriteStats);
        end
    else
        vsObj.ResponseWindow;   % load cached result only; no recompute
    end

    % Retrieve whatever is cached (identical API regardless of compute path)
    switch params.StatMethod
        case 'ObsWindow';        stats = vsObj.ShufflingAnalysis;
        case 'bootsrapRespBase'; stats = vsObj.BootstrapPerNeuron;
        case 'maxPermuteTest';   stats = vsObj.StatisticsPerNeuron;
    end
    rw = vsObj.ResponseWindow;
end


% -------------------------------------------------------------------------
function [zS, pV, spkR, spkDiff] = extractStimVals(stats, rw, stimKey, StatMethod)
% extractStimVals  Pull z-scores, p-values, and spike metrics from a
%   stats/rw struct pair.  All outputs are N×1 column vectors.
%
%   Handles three structural cases:
%     • Speed-based subfields  (MB, MBR) – uses Speed1; prefers Speed2
%     • Moving/Static subfields (SDGm, SDGs)
%     • Flat structure          (RG, NI, NV, FFF)
%
%   Falls back to the flat-structure extractor when SDG subfields are absent
%   (e.g. when vsG is a dummy rectGridAnalysis object).

    switch stimKey

        case {'MB', 'MBR'}
            % Use the fastest available speed (Speed2 preferred over Speed1)
            sp = 'Speed1';
            if isfield(stats, 'Speed2'), sp = 'Speed2'; end

            zS      = stats.(sp).ZScoreU;
            pV      = stats.(sp).pvalsResponse;
            spkR    = max(rw.(sp).NeuronVals(:,:,4), [], 2);   % max across directions
            spkDiff = max(rw.(sp).NeuronVals(:,:,5), [], 2);   % response – baseline
            if ~isequal(StatMethod, 'ObsWindow')
                try; spkR = mean(stats.(sp).ObsResponse)'; catch; end
            end

        case 'SDGm'
            % Drifting grating – moving condition
            try
                zS      = stats.Moving.ZScoreU;
                pV      = stats.Moving.pvalsResponse;
                spkR    = max(rw.Moving.NeuronVals(:,:,4), [], 2);
                spkDiff = max(rw.Moving.NeuronVals(:,:,5), [], 2);
                if ~isequal(StatMethod, 'ObsWindow')
                    spkR = mean(stats.Moving.ObsResponse, 1)';
                end
            catch
                % Dummy vsG object (absent SDG): extract from flat struct
                [zS, pV, spkR, spkDiff] = extractFlat(stats, rw, StatMethod);
            end

        case 'SDGs'
            % Drifting grating – static condition
            try
                zS      = stats.Static.ZScoreU;
                pV      = stats.Static.pvalsResponse;
                spkR    = max(rw.Static.NeuronVals(:,:,4), [], 2);
                spkDiff = max(rw.Static.NeuronVals(:,:,5), [], 2);
                if ~isequal(StatMethod, 'ObsWindow')
                    spkR = mean(stats.Static.ObsResponse, 1)';
                end
            catch
                [zS, pV, spkR, spkDiff] = extractFlat(stats, rw, StatMethod);
            end

        otherwise   % RG, NI, NV, FFF – flat (no subfields)
            [zS, pV, spkR, spkDiff] = extractFlat(stats, rw, StatMethod);
    end

    % Guarantee column vectors regardless of how the analysis object returns data
    zS = zS(:); pV = pV(:); spkR = spkR(:); spkDiff = spkDiff(:);
end


% -------------------------------------------------------------------------
function [zS, pV, spkR, spkDiff] = extractFlat(stats, rw, StatMethod)
% extractFlat  Extract from a stats struct with no Speed/Moving/Static nesting.
    zS      = stats.ZScoreU;
    pV      = stats.pvalsResponse;
    spkR    = max(rw.NeuronVals(:,:,4), [], 2);
    spkDiff = max(rw.NeuronVals(:,:,5), [], 2);
    if ~isequal(StatMethod, 'ObsWindow')
        try; spkR = mean(stats.ObsResponse, 1)'; catch; end
    end
end


% -------------------------------------------------------------------------
function pAdj = bhFDR(pVals)
% bhFDR  Benjamini-Hochberg false-discovery rate correction.
%   Operates only on non-NaN entries; NaN values are preserved unchanged.
%
%   INPUT  pVals – N×1 (or 1×N) vector of raw p-values (may contain NaN)
%   OUTPUT pAdj  – BH-adjusted p-values, same shape as input
%
%   Algorithm:
%     1. Sort non-NaN p-values in ascending order.
%     2. Multiply each by n/rank  (BH formula).
%     3. Enforce monotonicity by a reverse cumulative minimum pass.
%     4. Cap at 1.

    pAdj      = pVals;
    validMask = ~isnan(pVals);
    p         = pVals(validMask);
    n         = numel(p);
    if n == 0, return; end

    [sortedP, sortIdx] = sort(p);
    adjP               = sortedP .* n ./ (1:n)';   % BH: p_i * n / rank_i
    adjP               = min(flipud(cummin(flipud(adjP))), 1);   % monotone & ≤ 1

    result            = zeros(n, 1);
    result(sortIdx)   = adjP;                       % restore original order
    pAdj(validMask)   = result;
end


% -------------------------------------------------------------------------
function [diffs, insers, animals] = collectPairDiffs(T, pairs, pairIdx, colName)
% collectPairDiffs  Pool per-neuron differences (stim1 – stim2) across
%   insertions for one row of the pairs matrix.
%
%   Returns column vectors diffs, insers, animals of equal length,
%   suitable for direct input to hierBoot().

    diffs = []; insers = []; animals = [];
    for ins = unique(T.insertion)'
        idx1 = T.insertion == categorical(ins) & T.stimulus == pairs{pairIdx, 1};
        idx2 = T.insertion == categorical(ins) & T.stimulus == pairs{pairIdx, 2};
        V1   = T.(colName)(idx1);
        V2   = T.(colName)(idx2);
        an   = unique(T.animal(idx1));
        diffs   = [diffs;   V1 - V2];                          %#ok<AGROW>
        insers  = [insers;  double(repmat(ins, size(V1,1), 1))]; %#ok<AGROW>
        animals = [animals; double(repmat(an,  size(V1,1), 1))]; %#ok<AGROW>
    end
end


% -------------------------------------------------------------------------
function [p1, p2, colorAnimal] = getPairScatterData(T, pairs, colName)
% getPairScatterData  Extract aligned vectors for a scatter plot.
%   p1 and p2 are the column values for the first and second stimuli;
%   colorAnimal is the animal index for colouring markers.

    p1          = T.(colName)(T.stimulus == pairs{1});
    p2          = T.(colName)(T.stimulus == pairs{2});
    colorAnimal = T.animal(T.stimulus == pairs{1});
end


% -------------------------------------------------------------------------
function [probs, ps_out] = runGroupBoot(y, x, InsIndex, AnIndex, nStims)
% runGroupBoot  Hierarchical bootstrap comparing stimulus 1 against all
%   others in the swarm data.
%
%   INPUTS
%     y        – concatenated response values (all stimuli stacked)
%     x        – stimulus index label for each element of y (1…nStims)
%     InsIndex – insertion index per neuron (same length as y for stim 1)
%     AnIndex  – animal index per neuron
%     nStims   – total number of stimuli
%
%   OUTPUTS
%     probs  – cell of Bayesian overlap probabilities (stim2…stimN vs stim1)
%     ps_out – cell of frequentist p-values (mean(BootSec >= BootFirst))

    probs  = cell(1, nStims - 1);
    ps_out = cell(1, nStims - 1);

    FirstStim = y(x == 1);
    validF    = ~isnan(FirstStim);
    BootFirst = hierBoot(FirstStim(validF), 10000, InsIndex(validF), AnIndex(validF));

    for i = 2:nStims
        sec   = y(x == i);
        sec(isnan(sec)) = 0;                        % NaN → 0 (absent = no response)
        validMask = ~isinf(sec);                    % [SUGG-2] isinf replaces ==-inf
        sec   = sec(validMask);
        BootSec  = hierBoot(sec, 10000, InsIndex(validMask), AnIndex(validMask));
        probs{i-1}  = get_direct_prob(BootFirst, BootSec);
        ps_out{i-1} = mean(BootSec >= BootFirst);
    end
end


% -------------------------------------------------------------------------
function applyPaperAxes(ax)
% applyPaperAxes  Apply standard axis formatting for publication.
%   Centralising this in one function means font size / family changes
%   propagate everywhere with a single edit.
    ax.YAxis.FontSize = 8;  ax.YAxis.FontName = 'helvetica';
    ax.XAxis.FontSize = 8;  ax.XAxis.FontName = 'helvetica';
end