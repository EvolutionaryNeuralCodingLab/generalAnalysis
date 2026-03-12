
%% Run/load bombcell and confusion matrices

%
exp = [49:54,64:97];%
%tiledlayout(numel(exp),1)
for ex =  exp%GoodRecordingsPV%allGoodRec %GoodRecordings%GoodRecordingsPV%GoodRecordingsPV%selecN{1}(1,:) %1:size(data,1)
    %%%%%%%%%%%% Load data and data paremeters
    %1. Load NP class
    ex=53
    NP = loadNPclassFromTable(ex);
    vs = linearlyMovingBallAnalysis(NP,Session=1);
    KSversion =4;

    [qMetric,unitType]=NP.getBombCell(NP.recordingDir+"\kilosort4",0,KSversion,1);

    %convertPhySorting2tIc(obj,pathToPhyResults,tStart,BombCelled)

    %
    % goodUnits = unitType == 1;
    % muaUnits = unitType == 2;
    % noiseUnits = unitType == 0;
    % nonSomaticUnits = unitType == 3;

    % Concordance analysis
    % bc    load_manual_classifications(vs.spikeSortingFolder)
    % pMC = vs.dataObj.convertPhySorting2tIc(vs.spikeSortingFolder,0,0,1);
    % pBC = vs.dataObj.convertPhySorting2tIc(vs.spikeSortingFolder,0,1,1);

    bombcell_table = readtable([vs.spikeSortingFolder filesep 'cluster_bc_unitType.tsv'], 'FileType', 'text', 'Delimiter', '\t');
    manual_table = readtable([vs.spikeSortingFolder filesep 'cluster_info.tsv'],'FileType','delimitedtext');

    manual_table = manual_table(:,{'cluster_id','KSLabel','group'});
    sum(strcmp(pKS.label, 'good'))

   
    % Load and prepare data
    % Assume:
    % bombcell_table: Nx2 table, columns: [id, bc_label]  ("GOOD","MUA","NON-SOMA","NOISE")
    % manual_table:   Mx3 table, columns: [id, KS_label, group]  ("good","mua","noise")

    % Rename columns for clarity (adjust if yours differ)
    bombcell_table.Properties.VariableNames = {'id', 'bc_label'};
    manual_table.Properties.VariableNames   = {'id', 'KS_label', 'group'};

    % Remove NON-SOMA from bombcell
    bc = bombcell_table(~strcmp(bombcell_table.bc_label, 'NON-SOMA'), :);

    % Match IDs — keep only IDs present in both tables
    [~, ia, ib] = intersect(bc.id, manual_table.id);
    bc_matched  = bc(ia, :);
    man_matched = manual_table(ib, :);

    % Harmonize labels to lowercase for comparison
    bc_labels  = lower(bc_matched.bc_label);    % "good","mua","noise"
    ks_labels  = lower(man_matched.KS_label);   % "good","mua","noise"
    man_labels = lower(man_matched.group);      % "good","mua","noise"

    %%Define category order
    cats = {'good', 'mua', 'noise'};

    bc_cat  = categorical(bc_labels,  cats);
    ks_cat  = categorical(ks_labels,  cats);
    man_cat = categorical(man_labels, cats);

    % --- Confusion Matrix 1: Manual curation vs BombCell ---
    % figure('Position', [100, 100, 700, 600]);
    %
    % tiledlayout(3,2)
    % nexttile
    % cm1 = confusionchart(man_cat, bc_cat, ...
    %     'Title', sprintf('%s-Manual curation vs BombCell',NP.recordingName),...
    %     'XLabel', 'BombCell', ...
    %     'YLabel', 'Manual Curation', ...
    %     'RowSummary', 'row-normalized', ...
    %     'ColumnSummary', 'column-normalized');
    %
    % cm1.FontSize = 9;
    %
    % % Give the chart more room inside the figure
    % %cm1.Position = [10, 10, 680, 580];

    % --- Confusion Matrix 2: KS label vs BombCell ---
    fig = figure('Position', [100, 100, 700, 600]);
    %tl = nexttile;
    cm2 = confusionchart(ks_cat, bc_cat, ...       
        'XLabel', 'BombCell', ...
        'YLabel', 'KS Label', ...
        'RowSummary', 'row-normalized', ...
        'ColumnSummary', 'column-normalized');
    cm2.FontSize = 9;
    title(sprintf('%KS Label vs BombCell',NP.recordingName));



    % %% --- Confusion Matrix 3: KS label vs Manual curation ---
    % figure;
    % cm3 = confusionchart(ks_cat, man_cat, ...
    %     'Title', printf('KS Label vs Manual Curation',NP.recordingName), ...
    %     'XLabel', 'Manual Curation', ...
    %     'YLabel', 'KS Label', ...
    %     'RowSummary', 'row-normalized', ...
    %     'ColumnSummary', 'column-normalized');

    % --- Print mismatch summary ---
    % fprintf('\n=== Manual vs BombCell ===\n')
    % mismatch_man_bc = man_cat ~= bc_cat;
    % fprintf('Total mismatches: %d / %d (%.1f%%)\n', ...
    %     sum(mismatch_man_bc), numel(mismatch_man_bc), ...
    %     100*mean(mismatch_man_bc));

    fprintf('\n=== KS Label vs BombCell ===\n')
    mismatch_ks_bc = ks_cat ~= bc_cat;
    fprintf('Total mismatches: %d / %d (%.1f%%)\n', ...
        sum(mismatch_ks_bc), numel(mismatch_ks_bc), ...
        100*mean(mismatch_ks_bc));

    vs.printFig(fig,sprintf('%KS Label vs BombCell',NP.recordingName),PaperFig =1)

    close

    % fprintf('\n=== KS Label vs Manual Curation ===\n')
    % mismatch_ks_man = ks_cat ~= man_cat;
    % fprintf('Total mismatches: %d / %d (%.1f%%)\n', ...
    %     sum(mismatch_ks_man), numel(mismatch_ks_man), ...
    %     100*mean(mismatch_ks_man));

    imec = Neuropixel.ImecDataset(NP.recordingDir);
    ks = Neuropixel.KilosortDataset(vs.spikeSortingFolder,'imecDataset', imec);
    ks.load();

end
%I want to compare bombcell unit classification with manual classification in phy.



%% Plot raw waveforms of specific units:

% 1. Add to path: https://github.com/cortex-lab/spikes
%                 https://github.com/kwikteam/npy-matlab  (dependency)


ksDir = vs.spikeSortingFolder;
sp = loadKSdir(ksDir);   % loads all KS output into a struct

% Get waveforms
gwfparams.dataDir    = ksDir;
gwfparams.fileName   = NP.recordingDir;
gwfparams.dataType   = 'int16';
gwfparams.nCh        = 385;
gwfparams.wfWin      = [-40 41];   % samples around spike
gwfparams.nWf        = 100;        % waveforms per unit
gwfparams.spikeTimes = sp.st;      % spike times
gwfparams.spikeClusters = sp.clu;  % cluster IDs

wf = getWaveForms(gwfparams);      % wf.waveForms: [units x waveforms x channels x samples]

% Plot mean waveform for unit 1, best channel
figure;
plot(squeeze(mean(wf.waveFormsMean(1,:,:), 2)));

%% Check low amp waveforms 10 neurons per experiment

PVexps = [49:54,64:97];
idx = randi(length(PVexps), 1, 4);
selected = PVexps(idx);



for i = selected
    NP = loadNPclassFromTable(53);
    vs = linearlyMovingBallAnalysis(NP,Session=1);

    p = vs.dataObj.convertPhySorting2tIc(vs.spikeSortingFolder,0,1,1);
    phy_IDg = p.phy_ID(string(p.label') == 'good');


    plotRawWaveforms(vs, [47:50], showCorr=true, corrWin=50, corrBin=0.5)

end
