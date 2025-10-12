function fig = PlotZScoreComparison(expList, Stims2Comp,params)

arguments
    expList  (1,:) double  %%Number of experiment from excel list
    Stims2Comp cell %% Comparison order {'MB','RG','MBR'} would select neurons responsive to moving ball and 
    % compare this neurons responses to other stimuli. 
    params.threshold = 0.005;
    params.diffResp = false;
    params.overwrite = false;
    params.StimsPresent = {'MB'}; %assumes that at list moving ball is present
    params.StimsNotPresent = {};
end

% Compare z-scores and p-values between moving ball and rect grid analyses

animal = 0;
animalVector = cell(1,numel(expList));
zScoresMB = cell(1,numel(expList));
zScoresRG = cell(1,numel(expList));
spKrMB = cell(1,numel(expList));
spKrRG = cell(1,numel(expList));
diffSpkMB = cell(1,numel(expList));
diffSpkRG = cell(1,numel(expList));

zScoresSDGm = cell(1,numel(expList));
zScoresMBR = cell(1,numel(expList));
zScoresFFF = cell(1,numel(expList));
spKrMBR = cell(1,numel(expList));
spKrFFF = cell(1,numel(expList));
spKrSDGm = cell(1,numel(expList));
diffSpkMBR = cell(1,numel(expList));
diffSpkFFF = cell(1,numel(expList));
diffSpkSDGm = cell(1,numel(expList));

zScoresNI = cell(1,numel(expList));
% zScoresNV = cell(1,numel(expList));
spKrNI = cell(1,numel(expList));
spKrNV = cell(1,numel(expList));
diffSpkNI = cell(1,numel(expList));
diffSpkNV = cell(1,numel(expList));

j = 1;
AnimalI = "";

NP = loadNPclassFromTable(expList(1)); %73 81
vs = linearlyMovingBallAnalysis(NP);

%%% Asumes all experiments were analyzed using the same window
vs.ResponseWindow;
MBvs = vs.ResponseWindow;
%%%

nameOfFile = sprintf('\\Ex_%d-%d_Combined_Neural_responses_%s_filtered.mat',expList(1),expList(end),Stims2Comp{1});
p = extractBefore(vs.getAnalysisFileName,'lizards');
p = [p 'lizards'];

if ~exist([p '\Combined_lizard_analysis'],'dir')
    cd(p)    
    mkdir Combined_lizard_analysis
end
saveDir = [p '\Combined_lizard_analysis'];

if exist([saveDir nameOfFile],'file') == 2 && ~params.overwrite 

    S = load([saveDir nameOfFile]);

    expList2 = S.expList;

    if isequal(expList2,expList)

        spKrMB = S.spKrMB;
        spKrRG = S.spKrRG;
        diffSpkMB = S.diffSpkMB;
        diffSpkRG = S.diffSpkRG;
        animalVector = S.animalVector;
        zScoresMB = S.zScoresMB;
        zScoresRG = S.zScoresRG;

        if isfield(S,'spKrMBR')
            spKrMBR = S.spKrMBR;
            diffSpkMBR = S.diffSpkMBR;
            zScoresMBR = S.zScoresMBR;
        end

        if isfield(S,'spKrSDGm')
            spKrSDGm = S.spKrSDGm; %Moving
            diffSpkSDGm = S.diffSpkSDGm;
            zScoresSDGm = S.zScoresSDGm;

            spKrSDGs = S.spKrSDGs; %Static
            diffSpkSDGs = S.diffSpkSDGs;
            zScoresSDGs = S.zScoresSDGs;
        end

        if isfield(S,'spKrFFF')
            spKrFFF = S.spKrFFF;
            diffSpkFFF = S.diffSpkFFF;
            zScoresFFF = S.zScoresFFF;
        end

        if isfield(S,'spKrNI')
            spKrNI = S.spKrNI;
            diffSpkNI = S.diffSpkNI;
            zScoresNI = S.zScoresNI;
        end

        if isfield(S,'spKrNV')
            spKrNV = S.spKrNV;
            diffSpkNV = S.diffSpkNV;
            zScoresNV = S.zScoresNV;
        end

        forloop = false;
    else
        forloop = true;
    end
else
    forloop = true;
end

if forloop
    for ex = expList

        NP = loadNPclassFromTable(ex); %73 81
        vs = linearlyMovingBallAnalysis(NP);
        vsR = rectGridAnalysis(NP);
        try
            vsBr = linearlyMovingBarAnalysis(NP);
            params.StimsPresent{3} = 'MBR';
        catch
            params.StimsPresent{3} = '';
            fprintf('Moving Bar stimulus not found.\n')
            vsBr = linearlyMovingBallAnalysis(NP); %use rectGrid here to avoid puting lots of ifs.
        end
        try
            vsG = StaticDriftingGratings(NP);
            params.StimsPresent{4} = 'SDG';
        catch
           params.StimsPresent{4} = '';
            fprintf('Gratings stimulus not found.\n')
            vsG = rectGridAnalysis(NP); %use rectGrid here to avoid puting lots of ifs.
        end
        try
            vsNI = NaturalImage(NP);
            params.StimsPresent{5} = 'NI';
        catch
            params.StimsPresent{5} = '';
            fprintf('Natural images stimulus not found.\n')
            vsNI = rectGridAnalysis(NP); %use rectGrid here to avoid puting lots of ifs.
        end
        try
            vsNV = NaturalVideo(NP);
            params.StimsPresent{6} = 'NV';
        catch
            params.StimsPresent{6} = '';
            fprintf('Natural video stimulus not found.\n')
            vsNV = rectGridAnalysis(NP); %use rectGrid here to avoid puting lots of ifs.
        end

        try
            vsFFF = FullFieldFlashAnalysis(NP);
            params.StimsPresent{7} = 'FFF';
        catch
            params.StimsPresent{7} = '';
            fprintf('FFF stimulus not found.\n')
            vsFFF = rectGridAnalysis(NP); %use moving ball here to avoid puting lots of ifs.
        end
        

        %%Load pvals and zscore from rect grid and moving ball
        vs.ResponseWindow;
        vsR.ResponseWindow;
        vsBr.ResponseWindow;
        vsG.ResponseWindow;
        vsNI.ResponseWindow;
        vsNV.ResponseWindow;
        vsFFF.ResponseWindow;

        try
            statsMB = vs.ShufflingAnalysis;
        catch
            vs.ShufflingAnalysis;
            statsMB = vs.ShufflingAnalysis;
        end

        try
            statsRG = vsR.ShufflingAnalysis;
        catch
            vsR.ShufflingAnalysis;
            statsRG =  vsR.ShufflingAnalysis;
        end

        try
            statsMBR = vsBr.ShufflingAnalysis;
        catch
            vsBr.ShufflingAnalysis;
            statsMBR =  vsBr.ShufflingAnalysis;
        end

        try
            statsSDG = vsG.ShufflingAnalysis;
        catch
            vsG.ShufflingAnalysis;
            statsSDG =  vsG.ShufflingAnalysis;
        end

        try
            statsFFF = vsFFF.ShufflingAnalysis;
        catch
            vsFFF.ShufflingAnalysis;
            statsFFF =  vsFFF.ShufflingAnalysis;
        end


        try
            statsNI = vsNI.ShufflingAnalysis;
        catch
            vsNI.ShufflingAnalysis;
            statsNI =  vsNI.ShufflingAnalysis;
        end

        try
            statsNV = vsNV.ShufflingAnalysis;
        catch
            vsNV.ShufflingAnalysis;
            statsNV =  vsNV.ShufflingAnalysis;
        end

        rwRG = vsR.ResponseWindow;
        rwMB = vs.ResponseWindow;
        rwMBR = vsBr.ResponseWindow;
        rwFFF = vsFFF.ResponseWindow;
        rwSDG = vsG.ResponseWindow;
        rwNI = vsNI.ResponseWindow;
        rwNV = vsNV.ResponseWindow;

        %Load stats of Moving Ball, select fastest speed if there are several
        zScores_MB = statsMB.Speed1.ZScoreU;
        pValuesMB = statsMB.Speed1.pvalsResponse;
        spkR_MB = max(rwMB.Speed1.NeuronVals(:,:,1),[],2);
        spkDiff_MB = max(rwMB.Speed1.NeuronVals(:,:,4),[],2);

        if isfield(statsMB, 'Speed2') %If
            zScores_MB = statsMB.Speed2.ZScoreU;
            pValuesMB = statsMB.Speed2.pvalsResponse;
            spkR_MB = max(rwMB.Speed2.NeuronVals(:,:,1),[],2);
            spkDiff_MB = max(rwMB.Speed2.NeuronVals(:,:,4),[],2);
        end

        %Load stats of Rect Grid.
        zScores_RG = statsRG.ZScoreU;
        pValuesRG = statsRG.pvalsResponse;
        spkR_RG = max(rwRG.NeuronVals(:,:,1),[],2);
        spkDiff_RG = max(rwRG.NeuronVals(:,:,4),[],2);

        %Load stats of Moving bar.
        zScores_MBR = statsMBR.Speed1.ZScoreU;
        pValuesMBR = statsMBR.Speed1.pvalsResponse;
        spkR_MBR = max(rwMBR.Speed1.NeuronVals(:,:,1),[],2);
        spkDiff_MBR = max(rwMBR.Speed1.NeuronVals(:,:,4),[],2);

        %Load stats of FFF
        zScores_FFF = statsFFF.ZScoreU;
        pValuesFFF = statsFFF.pvalsResponse;
        spkR_FFF = max(rwFFF.NeuronVals(:,:,1),[],2);
        spkDiff_FFF = max(rwFFF.NeuronVals(:,:,4),[],2);

        %Load stats of SDG moving
        zScores_SDGm = statsSDG.ZScoreU;
        pValuesSDGm = statsSDG.pvalsResponse;
        spkR_SDGm = max(rwSDG.NeuronVals(:,:,1),[],2);
        spkDiff_SDGm = max(rwSDG.NeuronVals(:,:,4),[],2);

        %Load stats of SDG static
        zScores_SDGs = statsSDG.ZScoreU;
        pValuesSDGs = statsSDG.pvalsResponse;
        spkR_SDGs = max(rwSDG.NeuronVals(:,:,1),[],2);
        spkDiff_SDGs = max(rwSDG.NeuronVals(:,:,4),[],2);


        %Load stats of SDG moving
        zScores_NI = statsNI.ZScoreU;
        pValuesNI = statsNI.pvalsResponse;
        spkR_NI = max(rwNI.NeuronVals(:,:,1),[],2);
        spkDiff_NI = max(rwNI.NeuronVals(:,:,4),[],2);

        %Load stats of SDG static
        zScores_NV = statsNV.ZScoreU;
        pValuesNV = statsNV.pvalsResponse;
        spkR_NV = max(rwNV.NeuronVals(:,:,1),[],2);
        spkDiff_NV = max(rwNV.NeuronVals(:,:,4),[],2);

        pvals = {'pValuesMB','pValuesRG','pValuesMBR','pValuesFFF','pValuesSDGm','pValuesSDGs','pValuesNI','pValuesNV'...
            ;pValuesMB,pValuesRG,pValuesMBR,pValuesFFF,pValuesSDGm,pValuesSDGs,pValuesNI,pValuesNV};

        [row, col] = find(cellfun(@(x) ischar(x) && endsWith(x, Stims2Comp{1}), pvals));

        pvalsStimSelected = pvals{2,col};


        zScores_MB = zScores_MB(pvalsStimSelected<=params.threshold);
        spkR_MB =  spkR_MB(pvalsStimSelected<=params.threshold);
        spkDiff_MB =  spkDiff_MB(pvalsStimSelected<=params.threshold);

        zScores_RG = zScores_RG(pvalsStimSelected<=params.threshold);
        spkR_RG = spkR_RG(pvalsStimSelected<=params.threshold);
        spkDiff_RG = spkDiff_RG(pvalsStimSelected<=params.threshold);

        zScores_MBR = zScores_MBR(pvalsStimSelected<=params.threshold);
        spkR_MBR = spkR_MBR(pvalsStimSelected<=params.threshold);
        spkDiff_MBR = spkDiff_MBR(pvalsStimSelected<=params.threshold);

        zScores_SDGm = zScores_SDGm(pvalsStimSelected<=params.threshold);
        spkR_SDGm = spkR_SDGm(pvalsStimSelected<=params.threshold);
        spkDiff_SDGm = spkDiff_SDGm(pvalsStimSelected<=params.threshold);

        zScores_SDGs = zScores_SDGs(pvalsStimSelected<=params.threshold);
        spkR_SDGs = spkR_SDGs(pvalsStimSelected<=params.threshold);
        spkDiff_SDGs = spkDiff_SDGs(pvalsStimSelected<=params.threshold);

        zScores_FFF = zScores_FFF(pvalsStimSelected<=params.threshold);
        spkR_FFF = spkR_FFF(pvalsStimSelected<=params.threshold);
        spkDiff_FFF = spkDiff_FFF(pvalsStimSelected<=params.threshold);

        zScores_NI= zScores_NI(pvalsStimSelected<=params.threshold);
        spkR_NI = spkR_NI(pvalsStimSelected<=params.threshold);
        spkDiff_NI = spkDiff_NI(pvalsStimSelected<=params.threshold);

        zScores_NV = zScores_NV(pvalsStimSelected<=params.threshold);
        spkR_NV = spkR_NV(pvalsStimSelected<=params.threshold);
        spkDiff_NV = spkDiff_NV(pvalsStimSelected<=params.threshold);

        %Check animal name
        Animal = string(regexp( vs.getAnalysisFileName, 'PV\d+', 'match', 'once'));

        if Animal ~= AnimalI %wont work if you start with the first animal (noisy animal)
            animal = animal+1;
            AnimalNames{animal} = Animal;
            AnimalI = Animal;
        end

        animalVector{j} = repmat(animal,[1, numel(zScores_MB)]);
        zScoresMB{j} = zScores_MB;
        zScoresRG{j} = zScores_RG;

        spKrMB{j} = spkR_MB';
        spKrRG{j} = spkR_RG';
        diffSpkMB{j} = spkDiff_MB;
        diffSpkRG{j} = spkDiff_RG;

        zScoresFFF{j} = zScores_FFF;
        spKrFFF{j} = spkR_FFF';
        diffSpkFFF{j} = spkDiff_FFF;

        zScoresMBR{j} = zScores_MBR;
        spKrMBR{j} = spkR_MBR';
        diffSpkMBR{j} = spkDiff_MBR;

        zScoresSDGm{j} = zScores_SDGm;
        spKrSDGm{j} = spkR_SDGm';
        diffSpkSDGm{j} = spkDiff_SDGm;

        zScoresSDGs{j} = zScores_SDGs;
        spKrSDGs{j} = spkR_SDGs';
        diffSpkSDGs{j} = spkDiff_SDGs;

        zScoresNI{j} = zScores_NI;
        spKrNI{j} = spkR_NI';
        diffSpkNI{j} = spkDiff_NI;

        zScoresNV{j} = zScores_NV;
        spKrNV{j} = spkR_NV';
        diffSpkNV{j} = spkDiff_NV;
        % Create a figure for comparison

        j = j +1;

    end

    S.spKrMB = spKrMB;
    S.spKrRG = spKrRG;
    S.diffSpkMB = diffSpkMB;
    S.diffSpkRG = diffSpkRG;
    S.zScoresMB =zScoresMB;
    S.zScoresRG = zScoresRG;
    S.expList = expList;
    S.animalVector = animalVector;
    S.params = params;

    S.spKrMBR = spKrMBR;
    S.spKrFFF = spKrFFF;
    S.diffSpkMBR = diffSpkMBR;
    S.diffSpkFFF = diffSpkFFF;
    S.zScoresMBR =zScoresMBR;
    S.zScoresFFF = zScoresFFF;

    S.spKrSDGm = spKrSDGm;
    S.spKrSDGs = spKrSDGs;
    S.diffSpkSDGm = diffSpkSDGm;
    S.diffSpkSDGs = diffSpkSDGs;
    S.zScoresSDGm =zScoresSDGm;
    S.zScoresSDGs = zScoresSDGs;

    S.spKrNI = spKrNI;
    S.spKrNV = spKrNV;
    S.diffSpkNI= diffSpkNI;
    S.diffSpkNV = diffSpkNV;
    S.zScoresNI =zScoresNI;
    S.zScoresNV = zScoresNV;

    S.expList = expList;
    S.params = params;

    save([saveDir nameOfFile],'-struct', 'S');
end

fig=figure;

tiledlayout(2,2,"TileSpacing","compact");

fn = fieldnames(S);
StimZS = cell(numel(Stims2Comp),1);
stimRSP = cell(numel(Stims2Comp),1);

x = [];

for i = 1:numel(Stims2Comp)

    ending = Stims2Comp{i};
    pattern = ['^zS.*' ending '$'];
    matches = fn(~cellfun('isempty', regexp(fn, pattern)));

    StimZS{i} = cell2mat(S.(matches{1}))';

    if  ~params.diffResp
        pattern = ['^spKr.*' ending '$'];
    else
        pattern = ['^diffSpk.*' ending '$'];
    end
    
    matches = fn(~cellfun('isempty', regexp(fn, pattern)));
    stimRSP{i} = cell2mat(S.(matches{1}))';

    x = [x;ones(size(cell2mat(S.(matches{1}))'))*i];

end
% 
% MB = cell2mat(zScoresMB)';
% RG = cell2mat(zScoresRG)';
AnIndex = cell2mat(animalVector)';
colormapUsed = parula(max(AnIndex)).*0.6;

% eMB = cell2mat(EntropMB);
% eRG = cell2mat(EntropRG);

%y =log10([MB; RG; SDGm; SDGs; MB-RG; MB-SDGm; MB-SDGs]);

%y =([MB; RG; MB-RG]);
y = cell2mat(StimZS);

y = cell2mat(stimRSP);

% x =[];
% % Combine data
% for i =1:length(y)/length(MB)
%     x = [x;ones(size(MB))*i];
% end

% Combine all color indices
allColorIndices = repmat(AnIndex,numel(Stims2Comp),1);

% ---- Swarmchart (Larger Left Subplot) ----
nexttile % Takes most of the space
swarmchart(x, y, 5, [colormapUsed(allColorIndices,:)], 'filled','MarkerFaceAlpha',0.7); % Marker size 50
xticks(1:7);
xticklabels(Stims2Comp);
ylabel('SpikeRate');
set(fig,'Color','w')
%set(gca, 'YScale', 'log')
yline([0],'LineWidth',2)
ylim([-5 40])


nexttile
scatter(StimZS{1},StimZS{2},10,AnIndex,"filled","MarkerFaceAlpha",0.5)
colormap(colormapUsed)
hold on
axis equal
lims = [0 40];
plot(lims, lims, 'k--', 'LineWidth', 1.5)
ylim(lims)
xlim(lims)
xlabel('MB')
ylabel('SB')

%%% Spike rate comparison


if ~params.diffResp
    MB = cell2mat(spKrMB)'.*(1000/MBvs.params.durationWindow);
    RG = cell2mat(spKrRG)'.*(1000/MBvs.params.durationWindow);
else
    MB = cell2mat(diffSpkMB').*(1000/MBvs.params.durationWindow);
    RG = cell2mat(diffSpkRG')*(1000/MBvs.params.durationWindow);
end
AnIndex = cell2mat(animalVector)';

y =([MB; RG; MB-RG]);


x =[];
% Combine data
for i =1:length(y)/length(MB)
    x = [x;ones(size(MB))*i];
end

% Combine all color indices
allColorIndices = repmat(AnIndex,length(y)/length(MB),1);

% ---- Swarmchart (Larger Left Subplot) ----
nexttile% Takes most of the space
swarmchart(x, y, 5, [colormapUsed(allColorIndices,:)], 'filled','MarkerFaceAlpha',0.7); % Marker size 50
xticks(1:7);
xticklabels({'MB', 'SB', 'MB-SB'});
ylabel('Spikes/sec');
if params.diffResp
    ylabel('Resp-Base spkR');
end
%set(gca, 'YScale', 'log')
%ylim([-0.01,0.1])
yline([0],'LineWidth',2)

%lims = [-5 50];
nexttile
scatter(MB,RG,10,AnIndex,"filled","MarkerFaceAlpha",0.5)
colormap(colormapUsed)
xlims = xlim;
ylims = ylim;
lims = [0 max([xlims ylims])];
hold on
plot(lims, lims, 'k--', 'LineWidth', 1.5)
axis equal
ylim(lims)
xlim(lims)
xlabel('MB')
ylabel('SB')
fig.Position = [454   147   776   831];

end