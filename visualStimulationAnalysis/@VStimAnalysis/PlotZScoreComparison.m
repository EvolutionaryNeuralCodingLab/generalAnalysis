function fig = PlotZScoreComparison(expList, Stims2Comp,params)

arguments
    expList  (1,:) double  %%Number of experiment from excel list
    Stims2Comp cell %% Comparison order {'MB','RG','MBR'} would select neurons responsive to moving ball and 
    % compare this neurons responses to other stimuli. 
    params.threshold = 0.005;
    params.diffResp = false;
    params.overwrite = false;
    % params.MBRpresent = true;
    % params.SDGpresent = false;
    % params.FFFpresent = false;
end

% Compare z-scores and p-values between moving ball and rect grid analyses

animal = 0;
animalVector = cell(1,numel(expList));
zScoresMBall = cell(1,numel(expList));
zScoresRGall = cell(1,numel(expList));
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
zScoresNV = cell(1,numel(expList));
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


p = extractBefore(vs.getAnalysisFileName,'lizards');
p = [p 'lizards'];

if ~exist([p '\Combined_lizard_analysis'],'dir')
    cd(p)    
    mkdir Combined_lizard_analysis
end
saveDir = [p '\Combined_lizard_analysis'];

if exist(saveDir,'dir') && ~params.overwrite 

    S = load([saveDir '\Combined_Neural_responses.mat']);

    expList2 = S.expList;

    if isequal(expList2,expList)

        spKrMB = S.spKrMB;
        spKrRG = S.spKrRG;
        diffSpkMB = S.diffSpkMB;
        diffSpkRG = S.spKrdiffSpkRGMB;
        animalVector = S.animalVector;
        zScoresMBall = S.zScoresMBall;
        zScoresRGall = S.zScoresRGall;

        if isfield(S,'spKrMBR')
            spKrMBR = S.spKrMBR;
            diffSpkMBR = S.diffSpkMBR;
            zScoresMBar = S.zScoresMBar;
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
        catch
            MBRpresent = false;
            fprintf('Moving Bar stimulus not found.\n')
            vsBr = linearlyMovingBallAnalysis(NP); %use moving ball here to avoid puting lots of ifs.
        end
        try
            vsG = StaticDriftingGratings(NP);
        catch
            vsGpresent = false;
            fprintf('Gratings stimulus not found.\n')
            vsG = linearlyMovingBallAnalysis(NP); %use moving ball here to avoid puting lots of ifs.
        end
        try
            vsNI = NaturalImage(NP);
        catch
            vsGpresent = false;
            fprintf('Natural images stimulus not found.\n')
            vsG = linearlyMovingBallAnalysis(NP); %use moving ball here to avoid puting lots of ifs.
        end
        try
            vsNV = NaturalVideo(NP);
        catch
            vsGpresent = false;
            fprintf('Natural video stimulus not found.\n')
            vsG = linearlyMovingBallAnalysis(NP); %use moving ball here to avoid puting lots of ifs.
        end

        %%Load pvals and zscore from rect grid and moving ball
        vs.ResponseWindow;
        vsR.ResponseWindow;
        sBr.ResponseWindow;
        sG.ResponseWindow;
        sNI.ResponseWindow;
        sNV.ResponseWindow;

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
        rwSDG = sG.ResponseWindow;
        rwNI = NI.ResponseWindow;
        rwNV = NV.ResponseWindow;

        %Load stats of Moving Ball, select fastest speed if there are several
        zScoresMB = statsMB.Speed1.ZScoreU;
        pValuesMB = statsMB.Speed1.pvalsResponse;
        spkR_MB = max(rwMB.Speed1.NeuronVals(:,:,1),[],2);
        spkDiff_MB = max(rwMB.Speed1.NeuronVals(:,:,4),[],2);

        if isfield(statsMB, 'Speed2') %If
            zScoresMB = statsMB.Speed2.ZScoreU;
            pValuesMB = statsMB.Speed2.pvalsResponse;
            spkR_MB = max(rwMB.Speed2.NeuronVals(:,:,1),[],2);
            spkDiff_MB = max(rwMB.Speed2.NeuronVals(:,:,4),[],2);
        end

        %Load stats of Rect Grid.
        zScoresRG = statsRG.ZScoreU;
        pValuesRG = statsRG.pvalsResponse;
        spkR_RG = max(rwRG.NeuronVals(:,:,1),[],2);
        spkDiff_RG = max(rwRG.NeuronVals(:,:,4),[],2);


        %Load stats of Rect Grid.
        zScoresMBR = statsMBR.ZScoreU;
        pValuesMBR = statsMBR.pvalsResponse;
        spkR_MBR = max(rwMBR.NeuronVals(:,:,1),[],2);
        spkDiff_MBR = max(rwMBR.NeuronVals(:,:,4),[],2);

        if strcmp(Stims2Comp{1}, 'MB') %%Check if MB goes first, if so select responsive units to MB
            zScoresMB = zScoresMB(pValuesMB<=params.threshold);
            spkR_MB =  spkR_MB(pValuesMB<=params.threshold);
            spkDiff_MB =  spkDiff_MB(pValuesMB<=params.threshold);

            zScoresRG = zScoresRG(pValuesMB<=params.threshold);
            spkR_RG = spkR_RG(pValuesMB<=params.threshold);
            spkDiff_RG = spkDiff_RG(pValuesMB<=params.threshold);

            zScoresMBR = zScoresMBR(pValuesMB<=params.threshold);
            spkR_MBR = spkR_MBR(pValuesMB<=params.threshold);
            spkDiff_MBR = spkDiff_MBR(pValuesMB<=params.threshold);

        elseif strcmp(Stims2Comp{1}, 'RG')
            zScoresMB = zScoresMB(pValuesRG<=params.threshold);
            spkR_MB =  spkR_MB(pValuesRG<=params.threshold);
            spkDiff_MB =  spkDiff_MB(pValuesRG<=params.threshold);

            zScoresRG = zScoresRG(pValuesRG<=params.threshold);
            spkR_RG = spkR_RG(pValuesRG<=params.threshold);
            spkDiff_RG = spkDiff_RG(pValuesRG<=params.threshold);

            zScoresMBR = zScoresMBR(pValuesRG<=params.threshold);
            spkR_MBR = spkR_MBR(pValuesRG<=params.threshold);
            spkDiff_MBR = spkDiff_MBR(pValuesRG<=params.threshold);

        elseif strcmp(Stims2Comp{1}, 'MBR')
            zScoresMB = zScoresMB(pValuesMBR<=params.threshold);
            spkR_MB =  spkR_MB(pValuesMBR<=params.threshold);
            spkDiff_MB =  spkDiff_MB(pValuesMBR<=params.threshold);

            zScoresRG = zScoresRG(pValuesMBR<=params.threshold);
            spkR_RG = spkR_RG(pValuesMBR<=params.threshold);
            spkDiff_RG = spkDiff_RG(pValuesMBR<=params.threshold);

            zScoresMBR = zScoresMBR(pValuesMBR<=params.threshold);
            spkR_MBR = spkR_MBR(pValuesMBR<=params.threshold);
            spkDiff_MBR = spkDiff_MBR(pValuesMBR<=params.threshold);
        end

        %Check animal name
        Animal = string(regexp( vs.getAnalysisFileName, 'PV\d+', 'match', 'once'));

        if Animal ~= AnimalI %wont work if you start with the first animal (noisy animal)
            animal = animal+1;
            AnimalNames{animal} = Animal;
            AnimalI = Animal;
        end

        animalVector{j} = repmat(animal,[1, numel(zScoresMB)]);
        zScoresMBall{j} = zScoresMB;
        zScoresRGall{j} = zScoresRG;

        spKrMB{j} = spkR_MB';
        spKrRG{j} = spkR_RG';
        diffSpkMB{j} = spkDiff_MB;
        diffSpkRG{j} = spkDiff_RG;

        zScoresMBR{j} = zScoresMBR;
        spKrMBR{j} = spkR_MBR';
        diffSpkMBR{j} = spkDiff_MBR;


        % Create a figure for comparison

        j = j +1;

    end

    S.spKrMB = spKrMB;
    S.spKrRG = spKrRG;
    S.diffSpkMB = diffSpkMB;
    S.diffSpkRG = diffSpkRG;
    S.zScoresMBall =zScoresMBall;
    S.zScoresRGall = zScoresRGall;
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

    save([saveDir '\Combined_Neural_responses.mat'],'-struct', 'S');
end

fig=figure;

tiledlayout(2,2,"TileSpacing","compact");

MB = cell2mat(zScoresMBall)';
RG = cell2mat(zScoresRGall)';
AnIndex = cell2mat(animalVector)';
colormapUsed = parula(max(AnIndex)).*0.8;

% eMB = cell2mat(EntropMB);
% eRG = cell2mat(EntropRG);

%y =log10([MB; RG; SDGm; SDGs; MB-RG; MB-SDGm; MB-SDGs]);

y =([MB; RG; MB-RG]);


x =[];
% Combine data
for i =1:length(y)/length(MB)
    x = [x;ones(size(MB))*i];
end

% Combine all color indices
allColorIndices = repmat(AnIndex,length(y)/length(MB),1);

% ---- Swarmchart (Larger Left Subplot) ----
nexttile % Takes most of the space
swarmchart(x, y, 5, [colormapUsed(allColorIndices,:)], 'filled','MarkerFaceAlpha',0.7); % Marker size 50
xticks(1:7);
xticklabels({'MB', 'SB', 'MB-SB'});
ylabel('Z-score');
set(fig,'Color','w')
%set(gca, 'YScale', 'log')
yline([0],'LineWidth',2)
ylim([-5 30])



nexttile
scatter(MB,RG,10,AnIndex,"filled","MarkerFaceAlpha",0.5)
colormap(colormapUsed)
hold on
axis equal
lims = [-5 50];
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