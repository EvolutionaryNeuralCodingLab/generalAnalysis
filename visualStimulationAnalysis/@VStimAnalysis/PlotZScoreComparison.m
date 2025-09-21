function fig = PlotZScoreComparison(expList, Stims2Comp,params)

arguments
    expList  (1,:) double  %%Number of experiment from excel list
    Stims2Comp cell %% Comparison order
    params.threshold = 0.005;
end

% Compare z-scores and p-values between moving ball and rect grid analyses

animal = 0;
animalVector = cell(1,numel(expList));
zScoresMBall = cell(1,numel(expList));
zScoresRGall = cell(1,numel(expList));
j = 1;
AnimalI = "";

for ex = expList

    NP = loadNPclassFromTable(ex); %73 81
    vs = linearlyMovingBallAnalysis(NP);

    %%Load pvals and zscore from rect grid and moving ball
    vsR = rectGridAnalysis(NP);
    vs.ResponseWindow;
    vsR.ResponseWindow;
    
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

    rwRG = vsR.ResponseWindow;
    rwMB = vs.ResponseWindow;

    %Load stats of Moving Ball, select fastest speed if there are several
    zScoresMB = statsMB.Speed1.ZScoreU;
    pValuesMB = statsMB.Speed1.pvalsResponse;
    if isfield(statsMB, 'Speed2') 
        zScoresMB = statsMB.Speed2.ZScoreU;
        pValuesMB = statsMB.Speed2.pvalsResponse;
    end

    %Load stats of Rect Grid. 
    zScoresRG = statsRG.ZScoreU;
    pValuesRG = statsRG.pvalsResponse;

    if strcmp(Stims2Comp{1}, 'MB') %%Check if MB goes first, if so select responsive units to MB
        zScoresMB = zScoresMB(pValuesMB<=params.threshold);
        zScoresRG = zScoresRG(pValuesMB<=params.threshold);
    else
        zScoresMB = zScoresMB(pValuesRG<=params.threshold);
        zScoresRG = zScoresRG(pValuesRG<=params.threshold);
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

    % Create a figure for comparison

    j = j +1;
    
end

fig=figure;

MB = cell2mat(zScoresMBall)';
RG = cell2mat(zScoresRGall)';
AnIndex = cell2mat(animalVector)';
colormapUsed = parula(animal).*0.8;

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
subplot(1, 15, 1:13); % Takes most of the space
swarmchart(x, y, 5, [colormapUsed(allColorIndices,:)], 'filled','MarkerFaceAlpha',0.7); % Marker size 50
xticks(1:7);
xticklabels({'MB', 'SB', 'MB-SB'});
ylabel('Z-score');
set(fig,'Color','w')
%set(gca, 'YScale', 'log')
yline([0],'LineWidth',2)
ylim([-5 100])
fig.Position = [680   680   359   298];

lims = [-5 50];
figure;
scatter(MB,RG)
ylim(lims)
xlim(lims)
hold on
plot(lims, lims, 'k--', 'LineWidth', 1.5)
axis equal
end