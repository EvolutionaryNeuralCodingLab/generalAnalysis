function fig = PlotZScoreComparison(expList, Stims2Comp,params)

arguments
    expList  (1,:) double  %%Number of experiment from excel list
    Stims2Comp cell %% Comparison order {'MB','RG','MBR'} would select neurons responsive to moving ball and 
    % compare this neurons responses to other stimuli. 
    params.threshold = 0.05;
    params.diffResp = false;
    params.overwrite = false;
    params.StimsPresent = {'MB','RG'}; %assumes that at least moving ball is present
    params.StimsNotPresent = {};
    params.StimsToCompare = {}; %Select 2 stims to compare scatter plots (default: 1st and 2nd stim are compared from the Stims2Comp cell array)
    params.overwriteResponse = false;
    params.RespDurationWin = 100; %same as default
    params.shuffles = 2000; %same as default
    params.ignoreNonSignif = false; %when comparing first stim, ignore neurons non responsive to other stim
    params.EachStimSignif = false; %resposnive neurons for each stim are selected (default: responsive neurons of first stime are selected)
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

        forloop = false;
    else
        forloop = true;
    end
else
    forloop = true;
end

if forloop
    for ex = expList

        fprintf('Processing recording: %s .\n',NP.recordingName)
        NP = loadNPclassFromTable(ex); %73 81
        vs = linearlyMovingBallAnalysis(NP);
        vsR = rectGridAnalysis(NP);

        try
            vsBr = linearlyMovingBarAnalysis(NP);
            params.StimsPresent{3} = 'MBR';
            if isempty(vsBr.VST)
                error('Moving Bar stimulus not found.\n')
            end
        catch
            params.StimsPresent{3} = '';
            fprintf('Moving Bar stimulus not found.\n')
            vsBr = linearlyMovingBallAnalysis(NP); %use rectGrid here to avoid puting lots of ifs.
        end
        try
            vsG = StaticDriftingGratingAnalysis(NP);
            params.StimsPresent{4} = 'SDG';
             if isempty(vsG.VST)
                error('Gratings stimulus not found.\n')
            end
        catch
           params.StimsPresent{4} = '';
            fprintf('Gratings stimulus not found.\n')
            vsG = rectGridAnalysis(NP); %use rectGrid here to avoid puting lots of ifs.
        end
        try
            vsNI = imageAnalysis(NP);
            params.StimsPresent{5} = 'NI';
             if isempty(vsNI.VST)
                error('Gratings stimulus not found.\n')
            end
        catch
            params.StimsPresent{5} = '';
            fprintf('Natural images stimulus not found.\n')
            vsNI = rectGridAnalysis(NP); %use rectGrid here to avoid puting lots of ifs.
        end
        try
            vsNV = movieAnalysis(NP);
            params.StimsPresent{6} = 'NV';
            if isempty(vsNV.VST)
                error('Gratings stimulus not found.\n')
            end
        catch
            params.StimsPresent{6} = '';
            fprintf('Natural video stimulus not found.\n')
            vsNV = rectGridAnalysis(NP); %use rectGrid here to avoid puting lots of ifs.
        end

        try
            vsFFF = fullFieldFlashAnalysis(NP);
            params.StimsPresent{7} = 'FFF';
            if isempty(vsFFF.VST)
                error('Gratings stimulus not found.\n')
            end
        catch
            params.StimsPresent{7} = '';
            fprintf('FFF stimulus not found.\n')
            vsFFF = rectGridAnalysis(NP); %use moving ball here to avoid puting lots of ifs.
        end
        

        %%Load pvals and zscore from rect grid and moving ball
        vs.ResponseWindow('overwrite',params.overwriteResponse,'durationWindow',params.RespDurationWin);
        if isequal(params.StimsPresent{2},'')
            vsR.ResponseWindow;
        else
             vsR.ResponseWindow('overwrite',params.overwriteResponse,'durationWindow',params.RespDurationWin);
        end
        if isequal(params.StimsPresent{3},'')
            vsBr.ResponseWindow;
        else
             vsBr.ResponseWindow('overwrite',params.overwriteResponse,'durationWindow',params.RespDurationWin);      
        end
        if isequal(params.StimsPresent{4},'')
            vsG.ResponseWindow;
        else
            vsG.ResponseWindow('overwrite',params.overwriteResponse,'durationWindow',params.RespDurationWin);
        end
        if isequal(params.StimsPresent{5},'')
            vsNI.ResponseWindow;
        else
            vsNI.ResponseWindow('overwrite',params.overwriteResponse,'durationWindow',params.RespDurationWin);
        end
        if isequal(params.StimsPresent{6},'')
            vsNV.ResponseWindow;
        else
            vsBr.ResponseWindow('overwrite',params.overwriteResponse,'durationWindow',params.RespDurationWin);
        end
        if isequal(params.StimsPresent{7},'')
            vsFFF.ResponseWindow;
        else
            vsFFF.ResponseWindow('overwrite',params.overwriteResponse,'durationWindow',params.RespDurationWin);
        end

        vs.ShufflingAnalysis('overwrite',params.overwriteResponse,"N_bootstrap", params.shuffles);
        statsMB = vs.ShufflingAnalysis;

        vsR.ShufflingAnalysis('overwrite',params.overwriteResponse,"N_bootstrap", params.shuffles);
        statsRG =  vsR.ShufflingAnalysis;


        vsBr.ShufflingAnalysis('overwrite',params.overwriteResponse,"N_bootstrap", params.shuffles);
        statsMBR =  vsBr.ShufflingAnalysis;

        vsG.ShufflingAnalysis('overwrite',params.overwriteResponse,"N_bootstrap", params.shuffles);
        statsSDG =  vsG.ShufflingAnalysis;

        vsFFF.ShufflingAnalysis('overwrite',params.overwriteResponse,"N_bootstrap", params.shuffles);
        statsFFF =  vsFFF.ShufflingAnalysis;

        vsNI.ShufflingAnalysis('overwrite',params.overwriteResponse,"N_bootstrap", params.shuffles);
        statsNI =  vsNI.ShufflingAnalysis;

        vsNV.ShufflingAnalysis('overwrite',params.overwriteResponse,"N_bootstrap", params.shuffles);
        statsNV =  vsNV.ShufflingAnalysis;


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
        spkR_MB = max(rwMB.Speed1.NeuronVals(:,:,4),[],2);
        spkDiff_MB = max(rwMB.Speed1.NeuronVals(:,:,5),[],2);

        if isfield(statsMB, 'Speed2') %If
            zScores_MB = statsMB.Speed2.ZScoreU;
            pValuesMB = statsMB.Speed2.pvalsResponse;
            spkR_MB = max(rwMB.Speed2.NeuronVals(:,:,4),[],2);
            spkDiff_MB = max(rwMB.Speed2.NeuronVals(:,:,5),[],2);
        end

        totalU{j} = numel(zScores_MB);
        %Load stats of Rect Grid.
        zScores_RG = statsRG.ZScoreU;
        pValuesRG = statsRG.pvalsResponse;
        spkR_RG = max(rwRG.NeuronVals(:,:,4),[],2);
        spkDiff_RG = max(rwRG.NeuronVals(:,:,5),[],2);

        %Load stats of Moving bar.
        zScores_MBR = statsMBR.Speed1.ZScoreU;
        pValuesMBR = statsMBR.Speed1.pvalsResponse;
        spkR_MBR = max(rwMBR.Speed1.NeuronVals(:,:,4),[],2);
        spkDiff_MBR = max(rwMBR.Speed1.NeuronVals(:,:,5),[],2);

        %Load stats of FFF
        zScores_FFF = statsFFF.ZScoreU;
        pValuesFFF = statsFFF.pvalsResponse;
        spkR_FFF = max(rwFFF.NeuronVals(:,:,4),[],2);
        spkDiff_FFF = max(rwFFF.NeuronVals(:,:,5),[],2);

        %Load stats of SDG moving

        if isequal(params.StimsPresent{4},'')

            zScores_SDGm = statsSDG.ZScoreU;
            pValuesSDGm = statsSDG.pvalsResponse;
            spkR_SDGm = max(rwSDG.NeuronVals(:,:,4),[],2);
            spkDiff_SDGm = max(rwSDG.NeuronVals(:,:,5),[],2);

            %Load stats of SDG static
            zScores_SDGs = statsSDG.ZScoreU;
            pValuesSDGs = statsSDG.pvalsResponse;
            spkR_SDGs = max(rwSDG.NeuronVals(:,:,4),[],2);
            spkDiff_SDGs = max(rwSDG.NeuronVals(:,:,5),[],2);

        else
            zScores_SDGm = statsSDG.Moving.ZScoreU;
            pValuesSDGm = statsSDG.Moving.pvalsResponse;
            spkR_SDGm = max(rwSDG.Moving.NeuronVals(:,:,4),[],2);
            spkDiff_SDGm = max(rwSDG.Moving.NeuronVals(:,:,5),[],2);

            %Load stats of SDG static
            zScores_SDGs = statsSDG.Static.ZScoreU;
            pValuesSDGs = statsSDG.Static.pvalsResponse;
            spkR_SDGs = max(rwSDG.Static.NeuronVals(:,:,4),[],2);
            spkDiff_SDGs = max(rwSDG.Static.NeuronVals(:,:,5),[],2);
        end

        %Load stats of Natural images
        zScores_NI = statsNI.ZScoreU;
        pValuesNI = statsNI.pvalsResponse;
        spkR_NI = max(rwNI.NeuronVals(:,:,4),[],2);
        spkDiff_NI = max(rwNI.NeuronVals(:,:,5),[],2);

        %Load stats of video
        zScores_NV = statsNV.ZScoreU;
        pValuesNV = statsNV.pvalsResponse;
        spkR_NV = max(rwNV.NeuronVals(:,:,4),[],2);
        spkDiff_NV = max(rwNV.NeuronVals(:,:,5),[],2);

        if params.ignoreNonSignif

            zScores_NV(pValuesNV>params.threshold) = -1000;
            zScores_NI(pValuesNI>params.threshold) = -1000;
            zScores_SDGs(pValuesSDGs>params.threshold) = -1000;
            zScores_SDGm(pValuesSDGm>params.threshold) = -1000;
            zScores_FFF(pValuesFFF>params.threshold) = -1000;
            zScores_MBR(pValuesMBR>params.threshold) = -1000;
            zScores_RG(pValuesRG>params.threshold) = -1000;
            zScores_MB(pValuesMB>params.threshold) = -1000;

        end

        pvals = {'pValuesMB','pValuesRG','pValuesMBR','pValuesFFF','pValuesSDGm','pValuesSDGs','pValuesNI','pValuesNV'...
            ;pValuesMB,pValuesRG,pValuesMBR,pValuesFFF,pValuesSDGm,pValuesSDGs,pValuesNI,pValuesNV};

        [row, col] = find(cellfun(@(x) ischar(x) && endsWith(x, Stims2Comp{1}), pvals));

        pvalsStimSelected = pvals{2,col};


         %%% Responsive to one specific stimuli: 

        zScores_MBs = zScores_MB(pvalsStimSelected<=params.threshold);
        spkR_MBs =  spkR_MB(pvalsStimSelected<=params.threshold);
        spkDiff_MBs =  spkDiff_MB(pvalsStimSelected<=params.threshold);
        pvals_MB = pValuesMB(pvalsStimSelected<=params.threshold);
        sumNeurMB = numel(zScores_MB);

        zScores_RGs = zScores_RG(pvalsStimSelected<=params.threshold);
        spkR_RGs = spkR_RG(pvalsStimSelected<=params.threshold);
        spkDiff_RGs = spkDiff_RG(pvalsStimSelected<=params.threshold);
        pvals_RG = pValuesRG(pvalsStimSelected<=params.threshold);
        sumNeurRG = numel(zScores_RG);

        if isequal(params.StimsPresent{2},'') %Asign - inf values if stim is not present in recording
            zScores_RGs = zScores_RG-inf;
            spkR_RGs = zScores_RG-inf;
            spkDiff_RGs = zScores_RG-inf;
            pvals_RG = zScores_RG-inf;
            sumNeurRG = 0;
        end

        zScores_MBRs = zScores_MBR(pvalsStimSelected<=params.threshold);
        spkR_MBRs = spkR_MBR(pvalsStimSelected<=params.threshold);
        spkDiff_MBRs = spkDiff_MBR(pvalsStimSelected<=params.threshold);
        pvals_MBR = pValuesMBR(pvalsStimSelected<=params.threshold);
        sumNeurMBR = numel(zScores_MBR);

        if isequal(params.StimsPresent{3},'') %Asign - inf values if stim is not present in recording
            zScores_MBRs = zScores_MBRs-inf;
            spkR_MBRs = zScores_MBRs-inf;
            spkDiff_MBRs = zScores_MBRs-inf;
            pvals_MBR = zScores_MBRs-inf;
            sumNeurMBR = 0;
        end

        zScores_SDGms = zScores_SDGm(pvalsStimSelected<=params.threshold);
        spkR_SDGms = spkR_SDGm(pvalsStimSelected<=params.threshold);
        spkDiff_SDGms = spkDiff_SDGm(pvalsStimSelected<=params.threshold);
        pvals_SDGm = pValuesSDGm(pvalsStimSelected<=params.threshold);

        zScores_SDGss = zScores_SDGs(pvalsStimSelected<=params.threshold);
        spkR_SDGss = spkR_SDGs(pvalsStimSelected<=params.threshold);
        spkDiff_SDGss = spkDiff_SDGs(pvalsStimSelected<=params.threshold);
        pvals_SDGs = pValuesSDGs(pvalsStimSelected<=params.threshold);
        sumNeurSDG = numel(zScores_SDGs);

        if isequal(params.StimsPresent{4},'') %Asign - inf values if stim is not present in recording
            zScores_SDGss = zScores_SDGss-inf;
            spkR_SDGss = spkR_SDGss-inf;
            spkDiff_SDGss = spkDiff_SDGss-inf;
            pvals_SDGs = pvals_SDGs-inf;

            zScores_SDGms = zScores_SDGms-inf;
            spkR_SDGms = spkR_SDGms-inf;
            spkDiff_SDGms = spkDiff_SDGms-inf;
            pvals_SDGm = pvals_SDGm-inf;
            sumNeurSDG = 0;
        end

        zScores_FFFs = zScores_FFF(pvalsStimSelected<=params.threshold);
        spkR_FFFs = spkR_FFF(pvalsStimSelected<=params.threshold);
        spkDiff_FFFs = spkDiff_FFF(pvalsStimSelected<=params.threshold);
        pvals_FFF = pValuesFFF(pvalsStimSelected<=params.threshold);
        sumNeurFFF = numel(zScores_FFF);

        if isequal(params.StimsPresent{7},'') %Asign - inf values if stim is not present in recording
            zScores_FFFs = zScores_FFFs-inf;
            spkR_FFFs = spkR_FFFs-inf;
            spkDiff_FFFs = spkDiff_FFFs-inf;
            pvals_FFF = pvals_FFF-inf;
            sumNeurFFF = 0;
        end

        zScores_NIs= zScores_NI(pvalsStimSelected<=params.threshold);
        spkR_NIs = spkR_NI(pvalsStimSelected<=params.threshold);
        spkDiff_NIs = spkDiff_NI(pvalsStimSelected<=params.threshold);
        pvals_NI = pValuesNI(pvalsStimSelected<=params.threshold);
        sumNeurNI = numel(zScores_NI);


        if isequal(params.StimsPresent{5},'') %Asign - inf values if stim is not present in recording
            zScores_NIs = zScores_NIs-inf;
            spkR_NIs = spkR_NIs-inf;
            spkDiff_NIs = spkDiff_NIs-inf;
            pvals_NI = pvals_NI-inf;
            sumNeurNI = 0;
        end

        zScores_NVs = zScores_NV(pvalsStimSelected<=params.threshold);
        spkR_NVs = spkR_NV(pvalsStimSelected<=params.threshold);
        spkDiff_NVs = spkDiff_NV(pvalsStimSelected<=params.threshold);
        pvals_NV = pValuesNV(pvalsStimSelected<=params.threshold);
        sumNeurNV = numel(zScores_NV);
        
        if isequal(params.StimsPresent{6},'') %Asign - inf values if stim is not present in recording
            zScores_NVs = zScores_NV-inf;
            spkR_NVs = spkR_NVs-inf;
            spkDiff_NVs = spkDiff_NVs-inf;
            pvals_NV = pvals_NV-inf;
            sumNeurNV = 0;
        end
        


        %%% Responsive in general:
        respIndexes = [];
        zScores_MBg = zScores_MB(pValuesMB<=params.threshold);
        spkR_MBg =  spkR_MB(pValuesMB<=params.threshold);
        spkDiff_MBg =  spkDiff_MB(pValuesMB<=params.threshold);
        respIndexes = [respIndexes find(pValuesMB<=params.threshold)];

        zScores_RGg = zScores_RG(pValuesRG<=params.threshold);
        spkR_RGg = spkR_RG(pValuesRG<=params.threshold);
        spkDiff_RGg = spkDiff_RG(pValuesRG<=params.threshold);
        respIndexes = [respIndexes find(pValuesRG<=params.threshold)];

        zScores_MBRg = zScores_MBR(pValuesMBR<=params.threshold);
        spkR_MBRg = spkR_MBR(pValuesMBR<=params.threshold);
        spkDiff_MBRg = spkDiff_MBR(pValuesMBR<=params.threshold);
        respIndexes = [respIndexes find(pValuesMBR<=params.threshold)];

        zScores_SDGmg = zScores_SDGm(pValuesSDGm<=params.threshold);
        spkR_SDGmg = spkR_SDGm(pValuesSDGm<=params.threshold);
        spkDiff_SDGmg = spkDiff_SDGm(pValuesSDGm<=params.threshold);
        respIndexes = [respIndexes find(pValuesSDGm<=params.threshold)];

        zScores_SDGsg = zScores_SDGs(pValuesSDGs<=params.threshold);
        spkR_SDGsg = spkR_SDGs(pValuesSDGs<=params.threshold);
        spkDiff_SDGsg = spkDiff_SDGs(pValuesSDGs<=params.threshold);
        respIndexes = [respIndexes find(pValuesSDGs<=params.threshold)];

        zScores_FFFg = zScores_FFF(pValuesFFF<=params.threshold);
        spkR_FFFg = spkR_FFF(pValuesFFF<=params.threshold);
        spkDiff_FFFg = spkDiff_FFF(pValuesFFF<=params.threshold);
        respIndexes = [respIndexes find(pValuesFFF<=params.threshold)];

        zScores_NIg = zScores_NI(pValuesNI<=params.threshold);
        spkR_NIg = spkR_NI(pValuesNI<=params.threshold);
        spkDiff_NIg = spkDiff_NI(pValuesNI<=params.threshold);
        respIndexes = [respIndexes find(pValuesNI<=params.threshold)];

        zScores_NVg = zScores_NV(pValuesNV<=params.threshold);
        spkR_NVg = spkR_NV(pValuesNV<=params.threshold);
        spkDiff_NVg = spkDiff_NV(pValuesNV<=params.threshold);
        respIndexes = [respIndexes find(pValuesNV<=params.threshold)];

        responsiveNeuronsj = unique(respIndexes);

        %Check animal name
        Animal = string(regexp( vs.getAnalysisFileName, 'PV\d+', 'match', 'once'));

        if isequal(Animal,"")
              Animal = string(regexp( vs.getAnalysisFileName, 'SA\d+', 'match', 'once'));
        end

        if Animal ~= AnimalI %wont work if you start with the first animal (noisy animal)
            animal = animal+1;
            AnimalNames{animal} = Animal;
            AnimalI = Animal;
        end

        animalVector{j} = repmat(animal,[1, numel(zScores_MBs)]);
        zScoresMB{j} = zScores_MBs;
        zScoresRG{j} = zScores_RGs;
        pvalsRG{j} = pvals_RG;
        sumNeurRGt{j} = sumNeurRG;
        pvalsMB{j} = pvals_MB;
        sumNeurMBt{j} = sumNeurMB;

        spKrMB{j} = spkR_MBs';
        spKrRG{j} = spkR_RGs';
        diffSpkMB{j} = spkDiff_MBs;
        diffSpkRG{j} = spkDiff_RGs;

        zScoresFFF{j} = zScores_FFFs;
        spKrFFF{j} = spkR_FFFs';
        diffSpkFFF{j} = spkDiff_FFFs;
        pvalsFFF{j} = pvals_FFF;
        sumNeurFFFt{j} = sumNeurFFFt;

        zScoresMBR{j} = zScores_MBRs;
        spKrMBR{j} = spkR_MBRs';
        diffSpkMBR{j} = spkDiff_MBRs;
        pvalsMBR{j} = pvals_MBR;
        sumNeurMBRt{j} = sumNeurMBR;

        zScoresSDGm{j} = zScores_SDGms;
        spKrSDGm{j} = spkR_SDGms';
        diffSpkSDGm{j} = spkDiff_SDGms;
        pvalsSDGm{j} = pvals_SDGm;

        zScoresSDGs{j} = zScores_SDGss;
        spKrSDGs{j} = spkR_SDGss';
        diffSpkSDGs{j} = spkDiff_SDGss;
        pvalsSDGs{j} = pvals_SDGs;
        sumNeurSDGt{j} = sumNeurSDG;

        zScoresNI{j} = zScores_NIs;
        spKrNI{j} = spkR_NIs';
        diffSpkNI{j} = spkDiff_NIs;
        pvalsNI{j} = pvals_NI;
        sumNeurNIt{j} = sumNeurNI;

        zScoresNV{j} = zScores_NVs;
        spKrNV{j} = spkR_NVs';
        diffSpkNV{j} = spkDiff_NVs;
        pvalsNV{j} = pvals_NV;
        sumNeurNVt{j} = sumNeurNV;
        
        %%% Responsive in general:

        zScoresMBg{j} = zScores_MBg;
        spkRMBg{j} =  spkR_MBg;
        spkDiffMBg{j} =  spkDiff_MBg;

        zScoresRGg{j} = zScores_RGg;
        spkRRGg{j} = spkR_RGg;
        spkDiffRGg{j} = spkDiff_RGg;

        zScoresMBRg{j} = zScores_MBRg;
        spkRMBRg{j} = spkR_MBRg;
        spkDiffMBRg{j} = spkDiff_MBRg;

        zScoresSDGmg{j} = zScores_SDGmg;
        spkRSDGmg{j} = spkR_SDGmg;
        spkDiffSDGmg{j} = spkDiff_SDGmg;

        zScoresSDGsg{j} = zScores_SDGsg;
        spkRSDGsg{j} = spkR_SDGsg;
        spkDiffSDGsg{j} = spkDiff_SDGsg;

        zScoresFFFg{j} = zScores_FFFg;
        spkRFFFg{j} = spkR_FFFg;
        spkDiffFFFg{j} = spkDiff_FFFg;

        zScoresNIg{j} = zScores_NIg;
        spkRNIg{j} = spkR_NIg;
        spkDiffNIg{j} = spkDiff_NIg;

        zScoresNVg{j} = zScores_NVg;
        spkRNVg{j} = spkR_NVg;
        spkDiffNVg{j} = spkDiff_NVg;

        responsiveNeurons{j} = responsiveNeuronsj; 

        % Create a figure for comparison

        j = j +1;

        fprintf('Finished recording: %s .\n',NP.recordingName)

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
    S.pvalsMB = pvalsMB;
    S.pvalsRG = pvalsRG;
    S.sumNeurRG = sumNeurRGt;
    S.sumNeurMB = sumNeurMBt;

    S.spKrMBR = spKrMBR;
    S.spKrFFF = spKrFFF;
    S.diffSpkMBR = diffSpkMBR;
    S.diffSpkFFF = diffSpkFFF;
    S.zScoresMBR =zScoresMBR;
    S.zScoresFFF = zScoresFFF;
    S.pvalsFFF = pvalsFFF;
    S.pvalsMBR = pvalsMBR;
    S.sumNeurMBR = sumNeurMBRt;

    S.spKrSDGm = spKrSDGm;
    S.spKrSDGs = spKrSDGs;
    S.diffSpkSDGm = diffSpkSDGm;
    S.diffSpkSDGs = diffSpkSDGs;
    S.zScoresSDGm =zScoresSDGm;
    S.zScoresSDGs = zScoresSDGs;
    S.pvalsSDGm = pvalsSDGm;
    S.pvalsSDGs = pvalsSDGs;
    S.sumNeurSDG = sumNeurSDGt;

    S.spKrNI = spKrNI;
    S.spKrNV = spKrNV;
    S.diffSpkNI= diffSpkNI;
    S.diffSpkNV = diffSpkNV;
    S.zScoresNI =zScoresNI;
    S.zScoresNV = zScoresNV;
    S.pvalsNI = pvalsNI;
    S.pvalsNV = pvalsNV;
    S.sumNeurNV = sumNeurNVt;
    S.sumNeurNI = sumNeurNIt;

    S.zScoresMBg=zScoresMBg;
    S.spkRMBg=spkRMBg;
    S.spkDiffMBg=spkDiffMBg;

    S.zScoresRGg=zScoresRGg;
    S.spkRRGg=spkRRGg;
    S.spkDiffRGg=spkDiffRGg;

    S.zScoresMBRg=zScoresMBRg;
    S.spkRMBRg=spkRMBRg;
    S.spkDiffMBRg=spkDiffMBRg;

    S.zScoresSDGmg=zScoresSDGmg;
    S.spkRSDGmg=spkRSDGmg;
    S.spkDiffSDGmg=spkDiffSDGmg;

    S.zScoresSDGsg=zScoresSDGsg;
    S.spkRSDGsg=spkRSDGsg;
    S.spkDiffSDGsg=spkDiffSDGsg;

    S.zScoresFFFg=zScoresFFFg;
    S.spkRFFFg=spkRFFFg;
    S.spkDiffFFFg=spkDiffFFFg;

    S.zScoresNIg=zScoresNIg;
    S.spkRNIg=spkRNIg;
    S.spkDiffNIg=spkDiffNIg;

    S.zScoresNVg=zScoresNVg;
    S.spkRNVg=spkRNVg;
    S.spkDiffNVg=spkDiffNVg;

    S.expList = expList;
    S.totalUnits = totalU;
    S.params = params;
    S.responsiveNeurons = responsiveNeurons;

    save([saveDir nameOfFile],'-struct', 'S');
end

fig=figure;

tiledlayout(2,2,"TileSpacing","compact");

fn = fieldnames(S);

StimZS = cell(numel(Stims2Comp),1);
stimRSP = cell(numel(Stims2Comp),1);
stimPvals = cell(numel(Stims2Comp),1);

x = [];
endingOpts = {'','g'};
if params.EachStimSignif
    ending2 = endingOpts{2};
else
    ending2 = endingOpts{1};
end

Stims2Comp2 = {};
for i = 1:numel(Stims2Comp)
    if strcmp(Stims2Comp{i}, 'SDG')
        % Replace 'SDG' with two new elements
        Stims2Comp2 = [Stims2Comp2, {'SDGs', 'SDGm'}];
    else
        Stims2Comp2 = [Stims2Comp2, Stims2Comp(i)];
    end
end


StimZS = cell(numel(Stims2Comp2),1);
stimRSP = cell(numel(Stims2Comp2),1);
stimPvals = cell(numel(Stims2Comp2),1);

for i = 1:numel(Stims2Comp2)

    ending = Stims2Comp2{i};
    pattern = ['^zS.*' ending ending2 '$'];
    matches = fn(~cellfun('isempty', regexp(fn, pattern)));

    StimZS{i} = cell2mat(S.(matches{1}))';

    if  ~params.diffResp
        pattern = ['^spKr.*' ending ending2 '$'];
    else
        pattern = ['^diffSpk.*' ending ending2 '$'];
    end
    
    matches = fn(~cellfun('isempty', regexp(fn, pattern)));

    if params.EachStimSignif
        matches = fn(~cellfun('isempty', regexp(fn, pattern,'ignorecase')));
        C = S.(matches{1});
        C = cellfun(@(x) x', C, 'UniformOutput', false);
        stimRSP{i} = cell2mat(C);
    else
        try
            stimRSP{i} =  cell2mat(S.(matches{1})');
        catch
            try
                stimRSP{i} =  cell2mat(S.(matches{1}));
            catch
                stimRSP{i} =  cell2mat(S.(matches{1})');
            end
        end
    end

    pattern = ['^pvals.*' ending '$'];

    matches = fn(~cellfun('isempty', regexp(fn, pattern)));
    stimPvals{i} = cell2mat(S.(matches{1}))';

    x = [x;ones(size(StimZS{i}))*i];

end
% 
% MB = cell2mat(zScoresMB)';
% RG = cell2mat(zScoresRG)';
AnIndex = cell2mat(S.animalVector)';
colormapUsed = parula(max(AnIndex)).*0.6;

% eMB = cell2mat(EntropMB);
% eRG = cell2mat(EntropRG);

%y =log10([MB; RG; SDGm; SDGs; MB-RG; MB-SDGm; MB-SDGs]);

%y =([MB; RG; MB-RG]);
% StimZS{numel(StimZS)+1} = [StimZS{1}-StimZS{2}];
y = cell2mat(StimZS);

% x = [x;ones(size(cell2mat(S.(matches{1}))'))*i+1];

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
xticks(1:8);
xticklabels(Stims2Comp2);
ylabel('Z-score');
set(fig,'Color','w')
%set(gca, 'YScale', 'log')
yline([0],'LineWidth',2)
ylim([-5 40])


nexttile
%stims to compare
% boxplot(y2,'Labels',Stims2Comp)


if isempty(params.StimsToCompare)
    ind1 = 1;
    ind2 = 2;
else
    
    ind1 = find(strcmp(Stims2Comp, params.StimsToCompare{1}));
    ind2 = find(strcmp(Stims2Comp, params.StimsToCompare{2}));
   
end

ValsToCompare = {StimZS{ind1},StimZS{ind2}};

scatter(ValsToCompare{1},ValsToCompare{2},10,AnIndex,"filled","MarkerFaceAlpha",0.5)
colormap(colormapUsed)
hold on
axis equal

lims =[min(y) max(y)]; 
plot(lims, lims, 'k--', 'LineWidth', 1.5)
lims = [-5 40];
ylim(lims)
xlim(lims)
xlabel(Stims2Comp(ind1))
ylabel(Stims2Comp(ind2))

%%% Spike rate comparison


y = cell2mat(stimRSP');

% ---- Swarmchart (Larger Left Subplot) ----
nexttile % Takes most of the space
swarmchart(x, y, 5, [colormapUsed(allColorIndices,:)], 'filled','MarkerFaceAlpha',0.7); % Marker size 50
xticks(1:8);
xticklabels(Stims2Comp2);
ylabel('Spike Rate');
set(fig,'Color','w')
%set(gca, 'YScale', 'log')
%yline([0],'LineWidth',2)


ValsToCompare = {stimRSP{ind1},stimRSP{ind2}};
nexttile
scatter(ValsToCompare{1},ValsToCompare{2},10,AnIndex,"filled","MarkerFaceAlpha",0.5)
colormap(colormapUsed)
hold on
axis equal
lims = [0 max(xlim)];
plot(lims, lims, 'k--', 'LineWidth', 1.5)
ylim(lims)
xlim(lims)
xlabel(Stims2Comp(ind1))
ylabel(Stims2Comp(ind2))

%%Number of responsive units:


for i = 1:numel(Stims2Comp2)

    ending = [Stims2Comp2{i} 'g'];
    pattern = ['^zS.*' ending '$'];
    matches = fn(~cellfun('isempty', regexp(fn, pattern)));

    RespNeurCount{i} = numel(cell2mat(S.(matches{1})))';

end

fig2=figure;

AllRespNeurs = numel(cell2mat(S.responsiveNeurons));

BarVals = [cell2mat(RespNeurCount)];

bar([AllRespNeurs BarVals],'FaceColor','k','FaceAlpha',0.7)

xticklabels({'Total resp. neurons',Stims2Comp2{:}})

ylabel('Responsive neuron count')


end