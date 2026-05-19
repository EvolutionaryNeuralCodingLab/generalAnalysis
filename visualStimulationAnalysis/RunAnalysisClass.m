cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

%%
%% Rect Grid
for ex = [69] %84:91
    NP = loadNPclassFromTable(ex); %73 81
    vsRe = rectGridAnalysis(NP);
    % vsRe.getSessionTime("overwrite",true);
    % %vsRe.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    % vsRe.getDiodeTriggers('overwrite',true);
    % vsRe.getSyncedDiodeTriggers("overwrite",true);
    % % vsRe.plotSpatialTuningSpikes;
    % % vsRe.plotSpatialTuningLFP;
     %vsRe.ResponseWindow('overwrite',true)
    % results = vsRe.ShufflingAnalysis('overwrite',true);
    %vsRe.plotRaster(MergeNtrials=1,overwrite=true,AllResponsiveNeurons = true, selectedLum=[],oneTrial = true,PaperFig = true) %43
    vsRe.plotRaster(MergeNtrials=1,overwrite=true,AllResponsiveNeurons=false,exNeurons=21, selectedLum=255,oneTrial = true,PaperFig = true) %43
    vsRe.CalculateReceptiveFields('overwrite',true,AllResponsiveNeurons=false)
    [colorbarLimsRG] = vsRe.PlotReceptiveFields(exNeurons=21,allStimParamsCombined=false,PaperFig=true,overwrite=true);
    %result = vsRe.BootstrapPerNeuron('overwrite',true);
    %result = vsRe.StatisticsPerNeuron('overwrite',true);

end
% vsRe.CalculateReceptiveFields
%vsRe.PlotReceptiveFields("exNeurons",18)

%% Moving ball

for ex =[88]%97  74:84  (Neurons, 96_74, )
    NP = loadNPclassFromTable(ex); %73 81
    vs = linearlyMovingBallAnalysis(NP,Multiplesizes=true);
    % vs.getSessionTime("overwrite",true);
    % vs.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    % % %vs.plotDiodeTriggers
    % vs.getSyncedDiodeTriggers("overwrite",true);
    % % % %vs.plotSpatialTuningSpikes;
    % r = vs.ResponseWindow('overwrite',true);
    % % % results = vs.ShufflingAnalysis('overwrite',true);
    % % % % vs.plotRaster('AllSomaticNeurons',true,'overwrite',true,'MergeNtrials',3)
    vs.plotRaster('AllResponsiveNeurons',true,'overwrite',true,'MergeNtrials',1,'bin',50,'GaussianLength',30,'MaxVal_1', false, oneLuminosity = "white", OneDirection="left", ...
        sortingOrder=["size","direction","luminosity","offset","speed"])
    %vs.plotRaster('AllSomaticNeurons',true,'overwrite',true,'MergeNtrials',3,PaperFig=true)
    % % %vs.plotRaster('exNeuronsPhyID',288,'overwrite',true,'MergeNtrials',3,'PaperFig',true)
    % % % % %vs.plotCorrSpikePattern
    % vs.plotRaster('exNeurons',82,'overwrite',true,'OneDirection','up','OneLuminosity','white','MergeNtrials',1,'PaperFig',false)
    % %vs.plotRaster('AllSomaticNeurons',true,'overwrite',true,'MergeNtrials',3,MaxVal_1=false,PaperFig=false)
    % vs.CalculateReceptiveFields('overwrite',true,testConvolution=false,AllResponsiveNeurons=false);   
    % colorbarLims=vs.PlotReceptiveFields('exNeurons',21,'overwrite',true,'OneDirection','up','OneLuminosity','white','PaperFig',true);
    % %result = vs.BootstrapPerNeuron('overwrite',true);%('overwrite',true);
    % % pvals0_6Filter =result.Speed2.pvalsResponse';
    % % compare = [pvals,pvalsNoFilt,pvals0_6Filter];
    % result = vs.StatisticsPerNeuron('overwrite',true);
    % result = vs.StatisticsPerNeuronPerCategory('compareCategory','sizes','overwrite',true);
end


%% AllExpAnalysis
%[49:54 57:81] MBR all experiments 'NV','NI'
%[44:56,64:88] All experiments
%[28:32,44,45,47,48,56,98] All SA experiments
%Check triggers 45, SA82 44,45,47:54,56,64:88 
% All stim: 'FFF','SDG','MBR','MB','RG','NI','NV'
%[49:54,64:97] %All PV good experiments [49:54,64:85 87:97]
% %%[89,90,92,93,95,96,97]  %Al NV and NI experiments 
%[49:54,84:90,92:96] %All SDG experiments
%solve MBR
%bootsrapRespBase

%% FIGURE 1 MOVING BALL VS STATIC BALL
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%

%% Compare MB vs RG, use gridmode true, selects maximum spatial category across directions
[tempTableMW] = AllExpAnalysis([49:54,64:66,68:85 87:97], overwrite=true,ComparePairs={'MB','RG'},PaperFig=true,...
    overwriteResponse=false,overwriteStats=true,useFDR=false,SpatialGridMode=true,maxCategory=true);

%% Calculate spatial tuning
results= SpatialTuningIndex([49:54,64:66, 68:85 87:97], indexType =  "L_amplitude_diff" ,overwrite=true...
    , topPercent = 30,useRF=true,onOff=1,unionResponsive = false,allResponsive=true, PaperFig=true, plotRFs=false, plotRFunion=false);
%% FIGURE 2 MOVING VS STATIC COMPARISON
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%

%% %% Compare SDGm vs SDGs, use gridmode true, selects maximum spatial category across directions
[tempTableMW] = AllExpAnalysis([49:54,64:66,68:85 87:97], overwrite=true,ComparePairs={'SDGm','SDGs'},PaperFig=true,...
    overwriteResponse=false,overwriteStats=true,useFDR=false,maxCategory=false,BaseRespWindow=1000);



%% FIGURE 3 SIZES AND LOCALITY COMPARISON
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
 
%% Compares MB across different sizes
[tempTableMW] = AllExpAnalysis([49:54,64:66,68:85 87:97], overwrite=false,ComparePairs={'MB'},CompareCategory="sizes",PaperFig=true,...
    overwriteResponse=false,overwriteStats=true, BaseRespWindow = 500);


%% Compares MB and SDGm across different categories, z-scores are computed with a moving window and responsive units are selected on the per category p value.
[tempTableMW] = AllExpAnalysis([49:54,64:66,68:85 87:97], overwrite=false,ComparePairs={'MB','SDGm'},CompareCategory={"directions","angles"},CompareLevels={[0,1.57,3.14, 4.71],[0,90,180,270]},PaperFig=true,...
    overwriteResponse=false,overwriteStats=true, BaseRespWindow = 1000,useGeneralFilter=true,SpatialGridMode = false,maxCategory = false);
%% Compares MB and MG, across all categories
[tempTableMW] = AllExpAnalysis([49:54,64:66,68:85 87:97], overwrite=false,ComparePairs={'MB','SDGm'},PaperFig=true,...
    overwriteResponse=false,overwriteStats=true, BaseRespWindow = 500, maxCategory = false, SpatialGridMode = false);

%% %% Compares MB and MG, across all categories
[tempTableMW] = AllExpAnalysis([49:54,64:66,68:85 87:97], overwrite=false,ComparePairs={'MB','MBR'},PaperFig=true,...
    overwriteResponse=false,overwriteStats=true, BaseRespWindow = 500, maxCategory = false, SpatialGridMode = false);

%% %% PSTH for all moving experiments
plotPSTH_MultiExp([49:54,64:66,68:85 87:97], overwrite=false, zScore=true,TakeTopPercentTrials=[], PaperFig=true, byDepth=false, smooth=50, stimTypes={"MB","MBR","SDGm"}); %stimTypes=["linearlyMovingBall"]


%% %% Compares SB, SG and FFF, across all categories
[tempTableMW] = AllExpAnalysis([49:54,64:66,68:85 87:97], overwrite=true,ComparePairs={'RG','SDGs','FFF'},PaperFig=true,...
    overwriteResponse=false,overwriteStats=true, BaseRespWindow = 300, maxCategory = false, SpatialGridMode = false);

%% PSTH for all static experiments
plotPSTH_MultiExp([49:54,64:66,68:85 87:97], overwrite=true, zScore=true,TakeTopPercentTrials=[], PaperFig=true, byDepth=false, smooth=50, stimTypes={"RG","SDGs","FFF"}); %stimTypes=["linearlyMovingBall"]

%% Plot changes in size for MB
plotPSTH_MultiExp([49:54,64:66,68:85 87:97], overwrite=true, zScore=true,TakeTopPercentTrials=[], PaperFig=true, byDepth=false, smooth=50, postStim= 2000, stimTypes={"MB"},splitBy="sizes",useCategoryPvals = true); %stimTypes=["linearlyMovingBall"]

%% SDGm spatial frequency
[tempTableMW] = AllExpAnalysis([49:54,64:66,68:85 87:97], overwrite=false,ComparePairs={'SDGs'},CompareCategory="spatFrequency",PaperFig=true,...
    overwriteResponse=false,overwriteStats=true, BaseRespWindow = 1000);

%% Plot different spatial freuqncies
plotPSTH_MultiExp([49:54,64:66,68:85 87:97], overwrite=true, zScore=true,TakeTopPercentTrials=[], PaperFig=true, byDepth=false, smooth=50, postStim= 2000, stimTypes={"SDGm"},splitBy="spatFrequency",useCategoryPvals = true); %stimTypes=["linearlyMovingBall"]


%% Plot MB raster sorted per size
plotRaster_MultiExp([49:54,64:66,68:85 87:97], overwrite=false, sortBy="preferredCategory", ...
    splitCategory="Sizes", stimTypes="MB",zScore=true,PaperFig=true)

%% Plot SDGm raster sorted per size
plotRaster_MultiExp([49:54,64:66,68:85 87:97], overwrite=true, sortBy="preferredCategory", ...
    splitCategory="spatFrequency", stimTypes="SDGm",zScore=true,PaperFig=true)

%% %% FIGURE 4 DIRECTION TUNING COMPARISON
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%

%%  Compares MB across differen directions
[tempTableMW] = AllExpAnalysis([49:54,64:66,68:85 87:97], overwrite=false,ComparePairs={'MB'},CompareCategory="sizes",PaperFig=true,...
    overwriteResponse=false,overwriteStats=true, BaseRespWindow = 500);

%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%

%% Raster for all experiment
plotRaster_MultiExp([49:54,64:66,68:85 87:97], sortBy = "spatialTuning",overwrite=true,TakeTopPercentTrials=[],PaperFig=true)

%% Calculate spatial tuning
results= SpatialTuningIndex([49:54,64:66, 68:85 87:97], indexType =  "L_amplitude_diff" ,overwrite=true...
    , topPercent = 30,useRF=true,onOff=1,unionResponsive = false,allResponsive=true, PaperFig=true, plotRFs=false, plotRFunion=false);

%% Get neuron depths
getNeuronDepths([49:54,64:97]) %[49:54,64:72,84:97] %% PV140 missing depth coordinates
%% Gratings

for ex = [97]
    NP = loadNPclassFromTable(ex); %73 81
    vs = StaticDriftingGratingAnalysis(NP);
    % vs.getSessionTime("overwrite",true);
    % vs.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    % dT = vs.getDiodeTriggers;
    % % vs.plotDiodeTriggers
    % vs.getSyncedDiodeTriggers("overwrite",true);
    %r = vs.ResponseWindow('overwrite',true);
    % results = vs.ShufflingAnalysis('overwrite',true);
    % result = vs.BootstrapPerNeuron('overwrite',true);
    % vs.StatisticsPerNeuron(overwrite=true)
    vs.plotRaster(MaxVal_1=true,OneAngle=270,exNeurons=28,AllResponsiveNeurons=false,PaperFig=false) %0.5208 %2.0833
    vs.plotRaster(MaxVal_1=false)
    close all
end

%% movie

for ex =  [92:97]
    NP = loadNPclassFromTable(ex); %73 81
    vs = movieAnalysis(NP);
    % vs.getSessionTime("overwrite",true);
    % vs.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    % dT = vs.getDiodeTriggers;
    % vs.plotDiodeTriggers
    %vs.getSyncedDiodeTriggers("overwrite",true);
    r = vs.ResponseWindow('overwrite',true);
    %results = vs.ShufflingAnalysis('overwrite',true);
    result = vs.StatisticsPerNeuron('overwrite',true);
    vs.plotRaster(AllResponsiveNeurons=true)
    close all
end


%% image

for ex = [97]
    NP = loadNPclassFromTable(ex); %73 81
    vs = imageAnalysis(NP);
    %vs.getSessionTime("overwrite",true);
    %vs.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    %dT = vs.getDiodeTriggers;
    % vs.plotDiodeTriggers
    %vs.getSyncedDiodeTriggers("overwrite",true);
    %r = vs.ResponseWindow('overwrite',true);
    %results = vs.ShufflingAnalysis('overwrite',true);
    vs.plotRaster('exNeurons',13,MergeNtrials=1,overwrite=true, selectCats =[], PaperFig=true)
    close all
    %result = vs.StatisticsPerNeuron('overwrite',true);

end


%% Moving bar
for ex = 81
    NP = loadNPclassFromTable(ex); %73 81
    vs = linearlyMovingBarAnalysis(NP);
    vs.getSessionTime("overwrite",true);
    vs.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    %vs.plotDiodeTriggers
    vs.getSyncedDiodeTriggers("overwrite",true);
    r = vs.ResponseWindow('overwrite',true);
    results = vs.ShufflingAnalysis('overwrite',true);
    if ~any(results.Speed1.pvalsResponse<0.05)
         fprintf('%d-No responsive neurons.\n',ex)
         continue
    end
    vs.CalculateReceptiveFields('overwrite',true,'nShuffle',20);
    vs.PlotReceptiveFields('overwrite',true,'RFsDivision',{'Directions','','Luminosities'},meanAllNeurons=true)
end

%% FFF
for ex = 56
    NP = loadNPclassFromTable(ex); %73 81
    vs = fullFieldFlashAnalysis(NP);
    vs.getSessionTime("overwrite",true);
    vs.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    %vs.plotDiodeTriggers
    vs.getSyncedDiodeTriggers("overwrite",true);
    r = vs.ResponseWindow('overwrite',true);
    results = vs.ShufflingAnalysis('overwrite',true);
    if ~any(results.Speed1.pvalsResponse<0.05)
         fprintf('%d-No responsive neurons.\n',ex)
         continue
    end
    vs.CalculateReceptiveFields('overwrite',true,'nShuffle',20);
    vs.PlotReceptiveFields('overwrite',true,'RFsDivision',{'Directions','','Luminosities'},meanAllNeurons=true)
end


%% Run for all
for ex = 85:88
    NP = loadNPclassFromTable(ex); %73 81
    vs = linearlyMovingBallAnalysis(NP);
    vs.getSessionTime("overwrite",true);
    vs.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true);
    %vs.plotDiodeTriggers
    vs.getSyncedDiodeTriggers("overwrite",true);
    r = vs.ResponseWindow('overwrite',true);
    results = vs.ShufflingAnalysis('overwrite',true);
    if ~any(results.Speed1.pvalsResponse<0.05)
         fprintf('%d-No responsive neurons.\n',ex)
         continue
    end
    vs.CalculateReceptiveFields('overwrite',true,'nShuffle',20);
    vs.PlotReceptiveFields('overwrite',true,'RFsDivision',{'Directions','','Luminosities'},meanAllNeurons=true)
end

%% Check experiments in timseseries viewer
timeSeriesViewer(NP)
t=NP.getTrigger;
data.VS_ordered(ex)

stimOn = t{3};
stimOff = t{4};

MBRtOn = stimOn(stimOn > t{1}(1) & stimOn < t{2}(1));
MBRtOff = stimOff(stimOff > t{1}(1) & stimOff < t{2}(1));

MBtOn = stimOn(stimOn > t{1}(2) & stimOn < t{2}(2));
MBtOff = stimOff(stimOff > t{1}(2) & stimOff < t{2}(2));

RGtOn = stimOn(stimOn > t{1}(3) & stimOn < t{2}(3));
RGtOff = stimOff(stimOff > t{1}(3) & stimOff < t{2}(3));

NGtOn = stimOn(stimOn > t{1}(4) & stimOn < t{2}(4));
NGtOff = stimOff(stimOff > t{1}(4) & stimOff < t{2}(4));

DtOn = stimOn(stimOn > t{1}(5) & stimOn < t{2}(5));
DtOff = stimOff(stimOff > t{1}(5) & stimOff < t{2}(5));

MovingBallTriggersDiode = d3.stimOnFlipTimes;



%% %% check neural data sync and analog data sync 

allTimes = [stimOn(:); stimOff(:); onSync(:); offSync(:)];  % concatenate as column

% Sort from earliest to latest
sortedTimesDiodeOldMethod = sort(allTimes);
