function results = hierarchicalBootsrapping(expList, Stims2Comp,params)

arguments
    expList  (1,:) double  %%Number of experiment from excel list
    Stims2Comp cell %% Comparison order {'MB','RG','MBR'} would select neurons responsive to moving ball and
    % compare this neurons responses to other stimuli.
    params.threshold = 0.05;
    params.nBoot = 10000; 
    params.categoriesWithinStim='Combined' %Could be direction, speed, etc
end

for ex = expList

    NP = loadNPclassFromTable(ex); %73 81
    fprintf('Processing recording: %s .\n',NP.recordingName)

    % Load Kilosort and phy results
    p = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
    label = string(p.label');
    goodU = p.ic(:,label == 'good');

    TempSpikeRates{size(goodU)}

    if ismember('MB',Stims2Comp)
        vs = linearlyMovingBallAnalysis(NP);
        A = [stimOn' obj.VST.directions' obj.VST.offsets' obj.VST.ballSizes' obj.VST.speeds' obj.VST.Luminosities'];
        [C indexS] = sortrows(A,[2 3 4 5]);

        B = [stimOff' obj.VST.directions' obj.VST.offsets' obj.VST.ballSizes' obj.VST.speeds' obj.VST.Luminosities'];
        [Coff indexSo] = sortrows(B,[2 3 4 5]);

        stimInter = obj.VST.interTrialDelay*1000;
    end

    if ismember('MB',Stims2Comp)
      
        
       
    end

    if ismember('MB',Stims2Comp)
        vs = linearlyMovingBallAnalysis(NP);
    end

    if ismember('MB',Stims2Comp)
        vs = linearlyMovingBallAnalysis(NP);
    end

    if ismember('MB',Stims2Comp)
        vs = linearlyMovingBallAnalysis(NP);
    end


    try
        DiodeCrossings = obj.getSyncedDiodeTriggers;
    catch
        obj.getSessionTime("overwrite",true);
        obj.getDiodeTriggers("extractionMethod",'digitalTriggerDiode','overwrite',true);
        DiodeCrossings = obj.getSyncedDiodeTriggers;
    end

    stimOn = DiodeCrossings.stimOnFlipTimes;
    stimOff = DiodeCrossings.stimOffFlipTimes;

    if  ~isfield(obj.VST,'Luminosities')
        obj.VST.Luminosities = zeros(1,numel(obj.VST.directions))+255;
    end

    A = [stimOn' obj.VST.directions' obj.VST.offsets' obj.VST.ballSizes' obj.VST.speeds' obj.VST.Luminosities'];
    [C indexS] = sortrows(A,[2 3 4 5]);

    B = [stimOff' obj.VST.directions' obj.VST.offsets' obj.VST.ballSizes' obj.VST.speeds' obj.VST.Luminosities'];
    [Coff indexSo] = sortrows(B,[2 3 4 5]);

    stimInter = obj.VST.interTrialDelay*1000;


end


end