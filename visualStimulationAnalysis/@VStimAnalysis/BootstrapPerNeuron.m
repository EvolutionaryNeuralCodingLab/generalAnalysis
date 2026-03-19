function results = BootstrapPerNeuron(obj, params)

arguments (Input)
    obj
    params.nBoot = 10000
    params.EmptyTrialPerc = 0.7 %If empty trials per category are higher than EmptyTrialPerc then filter
    params.FilterEmptyResponses = true
    params.overwrite = false
end
% Computes per-neuron z-scores of stimulus responses vs baseline using bootstrap


if isfile(obj.getAnalysisFileName) && ~params.overwrite
    if nargout==1
        fprintf('Loading saved results from file.\n');
        results=load(obj.getAnalysisFileName);
    else
        fprintf('Analysis already exists (use overwrite option to recalculate).\n');
    end

    return
end

p = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
label = string(p.label');
goodU = p.ic(:,label == 'good'); %somatic neurons
responseParams = obj.ResponseWindow;



if isempty(goodU)
    warning('%s has No somatic Neurons, skipping experiment/n',obj.dataObj.recordingName)
    results = [];
    fprintf('Saving results to file.\n');
     if isequal(obj.stimName, 'linearlyMovingBall')
        % S.(fieldName).BootResponse = respBoot;
        % S.(fieldName).BootBaseline = baseBoot;
        S.Speed1.BootDiff = [];
        S.Speed1.pvalsResponse = [];
        S.Speed1.ZScoreU = []; 
        S.Speed1.ObsDiff = [];
        S.Speed1.ObsReponse = [];
        S.Speed1.ObsBaseline = [];

        if isfield(responseParams, "Speed2")
            S.Speed2.BootDiff = [];
            S.Speed2.pvalsResponse = [];
            S.Speed2.ZScoreU = [];
            S.Speed2.ObsDiff = [];
            S.Speed2.ObsReponse = [];
            S.Speed2.ObsBaseline = [];
        end
    elseif isequal(obj.stimName,'StaticDriftingGrating')
        % S.(fieldName).BootResponse = respBoot;
        % S.(fieldName).BootBaseline = baseBoot;
        S.Moving.BootDiff = [];
        S.Moving.pvalsResponse = [];
        S.Moving.ZScoreU = [];
        S.Moving.ObsDiff = [];
        S.Moving.ObsReponse = [];
        S.Moving.ObsBaseline = [];

        S.Static.BootDiff = [];
        S.Static.pvalsResponse = [];
        S.Static.ZScoreU = [];
        S.Static.ObsDiff = [];
        S.Static.ObsReponse = [];
        S.Static.ObsBaseline = [];
    else
        % S.BootResponse = respBoot;
        % S.BootBaseline = baseBoot;
        S.BootDiff = [];
        S.pvalsResponse = [];
        S.ZScoreU = []; 
        S.ObsDiff = [];
        S.ObsReponse = [];
        S.ObsBaseline = [];
    end

    S.params = params;
    save(obj.getAnalysisFileName,'-struct', 'S');
    return
end

try
    DiodeCrossings = obj.getSyncedDiodeTriggers;
catch
    obj.getSessionTime("overwrite",true);
    obj.getDiodeTriggers("extractionMethod",'digitalTriggerDiode','overwrite',true);
    DiodeCrossings = obj.getSyncedDiodeTriggers;
end



%%If it is a moving stimulus with speed cathegories
if isfield(responseParams, "Speed1")

    x = length(obj.VST.speed);

    Times.Speed1 = responseParams.Speed1.C(:,1)';
    Durations.Speed1 = responseParams.Speed1.stimDur;
    trialsCats.Speed1 = numel(Times.Speed1) / size(responseParams.Speed1.NeuronVals, 2);
    
    if numel(unique(obj.VST.speed))>1
        Times.Speed2 =  responseParams.Speed2.C(:,1)';  
        Durations.Speed2 = responseParams.Speed2.stimDur;
        trialsCats.Speed2 = numel(Times.Speed2)/size(responseParams.Speed2.NeuronVals,2);
    end 

elseif isequal(obj.stimName,'StaticDriftingGrating')
    Times.Moving = responseParams.C(:,1)'+obj.VST.static_time*1000;
    Durations.Moving = responseParams.Moving.stimDur;
    trialsCat = numel(Times.Moving)/size(responseParams.Moving.NeuronVals,2);

    Times.Static = responseParams.C(:,1)';
    Durations.Static = responseParams.Static.stimDur;
    FieldNames ={'Static','Moving'};
    x=2;
else
    directimesSorted = responseParams.C(:,1)';
    stimDur = responseParams.stimDur;
    trialsCat = numel(directimesSorted)/size(responseParams.NeuronVals,2);
    x=1;
end

for s=1:x


    if isfield(responseParams, "Speed1")

        fieldName = sprintf('Speed%d', s);
        directimesSorted = Times.(fieldName);
        stimDur = Durations.(fieldName);
        trialsCat = trialsCats.(fieldName);

    end
    
    if isequal(obj.stimName,'StaticDriftingGrating')
        fieldName = FieldNames{s};
        directimesSorted = Times.(fieldName);
        stimDur = Durations.(fieldName);
 
    end

    %Mr = BuildBurstMatrix(goodU,round(p.t),round(directimesSorted),round(stimDur+ responseParams.params.durationWindow)); %response matrix
    Mr = BuildBurstMatrix(goodU,round(p.t),round(directimesSorted),round(stimDur)); %response matrix
    Mb = BuildBurstMatrix(goodU,round(p.t),round(directimesSorted-0.75*obj.VST.interTrialDelay*1000),round(0.75*obj.VST.interTrialDelay*1000)); %baseline matrix

    %Take categories where at least> 50% of trials have spikes.


    %Select top 3 trial categories.

    responses = mean(Mr,3);
    baselines = mean(Mb,3);

    if params.FilterEmptyResponses
        for i=1:trialsCat:size(Mr,1)

            for u = 1:size(goodU,2)
                tempM = responses(i:i+trialsCat-1,u);
                emptyRows = all(tempM == 0, 2);
                perc = sum(emptyRows) / size(tempM,1);

                if perc >= params.EmptyTrialPerc
                    responses(i:i+trialsCat-1, u) = zeros(1,trialsCat);
                    baselines(i:i+trialsCat-1, u) = zeros(1,trialsCat);% Store z-scores for neurons with sufficient trials
                end
            end
        end
    end


    Diff = responses - baselines;

    bootDiff = bootstrp(params.nBoot,@mean,Diff);

    pVal = mean(bootDiff <= 0); 
                    %Test the proportion of times the difference is greater or equal than 0
    bootBase = bootstrp(params.nBoot,@mean,baselines);
    stdDiff = std(bootDiff);
    
    stdBase = std(bootBase);

    z = mean(bootDiff,1) ./ stdDiff;
   

    if isfield(responseParams, "Speed1")
        % S.(fieldName).BootResponse = respBoot;
        % S.(fieldName).BootBaseline = baseBoot;
        S.(fieldName).BootDiff = bootDiff;
        S.(fieldName).pvalsResponse = pVal;
        S.(fieldName).ZScoreU = z; 
        S.(fieldName).ObsDiff = Diff;
        S.(fieldName).ObsReponse = responses;
        S.(fieldName).ObsBaseline = baselines;
    elseif isequal(obj.stimName,'StaticDriftingGrating')
        % S.(fieldName).BootResponse = respBoot;
        % S.(fieldName).BootBaseline = baseBoot;
        S.(fieldName).BootDiff = bootDiff;
        S.(fieldName).pvalsResponse = pVal;
        S.(fieldName).ZScoreU = z;
        S.(fieldName).ObsDiff = Diff;
        S.(fieldName).ObsReponse = responses;
        S.(fieldName).ObsBaseline = baselines;
    else
        % S.BootResponse = respBoot;
        % S.BootBaseline = baseBoot;
        S.BootDiff = bootDiff;
        S.pvalsResponse = pVal;
        S.ZScoreU = z; 
        S.ObsDiff = Diff;
        S.ObsReponse = responses;
        S.ObsBaseline = baselines;
    end

    S.params = params;

end
%save results in the right file
fprintf('Saving results to file.\n');
save(obj.getAnalysisFileName,'-struct', 'S');
results = S;

% respBoot_shuffled = respBoot(randperm(size(respBoot,1)),:);
% baseBoot_shuffled = baseBoot(randperm(size(baseBoot,1)),:);
% D_sh = respBoot_shuffled - baseBoot_shuffled;
% p_sh = mean(D_sh <= 0);

end
