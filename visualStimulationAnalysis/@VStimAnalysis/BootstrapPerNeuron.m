function results = BootstrapPerNeuron(obj, params)

arguments (Input)
    obj
    params.nBoot = 10000
    params.EmptyTrialPerc = 0.6
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

try
    DiodeCrossings = obj.getSyncedDiodeTriggers;
catch
    obj.getSessionTime("overwrite",true);
    obj.getDiodeTriggers("extractionMethod",'digitalTriggerDiode','overwrite',true);
    DiodeCrossings = obj.getSyncedDiodeTriggers;
end


responseParams = obj.ResponseWindow;

%%If it is a moving stimulus with speed cathegories
if isfield(responseParams, "Speed1")

    x = length(obj.VST.speed);

    % % Generate field names for each unique speed
    % topFields = arrayfun(@(n) sprintf('Speed%d', n), 1:x, 'UniformOutput', false);
    % 
    % % Initialize each as an empty struct with subfields
    % subFields = {'pvalsResponse','ZScoreU','boot_means'};  % subfield names
    % emptySubStruct = cell2struct(cell(1, numel(subFields)), subFields, 2);
    % S = cell2struct(repmat({emptySubStruct}, 1, x), topFields, 2);

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

    Mr = BuildBurstMatrix(goodU,round(p.t),round(directimesSorted),round(stimDur+ responseParams.params.durationWindow)); %response matrix
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

    % figure;imagesc(responses(:,60))
    %
    % figure;imagesc(baselines(:,60))
    %
    % figure;imagesc(squeeze(Mr(:,60,:)))
    %
    % figure;imagesc(baselines(:,60))


    % respBoot = bootstrp(params.nBoot,@mean,responses);
    % baseBoot =  bootstrp(params.nBoot,@mean,baselines);

    Diff = responses - baselines;

    bootDiff = bootstrp(params.nBoot,@mean,Diff);

    pVal = mean(bootDiff <= 0); 
                    %Test the proportion of times the difference is greater or equal than 0
    bootBase = bootstrp(params.nBoot,@mean,baselines);
    stdDiff = std(bootDiff);
    
    stdBase = std(bootBase);



    % nNeurons = size(responses,2);
    % 
    % p_hat = zeros(nNeurons,1);
    % p_dist_all = zeros(params.nBoot, nNeurons);   % optional
    % 
    % statFun2D = @(x) mean(x(:,2) > x(:,1));  % dominance statistic 
    % 
    % for nn = 1:nNeurons
    % 
    %     data = [baselines(:,nn), responses(:,nn)];
    % 
    %     p_dist = bootstrp(params.nBoot, statFun2D, data);
    %     p_dist_all(:,nn) = p_dist;  % store if needed
    % 
    %     p_hat(nn) = mean(p_dist);
    % end
    % params.alpha= 0.05;
    % % Decision rule like in the paragraph
    % alpha = params.alpha;
    % decision = strings(nNeurons,1);
    % 
    % for nn = 1:nNeurons
    %     if p_hat(nn) > 1 - alpha/2
    %         decision(nn) = "response > baseline";
    %     elseif p_hat(nn) < alpha/2
    %         decision(nn) = "baseline > response";
    %     else
    %         decision(nn) = "no difference";
    %     end
    % end
    % figure;
    % h = histogram(respBoot(:,1));
    % h.BinEdges
    % hold on;
    % histogram(baseBoot(:,1),'BinEdges',h.BinEdges)

    % D = respBoot - baseBoot; %hypothesis testing
    % 
    % pVal = mean(D <= 0); %by chance how often the baseline mean is bigger than the reponse.
    % 
    z = mean(bootDiff,1) ./ stdDiff;
    %z(isnan(z)) = 0;
 
    % 
    % compare2 = [p_hat,pHat', pVal'];



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
