function results = BootstrapPerNeuron(obj, params)

arguments
    obj
    params.nBoot = 10000;
end
% Computes per-neuron z-scores of stimulus responses vs baseline using bootstrap
%
% INPUT:
%   spikeRates{n,c}  = trial vectors (firing rates during stimulus)
%   baselineRates{n} = trial vector (baseline rates)
%   nBoot            = number of bootstrap iterations (e.g. 1000)
%
% OUTPUT:
%   zScores(n,c)     = bootstrapped z-score per neuron & condition

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

stimOn = DiodeCrossings.stimOnFlipTimes;
stimOff = DiodeCrossings.stimOffFlipTimes;

responseParams = obj.ResponseWindow;

%%Select fastest speed for moving ball to determine neuron statistic
if isfield(responseParams, "Speed1")
    directimesSorted = responseParams.Speed1.C(:,1)';
    stimDur = responseParams.Speed1.stimDur;
    trialsCat = numel(directimesSorted) / size(responseParams.Speed1.NeuronVals, 2);
    
    if numel(unique(obj.VST.speed))>1
        directimesSorted =  responseParams.Speed2.C(:,1)';  
        stimDur = responseParams.Speed2.stimDur;
        trialsCat = numel(directimesSorted)/size(responseParams.Speed2.NeuronVals,2);
    end 
else
    directimesSorted = responseParams.C(:,1)';
    stimDur = responseParams.stimDur;
    trialsCat = numel(directimesSorted)/size(responseParams.NeuronVals,2);
end

Mr = BuildBurstMatrix(goodU,round(p.t),round(directimesSorted),round(stimDur+ responseParams.params.durationWindow)); %response matrix
Mb = BuildBurstMatrix(goodU,round(p.t),round(directimesSorted-0.75*obj.VST.interTrialDelay*1000),round(0.75*obj.VST.interTrialDelay*1000)); %baseline matrix

%Take categories where at least> 50% of trials have spikes.


%Select top 3 trial categories. 

responses = mean(Mr,3);
baselines = mean(Mb,3);

for i=1:trialsCat:size(Mr,1)

    for u = 1:size(goodU,2)
        tempM = responses(i:i+trialsCat-1,u);
        emptyRows = all(tempM == 0, 2);
        perc = sum(emptyRows) / size(tempM,1);

        if perc > 0.5
            responses(i:i+trialsCat-1, u) = zeros(1,trialsCat); % Store z-scores for neurons with sufficient trials
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


respBoot = bootstrp(params.nBoot,@mean,responses);
baseBoot =  bootstrp(params.nBoot,@mean,baselines);


% figure;
% h = histogram(respBoot(:,1));
% h.BinEdges
% hold on;
% histogram(baseBoot(:,1),'BinEdges',h.BinEdges)

D = respBoot - baseBoot; %hypothesis testing 

pVal = mean(D <= 0); %by chance how often the baseline mean is bigger than the reponse.

z = mean(D,1) ./ std(D,[],1);

S.BootResponse = respBoot;
S.BootBaseline = baseBoot;
S.pValues = pVal;
S.Z.scores = z; 

%save results in the right file
fprintf('Saving results to file.\n');
save(obj.getAnalysisFileName,'-struct', 'S');
results = S;

% respBoot_shuffled = respBoot(randperm(size(respBoot,1)),:);
% baseBoot_shuffled = baseBoot(randperm(size(baseBoot,1)),:);
% D_sh = respBoot_shuffled - baseBoot_shuffled;
% p_sh = mean(D_sh <= 0);

end
