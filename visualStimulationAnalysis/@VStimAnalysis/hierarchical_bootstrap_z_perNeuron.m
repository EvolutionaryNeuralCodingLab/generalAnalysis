function zScores = hierarchical_bootstrap_z_perNeuron(obj, params)

arguments
    obj
    params.nBoot = 1000;
    params.baseline = 300;
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

if isfield(obj.VST.speed) %%Resolve speed conflict and take the mean of boths speeds
    if numel(unique(obj.VST.speed))>1
        DurationOne(obj.VST.speed)

    
    end

end

Mr = BuildBurstMatrix(goodU,round(p.t/1000),round(directimesSorted/1000),round((stimDur+ params.durationWindow)/params.binRaster)); %response matrix


nBoot = params.nBoot;




N = numel(baselineRates);
C = size(spikeRates,2);
zScores = nan(N,C);

for n = 1:N
    base = baselineRates{n};
    nBase = numel(base);
    if nBase < 2 || all(isnan(base))
        continue;
    end
    
    % Preallocate bootstrap samples for baseline
    bootMuBase = nan(nBoot,1);
    bootSdBase = nan(nBoot,1);
    
    for b = 1:nBoot
        idxB = randsample(nBase, nBase, true);
        bootVals = base(idxB);
        bootMuBase(b) = mean(bootVals, 'omitnan');
        bootSdBase(b) = std(bootVals, 'omitnan');
    end
    
    % Stimulus categories
    for c = 1:C
        stim = spikeRates{n,c};
        if isempty(stim)
            continue;
        end
        nStim = numel(stim);
        bootMuStim = nan(nBoot,1);
        for b = 1:nBoot
            idxS = randsample(nStim, nStim, true);
            bootMuStim(b) = mean(stim(idxS), 'omitnan');
        end
        
        % Compute bootstrap z-scores for this neuron & category
        bootZ = (bootMuStim - bootMuBase) ./ bootSdBase;
        zScores(n,c) = mean(bootZ, 'omitnan'); % or median(bootZ)
    end
end
end
