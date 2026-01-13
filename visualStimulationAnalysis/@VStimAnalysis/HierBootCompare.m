function results = HierBootCompare(obj, data, insertions, animals, params)

arguments (Input)
    obj
    data % Commbined 1D vector of measurements (z-score, spike rate, etc)
    insertions %Labels for insertions
    animals %Labels for animals
    params.nBoot = 10000
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

for i= 1:numel(data) %hier bootstrap stimuls
    bootstrapped = hierBootMatchFreq;
end

% S.BootResponse = respBoot;
% S.BootBaseline = baseBoot;
S.BootDiff = bootDiff;
S.pvalsResponse = pVal;
S.ZScoreU = z;


S.params = params;

%save results in the right file
fprintf('Saving results to file.\n');
save(obj.getAnalysisFileName,'-struct', 'S');
results = S;

end


