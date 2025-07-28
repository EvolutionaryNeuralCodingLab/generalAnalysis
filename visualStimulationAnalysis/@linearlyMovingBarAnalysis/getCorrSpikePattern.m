function result = getCorrSpikePattern(obj,varargin)
%Plots the correlation matrix between the responses of all pairs of stimlui
T = obj.getSyncedDiodeTriggers;

%order trials accoridng ot direction of movement and offset for each direction.
[~,pOrdered2]=sort(obj.VST.offsets);
[originalOrder,pOrdered1]=sort(obj.VST.directions(pOrdered2));
pOrdered=pOrdered2(pOrdered1);
trialCat="Dir="+ num2str(obj.VST.directions(pOrdered)',2)+",offset=" + num2str(obj.VST.offsets(pOrdered)',4);

result = getCorrSpikePattern@VStimAnalysis(obj,T.stimOnFlipTimes(pOrdered),trialCat,'win',obj.VST.stimDuration*1000,varargin{:});

end