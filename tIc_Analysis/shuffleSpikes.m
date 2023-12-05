function [tOut,icOut]=shuffleSpikes(t,ic)
networkLevelShuffle=0;
nNeu=size(ic,2);
tOut=t;
icOut=ic;
if networkLevelShuffle %shuffles the spikes between neurons
    tOut=t(randperm(numel(t)));
else %shuffle times within neurons
    for i=1:nNeu
        tmpT=diff(t(ic(3,i):ic(4,i)));
        tOut(ic(3,i):ic(4,i))=[t(ic(3,1)) cumsum(tmpT(randperm(numel(tmpT))))];
    end
end