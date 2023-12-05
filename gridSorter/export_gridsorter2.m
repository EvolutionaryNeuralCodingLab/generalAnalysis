function export_gridsorter2(varargin)
% Function which takes JRC results and outputs t,ic,annotations

load(varargin{1});

ic = tabulate(clusterSites);
ic(ic(:,2)==0,:)=[];
ic(:,3)=[];
ic1=[];
ic2=[];
for i=1:size(ic,1)
    ic1 = [ic1 repelem(ic(i,1),ic(i,2))];
    ic2 = [ic2 1:ic(i,2)];
end
ic = [ic1;ic2];

ind = cellfun(@(x) numel(x), spikesByCluster);
ic(3,:) = cumsum([0,ind(1:end-1)])+1;
ic(4,:) = cumsum(ind);
% t = cat(1,spikeTimes(spikesByCluster{:}));
t = cell2mat(spikesByCluster');
t = spikeTimes(t);
save([varargin{1}(1:end-4) '_gridSorter.mat'],'t','ic','-v7.3');

% for iClu=1:nClu
%     if ~isempty(S0.S_clu.cviSpk_clu{iClu})
%         viSpk_clu1 = S0.S_clu.cviSpk_clu{iClu};%
%         iSite_clu1 = S0.S_clu.viSite_clu(iClu);
%         if iSite_clu1~=prevSite
%             neuCounter=1;
%         end
%         idf
%         vlSpk_clu1 = iSite_clu1 == S0.viSite_spk(viSpk_clu1);
%         %viSite_clu1 = S0.S_clu.P.miSites(:,iSite_clu1);
%
%         ic=[ic [clusterSites{iClu};neuCounter;spkCounter+1;spkCounter+numel(viSpk_clu1)]];
%         t = [t S0.viTime_spk(viSpk_clu1(vlSpk_clu1));
%         spkCounter=spkCounter+numel(t{iClu});
%         neuCounter=neuCounter+1;
%         prevSite=iSite_clu1;
%     end
% end
% t=double(cell2mat(t)')./(P.sRateHz/1000);

end