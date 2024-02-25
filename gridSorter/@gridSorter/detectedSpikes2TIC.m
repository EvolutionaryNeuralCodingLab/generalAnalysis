function [gridObj] = detectedSpikes2TIC(gridObj,varargin)
%DETECTEDSPIKES2TIC retrieves all spike times from the spike detection
%and converts it into a t,IC format
%   detectedSpikes2TIC loads spike times from gridObj.sortingFileNames, and 
%   saves the output file to [gridObj.sortingDir '\GridSorterDetectedSpikes.m'] unless given 
%   otherwise in varargin (use 'ticPath' key,value to change this)
%   NOTICE: changing these variables was not QA-ed, so use with caution

ticPath=[gridObj.sortingDir filesep 'GridSorterDetectedSpikes.mat'];

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

channels=gridObj.selectedChannelSubset;
nCh=length(channels);

ic(1,:)=channels;
ic(2,:)=ones(1,nCh); %with no sorting - assuming 1 neuron per channel
t=[];

for i=1:nCh
    spikeTimes=load(gridObj.sortingFileNames.spikeDetectionFile{i},'spikeTimes');
    try
        spikeTimes=spikeTimes.spikeTimes;
    catch
        spikeTimes=[];
    end
    ic(3,i)=length(t)+1;
    ic(4,i)=length(t)+length(spikeTimes);
    t=[t spikeTimes];
end

%remove channels without spikes
ic(:,ic(3,:)==(ic(4,:)+1))=[];

save(ticPath,'t','ic')
