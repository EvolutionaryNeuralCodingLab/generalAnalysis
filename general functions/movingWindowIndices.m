function [startIdx,endIdx,allIdx,tWindowCenter_sec]=movingWindowIndices(Nsamples,win_sec,Fs_Hz,OL_sec)
% [startIdx,endIdx,tWindowCenter_sec]=movingWindowIndices(Nsamples,win_sec,Fs_Hz,OL_percent);
% Function purpose : calculate indices for moving window analysis
%
% Recives : Nsamples - the number of samples in the signal
%                       win_sec - the window duration [sec]
%                       Fs_Hz - the signal's sampling frequency
%                       OL_sec - the degree of overlap [sec]
%                                                                                                                                                        
% Function give back :  startIdx - the start indices of the windows
%                                                   endIdx - the end indices of the windows
%                                                   allIdx - all the indices of the moving window (each column corresponds to one moving window)
%                                                   tWindowCenter_sec - the times of window centers [sec]
%
% Recomended usage  : [startIdx,endIdx,allIdx,tWindowCenter_sec]=movingWindowIndices(1e6,2e-6,800e6,400e6)
% Last updated : 13/11/12

%translate window size and overlap to samples and check that the window size corresponds to the sampling rate 
win_samples=win_sec*Fs_Hz;
OL_samples=OL_sec*Fs_Hz;
if win_samples~=round(win_samples) 
    win_samples=round(win_samples);
    disp(['Due to mismatch between window size and sampling freq., the # of samples in window was rounded to: ' num2str(win_samples)]);
end

if Nsamples/Fs_Hz<win_sec %verify that the window is shorter than the signal
    disp(['Error: the entered window is longer than the data set']);
    startIdx=[];endIdx=[];allIdx=[];tWindowCenter_sec=[];
    return;
end
% determine the number of segments after parsing
nSegments = floor((Nsamples-OL_samples)/(win_samples-OL_samples));

% create parsing indices
startIdx = 1 + (0:(nSegments-1))*(win_samples-OL_samples);
endIdx = startIdx+win_samples-1;
tWindowCenter_sec=(startIdx+win_samples/2)/Fs_Hz; %real time [seconds] window centers

%checks that the chosen window enters an integer number of times in signal (if not the signal is trancated)
if (startIdx(end)+win_samples)<Nsamples
    disp(['The chosen win length and overlap do not correspond to signal length -> length cut to ' num2str(startIdx(end)+win_samples-1) ' samples']);
end

% allocate segmented matrices
allIdx = zeros(win_samples,nSegments);

% parse matrix
Idx = (1:win_samples)';
allIdx(:) = Idx(:,ones(1,nSegments))+startIdx(ones(win_samples,1),:)-1;
