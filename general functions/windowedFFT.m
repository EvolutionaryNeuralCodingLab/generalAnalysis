function [pSig,fSig,tSig]=windowedFFT(V,Fs,win,OL,varargin)
%[pSig,fSig,tSig]=windowedFFT(V,Fs,win,OL,varargin)
% Function purpose : Calculates the fourier transform in consecutive/overlapping time windows
%
% Recives :             V - the input signal
%                       Fs - sampling frequency [Hz]
%                       win - window in time units [s].
%                       OL - overlap [s]
%                       Options/varargin (format - 'option',value)
%                           FORCE_POWER2_SEQUENCES - optimize input windows to be an integer power of two in length
%                           NORMALIZATION - normalize spectrum in every window - 
%                                0=no normalization, 
%                                1=normalize to the fft sum, 
%                                2=normalize to the singal square amplitude
%                                3=normalize to the numner of samples N
%                                4=normalize to the square of samples N^2
%                           METHOD - method for power spectrum estimation - 'FFT','MatlabPeriodogram','MatlabWelch'
%                           
% Gives back:       pSig - the power spectrum in time windows
%                   fSig - the frequencies corresponding to the power spectrum
%                   tSig - the timing of each window [s]
%                                                                                                                                                        
% Recomended usage  : 
% Last updated : 13/11/2012

%paramters
FORCE_POWER2_SEQUENCES=0;
NORMALIZATION=0;
METHOD='FFT';
%collect input arguments
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};']);
end

nSamples=length(V);
%check that window contains an integer number of samples
if win*Fs~=round(win*Fs)
    win=round(win*Fs)/Fs;
    disp(['!!!Window size was changed to: ' num2str(win)]);
end
samplesWin=win*Fs;
%check that window overlap contains an integer number of samples
if OL*Fs~=round(OL*Fs)
    OL=round(OL*Fs)/Fs;
    disp(['!!!Window overlap was changed to: ' num2str(OL)]);
end
samplesOL=OL*Fs;
%Examine is window is an equal power of two
if FORCE_POWER2_SEQUENCES
    win=2^round(log2(win));
    OL=2^round(log2(OL));
    disp(['!!!The window size was changed to: ' num2str(win) ' and overlap to ' num2str(OL)]);
else
    if (2^floor(log2(win)))~=win && round(nSamples/OL)~=(nSamples/OL)
        disp('!!!Input data vector is not of 2^N length (N=1,2,3,...)');
    end
end

%buffering data into windows
[bufferedV,residualSamples] = buffer(V,samplesWin,samplesOL,'nodelay');
nWindows=size(bufferedV,2);

%generating hamming window
HammingWindow=hamming(samplesWin,'periodic')';
HammingWindow=HammingWindow/sum(HammingWindow);

[pSig,fSig,m0]=FourierAnalysis((HammingWindow'*ones(1,nWindows)).*bufferedV,Fs,'plotResults',0,'method',METHOD);
nFFT=length(fSig);
switch NORMALIZATION
    case 0
    case 1
        pSig=pSig./(ones(nFFT,1)*sum(pSig));
    case 2
        pSig=pSig./(ones(nFFT,1)*sum(bufferedV.^2));
    case 3
        pSig=pSig/samplesWin;
    case 4
        pSig=pSig/(samplesWin.^2);
    case 5
        pSig=pSig./(ones(nFFT,1)*((m0)./(samplesWin.^2)));
end

if nargin>=3
    tSig=(round(samplesWin/2):(samplesWin-samplesOL):(nSamples-samplesWin+samplesOL))/Fs;
end

