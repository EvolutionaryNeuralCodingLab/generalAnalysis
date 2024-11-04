function [Psig,fsig,m0,m1,m2]=FourierAnalysisPlot(signal,Fs,varargin)
%[Psig,fsig,m0,m1,m2]=FourierAnalysisPlot(signal,Fs,varargin)
% Function purpose : Calculates the fourier analysis
%
% Recives :             signal - the input signal
%                                   Fs - sampling frequency [Hz]
%                       Options/varargin (format - 'option',value)
%                                   'method' - 'fft',periodogram','pwelch','pyulear'
%                                   'plotResults' - [0/1] (default=1) - if==1 generates a plot of the results
%                                   'avgSingals'- [0/1] (default=1) - if==1 and if input contained several channels (matrix) - the results is averaged over all channels.
%                                   'axisHandles' - a vector with two axis handles, one for the signal and one for the spectrum
% Gives back:       Psig - the power spectrum
%                                   fsig - the frequencies corresponding to the power spectrum
%                                   m0-the zero moment of the spectrum
%                                   m1-the first moment of the spectrum
%                                   m2-the second moment of the spectrum
%                                                                                                                                                        
% Recomended usage  : [Psig,fsig,m0,m1,m2]=FourierAnalysisPlot(X,1e-3,'plotResults',0)
% Last updated : 30/05/2011

plotResults=1;
avgSingals=1;
method='fft';
%collect option variables
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};']);
end

%invert one dimensional signals and extract signal parameters
if size(signal,2)==1
    signal=signal';
end
[nsignals nsamples]=size(signal);

%calculate the fourier spectrum
for i=1:nsignals
    switch method
        case 'fft'
            Ysig = fft(signal(i,:));
            Psig_tmp = Ysig.*conj(Ysig);
            Psig(i,:) = Psig_tmp(1:(nsamples/2+1));
            fsig = Fs*(0:(nsamples/2))/nsamples;
        case 'periodogram'
            [Ysig,fsig] = periodogram(signal(i,:),[],[],Fs);
        case 'pwelch'
            [Ysig,fsig] = pwelch(signal(i,:),[],[],[],Fs);
        case 'pyulear'
            p=100;
            [Ysig,fsig] = pyulear(signal(i,:),p,[],Fs);
    end
end

%calculate moments
if avgSingals
    Psig=mean(Psig,1);
    signal=mean(signal,1);
    m0=sum(Psig);
    m1=sum(Psig.*(fsig));
    m2=sum(Psig.*(fsig.^2));
else
    m0=sum(Psig,2);
    m1=sum(Psig.*(ones(nsignals,1)*fsig),2);
    m2=sum(Psig.*(ones(nsignals,1)*(fsig.^2)),2);
end

%plot results
if plotResults==1
    if exist('axisHandles','var')
        axes(axisHandles(1));
    else
        subplot(10,1,1:6);
    end
    plot(fsig,Psig(:,1:(nsamples/2+1))./((m0/2)*ones(1,nsamples/2+1)),'.-');
    xlabel('frequency [Hz]','FontSize',14);
    ylabel('PSD','FontSize',14);
    legend('show');
    
    if exist('axisHandles','var')
        axes(axisHandles(2));
    else
        subplot(10,1,8:10);
    end
    Times=(1:size(signal,2))/Fs;
    plot(Times,signal);
    xlim([Times(1) Times(end)]);
    xlabel('Time [s]','FontSize',14);
    ylabel('Signal','FontSize',14);    
end
