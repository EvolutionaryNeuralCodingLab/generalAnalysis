function [R,lags]=xcorrmat(X,Y,lag)
% R=xcorrmat(X,Y,lag)
% Function purpose : Calculated the peason normalized cross correlation for 2 sets of vectors
%
% Function recives :    X [N x M] - matrix 1 containning M vectors where N is the length of the vectors
%                       Y [N x M] - matrix 2 containning M vectors where N is the length of the vectors
%                       lag [1 x 1] - the lag for cross corr calculation
%
% Function give back :  R [1 x M]- cross correlation array
%
% Last updated : 23/07/14
[N,M]=size(X);
if nargin==2
    lag=N-1;
end
if nargout==2
    lags=-lag:lag;
end

R = real(ifft(fft(X,2^nextpow2(2*N-1)) .* conj(fft(Y,2^nextpow2(2*N-1))))); %convolution in fft space
R = [R(end-lag+1:end,:);R(1:lag+1,:)]; %fft shift
R = bsxfun(@rdivide,R,sqrt(sum(abs(X.^2)).*sum(abs(Y.^2)))); %according to Parseval's theorem

%R=xcorr(X(:,1),Y(:,1),lag,'coeff');
%plot(R);hold on;plot(R(:,1),'or');
