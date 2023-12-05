s=randn(1,256); %256 elements, Fs = 1Hz.
s=sin(2*pi*10*(0.01:0.01:2.56));
numlevels = wmaxlev(length(s), 'haar');
numlevels=4;
[C,L]= wavedec(s, numlevels, 'haar');
D = detcoef(C,L,'cells');
[Ea,Ed] = wenergy(C,L);

T = wptree(2,10,s,'haar');
plot(T)

[phi,psi,xval] = wavefun('coif3',8); 
plot(xval,psi);

[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters('coif1')

subplot(1.3,1,1);
PlotPacketBasis(ijk,abs(wpd).^0.5);
if nargin==2
    fsig = Fs*(0:0.5:(Lt/2))/Lt;%this freq. may not be accurate, should be checked
    jumpFactor=round(length(fsig)/10);
    set(gca,'YTick',0:jumpFactor:Lt,'YTickLabel',fsig(1:jumpFactor:end));
    ylabel('Freq [Hz]');
end
