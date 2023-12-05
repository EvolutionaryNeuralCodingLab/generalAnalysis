function filtsig = butfilt(sig, Fcp, Fsp, order)
% Fcp=200;  %cutoff frequency
% order = 2;
% Fsp  = 20000;
[z,p,k]=butter(order,Fcp/(Fsp/2),'high');
[sos,g] = zp2sos(z,p,k);
filtsig=filtfilt(sos,g,sig);
end