function [h]=RawDataShiftedPlot(M,Fs,h)
% [h]=RawDataShiftedPlot(M,Fs,h)
% Function purpose : Plots multiple traces on a 3d plot
%
% Recives :    M (m X n matrix) - with m data traces, n length each
%                           Fs - sampling frequency
%                           h - axes handle
%                                                                                                                                                        
% Function give back :  h - axes handle
% Last updated : 9/01/10

TRes=1000; %the units of the time axes: 1000-for [ms], 1-for [s]
if ~exist('h','var')
    figure;h=axes;hold on;
else
    axis(h);hold on;
end

cmap=colormap('Lines');
[NTraces TraceLen]=size(M);
T=(1:(TraceLen))./Fs.*TRes; % to show traces in [ms] for seconds use
O=ones(1,TraceLen);
for i=1:NTraces
    plot3(h,T,O.*i,M(i,:),'color',cmap(rem(i-1,5)+1,:));
end

set(gca,'CameraPosition',[-42  -14  205]);
xlabel('Time [ms]','FontSize',14);
ylabel('# SBE','FontSize',14);
zlabel('Voltage [\muV]','FontSize',14);

MinV=min(min(M));
MaxV=max(max(M));
zlim([MinV-0.01*MinV MaxV+0.01*MaxV]);