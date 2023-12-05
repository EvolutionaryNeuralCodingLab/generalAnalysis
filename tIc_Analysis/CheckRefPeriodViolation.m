function [t_new,ic_new,V_new]=CheckRefPeriodViolation(t,ic,ref,SeeRefHist,V)
% [t_new,ic_new,V_new]=CheckRefPeriodViolation(t,indexchannel,ref,SeeRefHist,V);
%
% Function purpose : 1)Calculates the number of times in which differences between spikes
%                      are more than ref and plot a histogram of refractory period violations.
%                    2)Eliminates violations from data.
%                     3)Elliminates NaNs from data
% Function recives :    t - firing timings ([ms]
%                       ic - indexchannel
%                       ref [ms] - the refracory period for spike elimination from t,ic
%                       in  (defult valure ref=0).
%                       SeeRefHist='y' to see histograms and 'n' to not see
%                       V - the voltages of another vector corresponding to t
% Function give back :  t_new [ms] - sorted t with no violation of refractory period.
%                       ic_new - non violating indexchannel
%                       V_new - non violating V vector
% *indexchannel must be in order in times.
%
% Last updated : 06/01/10
% Remark: The use of ifexist('V','var') should be optimized

if nargin==2
    ref=0;
end

fprintf('Calculating Refractory Peroid violations - set to %f ms \n ---- \n',ref);

for i=1:size(ic,2)
    fails=length(find(diff(sort(t(ic(3,i):ic(4,i))))<=ref));
    tot=length(t(ic(3,i):ic(4,i)));
    fprintf('Channel:%d neuron: %d - %d out of %d  (%f %%) \n',ic(1,i),ic(2,i),fails,tot,fails/tot*100);
end


if (ref~=0 && ~exist('SeeRefHist','var'))
    SeeRefHist=input('Do you want to see the violations histograms (y/n) ? ','s');
else
    if ~exist('SeeRefHist','var')
        SeeRefHist='n';
    end
end
if SeeRefHist=='y'
    for i=1:size(ic,2)
        diffs=diff(sort(t(ic(3,i):ic(4,i))));
        small_diffs=diffs(find(diffs<=6));
        if ~isempty(small_diffs)
            IntervalHist(i,:)=histc(small_diffs,0:0.25:6);
            bar(0:0.25:6,IntervalHist(i,:),'hist');
            title(['channel ', int2str(ic(1,i)),' neuron ', int2str(ic(2,i))]);
            pause;
        end
    end
    close;
end

count=1;
ic_new=ic(1:2,:);
t_new=[];
V_new=[];
count=1;
for i=1:size(ic,2)
    ic_new(3,i)=count;
    if exist('V','var')
        channel_V=V(ic(3,i):ic(4,i));
    end
    channel_t=t(ic(3,i):ic(4,i));
    if ~issorted(t(ic(3,i):ic(4,i)))
        [channel_t P]=sort(channel_t);
        disp(['Warning - the input t vector was not sorted in channel: ' num2str(ic_new(1,i))]);
        if exist('V','var')
            channel_V=channel_V(P);
        end
    end
    if isempty(find(channel_t>=0))
        channel_t=abs(channel_t);
        disp(['Warning - the input t vector had non positive values in channel: ' num2str(ic_new(1,i))]);
    end
    ToDelete=find(diff(channel_t)<=ref)+1;
    channel_t(ToDelete)=[];
    if exist('V','var')
        channel_V(ToDelete)=[];
    end
    NanPlaces=find(isnan(channel_t));
    if ~isempty(NanPlaces)
        channel_t(NanPlaces)=[];  %Elliminates NaN from Data
        disp(['Warning - the input t vector had NaN values in channel: ' num2str(ic_new(1,i))]);
        if exist('V','var')
            channel_V(NanPlaces)=[];  %Elliminates NaN from Data
        end
    end
    if exist('V','var')
        V_new=[V_new channel_V];
    end
    t_new=[t_new channel_t];
    ic_new(4,i)=length(t_new);
    count=ic_new(4,i)+1;
end