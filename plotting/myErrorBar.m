function [h,hE,hB]=myErrorBar(x,y,dy,width,varargin);
%[hE,hB]=myErrorBar(x,y,dy,width,varargin)
%width [0-1] - the width in normalize units of the relative fraction of the bar (e.g. width=1 -> error bar in size of bar), leave empty [] if needed.
%[h,hE,hB]=myErrorBar(1:8,Int,Int+IntStd,0.4,'Bfacecolor',[0.3 0.3 0.3],'BLinewidth',5,'Ecolor','k');
if isempty(width)
    width=0.25;
end
if isempty(x)
    x=1:length(y);
end
%distribute varargin between those for the bar plot and those for the errorbar plot
barVarArgIn={};
errorVarArgIn={};
axisVarArgIn={};
for i=1:2:length(varargin)
    if varargin{i}(1)=='B'
        barVarArgIn={barVarArgIn{:},varargin{i}(2:end),varargin{i+1}};
    elseif varargin{i}(1)=='E'
        errorVarArgIn={errorVarArgIn{:},varargin{i}(2:end),varargin{i+1}};
    else
        axisVarArgIn={axisVarArgIn{:},varargin{i}(1:end),varargin{i+1}};
    end
end

%checks if handle to axis exists and if not creates axis 
if exist('axisHandle','var')
    h=axisHandle;
else
    h=gca;
end
hold on;
nSamples=length(x);

%generate bar plot
hB=bar(h,x,y,barVarArgIn{:});

upperBars=(y>=0);
lowerBars=(y<0);
U=zeros(1,nSamples);L=zeros(1,nSamples);
U(upperBars)=dy(upperBars);L(lowerBars)=dy(lowerBars);

%generate error bars
hE=errorbar(h,x,y,L,U,'Marker','none','linestyle','none','MarkerSize',eps,errorVarArgIn{:});
ylim(ylim*1.02); %extend ylimits

%{
%detect the bar width from bar plot
xD = hB.XData;
w = (xD(3,1)-xD(2,1))/2;

% Retrieve info from errorbar plot	
xD = hE.XData;		% Get xdata from errorbar plot

% Change xdata in the error bar plot
N=length(xD);
upperBarIdx{1}=4:9:N;upperBarIdx{1}=upperBarIdx{1}(upperBars);
upperBarIdx{2}=5:9:N;upperBarIdx{2}=upperBarIdx{2}(upperBars);
upperBarIdx{3}=7:9:N;upperBarIdx{3}=upperBarIdx{3}(upperBars);
upperBarIdx{4}=8:9:N;upperBarIdx{4}=upperBarIdx{4}(upperBars);
lowerBarIdx{1}=4:9:N;lowerBarIdx{1}=lowerBarIdx{1}(lowerBars);
lowerBarIdx{2}=5:9:N;lowerBarIdx{2}=lowerBarIdx{2}(lowerBars);
lowerBarIdx{3}=7:9:N;lowerBarIdx{3}=lowerBarIdx{3}(lowerBars);
lowerBarIdx{4}=8:9:N;lowerBarIdx{4}=lowerBarIdx{4}(lowerBars);

xD(upperBarIdx{1}) = xD(upperBarIdx{1}-3)-w*width;
xD(upperBarIdx{2}) = xD(upperBarIdx{2}-4)+w*width;
xD(upperBarIdx{3}) = xD(upperBarIdx{3}-6);
xD(upperBarIdx{4}) = xD(upperBarIdx{4}-7);
xD(lowerBarIdx{1}) = xD(lowerBarIdx{1}-3);
xD(lowerBarIdx{2}) = xD(lowerBarIdx{2}-4);
xD(lowerBarIdx{3}) = xD(lowerBarIdx{3}-6)-w*width;
xD(lowerBarIdx{4}) = xD(lowerBarIdx{4}-7)+w*width;

set(hhE(2),'xdata',xD(:));	% update error bars on the figure
%}

%update axis arguments
set(h,'xtick',x);
set(h,axisVarArgIn{:});





        
