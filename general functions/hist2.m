function [D,h,cb]=hist2(X1,X2,varargin)
% [D,h,cb]=hist2(X1,X2,varargin)
% Function purpose : Calculates and plots a 2D histogram
%
% Function recives :    X1 - first axis coordinates
%                       X2 - second axis coordinates
%                       varargin - arguments in format: 'property',value,...
%                           logColorScale - [default=true] whether to use a logarithmic color scale
%                           h - axis handle for plotting
%                           nBins1 - [1x1] - number of bins for coordinate 1
%                           nBins2 - [1x1] - number of bins for coordinate 2
%                           dX1 - [1x1] - interval for coordinate 1
%                           dX2 - [1x1] - interval for coordinate 2
%                           edges1 - [1,M] - bin centers for coordinate 1
%                           edges2 - [1,M] - bin centers for coordinate 2
%                           plotProb [default=false] - display the count histrogram or the probability 
%
% Function give back :  D - the count histogram
%                       h - a handle to the density plot (if not entered only calculates and doesnt plot)
%
% Last updated : 27/06/12

%set default variables
nBins1=100;
nBins2=100;
dX1=[];
dX2=[];
edges1=[];
edges2=[];
h=[];
plotResults=1;
logColorScale=1;
plotColorBar=1;
plotProb=0;
limits1=[];
limits2=[];
removeZeros=false;


% Output list of default variables
%print out default arguments and values if no inputs are given
if nargin==0
    defaultArguments=who;
    for i=1:numel(defaultArguments)
        eval(['defaultArgumentValue=' defaultArguments{i} ';']);
        if numel(defaultArgumentValue)==1
            disp([defaultArguments{i} ' = ' num2str(defaultArgumentValue)]);
        else
            fprintf([defaultArguments{i} ' = ']);
            disp(defaultArgumentValue);
        end
    end
    return;
end

% Collects all input variables
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end

%check if histogram limits were provided as inputs
if isempty(limits1)
    limits1(1)=min(X1(:));
    limits1(2)=max(X1(:));
else %if limits were added as a parameter - reject points outside limits
    pReject=X1<limits1(1) | X1>limits1(2);
    X1(pReject)=[];
    X2(pReject)=[];
end
if isempty(limits2)
    limits2(1)=min(X2(:));
    limits2(2)=max(X2(:));
else %if limits were added as a parameter - reject points outside limits
    pReject=X2<limits2(1) | X2>limits2(2);
    X1(pReject)=[];
    X2(pReject)=[];
end

if isempty(edges1) %coordinate 1
    if isempty(dX1) %divide hist to n bins
        if limits1(1)~=limits1(2)
            edges1=linspace(limits1(1),limits1(2),nBins1);
        else
            edges1=linspace(limits1(1)-1,limits1(2)+1,nBins1);
            disp('There is only one value on the X dimension, cant plot a 2D hist');
        end
    else %divide hist to sections of length dX
        edges1=limits1(1):dX1:(limits1(2)+dX1);
        edges1=edges1+dX1/2;
    end
end
if isempty(edges2) %coordinate 1
    if isempty(dX2) %divide hist to n bins
        if limits2(1)~=limits2(2)
            edges2=linspace(limits2(1),limits2(2),nBins2);
        else
            edges2=linspace(limits2(1)-1,limits2(2)+1,nBins2);
            disp('There is only one value on the Y dimension, cant plot a 2D hist');
        end
    else %divide hist to sections of length dX
        edges2=limits2(1):dX2:(limits2(2)+dX2);
        edges2=edges2+dX2/2;
    end
end

nBinBins1=numel(edges1);
nBinBins2=numel(edges2);

%interpolate to grid
x1r=interp1(edges1,1:numel(edges1),X1(:),'nearest','extrap');
x2r=interp1(edges2,1:numel(edges2),X2(:),'nearest','extrap');

D=accumarray([x1r x2r],1,[nBinBins1 nBinBins2])';
if plotProb
    D=D./sum(D);
end

if removeZeros
    zeroX=find(edges1>0,1,'first');
    zeroY=find(edges2>0,1,'first');
    D(zeroY,zeroX)=0;
end

if plotResults %plot results
    if isempty(h)
        h=axes;
    else
        axes(h);
    end
    colormap(flipud(gray(256)));
    if logColorScale
        imagesc(edges1,edges2,log10(1+D));
    else
        imagesc(edges1,edges2,D);
    end
    if plotColorBar
        cb=colorbar('location','EastOutside');
        if logColorScale
            if plotProb
                ylabel(cb,'log_{10}(1+Prob.)');
            else
                ylabel(cb,'log_{10}(1+#)');
            end
        else
            if plotProb
                ylabel(cb,'Prob.');
            else
                ylabel(cb,'#');
            end
            
        end
    else
        cb=[];
    end
    set(h,'YDir','normal');
end