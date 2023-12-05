function [D,h]=hist3(X1,X2,X3,varargin)
%[D,h]=hist3(X,Y,varargin)
% Function purpose : Calculates and plots a 3D histogram
%
% Function recives :    X1 - first axis coordinates
%                       X2 - second axis coordinates
%                       X3 - third axis coordinates
%                       varargin - arguments in format: 'property',value,...
%                           n1 - [1x1] - number of bins for coordinate 1
%                           n2 - [1x1] - number of bins for coordinate 2
%                           n3 - [1x1] - number of bins for coordinate 3
%                           dX1 - [1x1] - interval for coordinate 1
%                           dX2 - [1x1] - interval for coordinate 2
%                           dX3 - [1x1] - interval for coordinate 3
%                           c1 - [1,M] - bin centers for coordinate 1
%                           c2 - [1,M] - bin centers for coordinate 2
%                           c3 - [1,M] - bin centers for coordinate 3
%
% Function give back :  D - the count histogram
%                       h - a handle to the density plot (if not entered only calculates and doesnt plot)
%
% Last updated : 27/06/12

%set default variables
n1=100;
n2=100;
n3=100;
dX1=[];
dX2=[];
dX3=[];
c1=[];
c2=[];
c3=[];

%Collects all arguments
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end

minX1=min(X1(:));
minX2=min(X2(:));
minX3=min(X3(:));
maxX1=max(X1(:));
maxX2=max(X2(:));
maxX3=max(X3(:));

if isempty(c1) %coordinate 1
    if isempty(dX1) %divide hist to n bins
        c1=linspace(minX1,maxX1,n1);
    else %divide hist to sections of length dX
        c1=minX1:dX1:(maxX1+dX1);
        c1=x1i+dX1/2;
    end
end
if isempty(c2) %coordinate 1
    if isempty(dX2) %divide hist to n bins
        c2=linspace(minX2,maxX2,n2);
    else %divide hist to sections of length dX
        c2=minX2:dX2:(maxX2+dX2);
        c2=x2i+dX2/2;
    end
end
if isempty(c3) %coordinate 1
    if isempty(dX3) %divide hist to n bins
        c3=linspace(minX3,maxX3,n3);
    else %divide hist to sections of length dX
        c3=minX3:dX3:(maxX3+dX3);
        c3=x3i+dX3/2;
    end
end

nBin1=numel(c1);
nBin2=numel(c2);
nBin3=numel(c3);

%interpolate to grid
x1r=interp1(c1,1:numel(c1),X1(:),'nearest');
x2r=interp1(c2,1:numel(c2),X2(:),'nearest');
x3r=interp1(c3,1:numel(c3),X3(:),'nearest');

D=accumarray([x1r x2r x3r],1,[nBin1 nBin2 nBin3]);

if nargout==2 %plot results
    h=axes;
    cmap=flipud(colormap(gray(256)));
    colormap(cmap);
    plot3(c1,c2,c3,log10(D));
    cb=colorbar;
    ylabel(cb,'log_{10}#'); 
end



