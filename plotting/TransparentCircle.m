%  h=TransparentCircle(X,Y,R,C,transparency,NPoints,varargin)
% Function purpose : Plots transparant circle by drawing a rectangular patch
%
% Function recives :   X - x coordinate of circle center
%                      Y - y coordinate of circle center
%                      R - Radius of circle (R can be entered as a two coordinates to [R(1) R(2)] for aspect ratio different than 1:1.
%                      C-color of line ([0.1 0.4 0.7])
%                      transparency - a number between [0 1] indicating the transparency of the line
%                      NPoints - the number of points in the circle(default = 36)
%                      varargin - other arguments of the patch.
%
% Function give back :  h - circle patch handle
%
% Last updated : 07/12/09
function h=TransparentCircle(X,Y,R,C,transparency,NPoints,h,varargin)
normToAxis=true;
if nargin<6
    NPoints=36;
    h=gca;
end
if nargin<7
    h=gca;
end
if length(R)==1 %in case only one R is entered
    R(2)=R(1);
end
if normToAxis
    xl=diff(xlim(h));
    yl=diff(ylim(h));
    %R(1)=xl*R(1);
    %R(2)=yl*R(2);
end

angle=0:(360/NPoints):360;
PX=bsxfun(@plus,X,R(1).*cos(angle.*pi/180)');
PY=bsxfun(@plus,Y,R(2).*sin(angle.*pi/180)');
h=patch(PX,PY,C,'Parent',h,'EdgeColor','none','FaceAlpha',transparency,'FaceColor','flat','FaceVertexCData',C,varargin{:});