%  h=TransparentLine(X,Y,C,W,transparency,varargin)
% Function purpose : Plots transparant line by drawing a rectangular patch
%
% Function recives :   X - x coordinates
%                      Y- y coordinates
%                      C-color of line (either 'r' or [0.1 0.4 0.7])
%                      W - line width
%                      transparency - a number between [0 1] indicating the transparency of the line
%                      varargin - other arguments of the patch.
%
% Function give back : h - line patch handle
%
% Last updated : 07/12/09
function h=TransparentLine(X,Y,C,W,transparency,varargin)
m=(Y(2)-Y(1))/(X(2)-X(1));
if m==0
    PX(1)=X(1);
    PX(2)=X(1);
    PX(3)=X(2);
    PX(4)=X(2);
    PY(1)=Y(1)+W/2;
    PY(2)=Y(1)-W/2;
    PY(3)=Y(2)-W/2;
    PY(4)=Y(2)+W/2;
else
    if isinf(m)
        dX=W;
    else
        dX=W*cos(pi/2-atan(m));
    end
    a=-1/m;

    b1=-a*X(1)+Y(1);
    PX(1)=X(1)+dX/2;
    PX(2)=X(1)-dX/2;
    PY(1)=a*PX(1)+b1;
    PY(2)=a*PX(2)+b1;

    b2=-a*X(2)+Y(2);
    PX(4)=X(2)+dX/2;
    PX(3)=X(2)-dX/2;
    PY(3)=a*PX(3)+b2;
    PY(4)=a*PX(4)+b2;
    
    %hold on;plot(X,Y,'.r','markerSize',5);
    %plot(PX,PY,'.');
end


h=patch(PX,PY,C,'EdgeColor','none','FaceAlpha',transparency,varargin{:});


%line([XL(1) XL(2)],[i*aY+bY i*aY+bY],'LineWidth',2,'Color',[.1 .5 .7]);