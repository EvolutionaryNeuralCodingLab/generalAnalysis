%  h=TransparentCircledLine(X,Y,C,W,transparency,varargin)
% Function purpose : Plots transparant line by drawing a rectangular patch
%
% Function recives :   X - x coordinates
%                      Y- y coordinatesnsity of the N connecting lines
%                      C-color of line (either 'r' or [0.1 0.4 0.7])
%                      W - the width of the line
%                      transparency - a number between [0 1] indicating the transparency of the line
%                      varargin - other arguments of the patch.
%
% Function give back : h - line patch handle
%
% Last updated : 07/12/09
function h=TransparentCircledLine(X,Y,C,W,transparency,varargin)
flipAngle=0;
if Y(2)<=Y(1)
    Y=[Y(2) Y(1)];
    X=[X(2) X(1)];
end
if X(2)<X(1)
    flipAngle=1;
end
m=(Y(2)-Y(1))/(X(2)-X(1));
if m==0
    PX(1)=X(1);
    PY(1)=Y(1)+W/2;

    angle=(pi/2+(pi/9)):(pi/9):(3*pi/2-pi/9);
    PX(2:9)=X(1)+(W/2)*cos(angle);
    PY(2:9)=Y(1)+(W/2)*sin(angle);

    PX(10)=X(1);
    PY(10)=Y(1)-W/2;

    PX(11)=X(2);
    PY(11)=Y(2)-W/2;

    angle=(3*pi/2+pi/9):(pi/9):(5*pi/2-pi/9);
    PX(12:19)=X(2)+(W/2)*cos(angle);
    PY(12:19)=Y(2)+(W/2)*sin(angle);

    PX(20)=X(2);
    PY(20)=Y(2)+W/2;
else
    if isinf(m)
        dX=W;
    else
        dX=W*cos(pi/2-atan(m));
    end
    a=-1/m;

    b1=-a*X(1)+Y(1);

    PX(1)=X(1)-dX/2;
    PY(1)=a*PX(1)+b1;

    angle=(atan(a)+pi+pi/9):(pi/9):(atan(a)+2*pi-pi/9);
    if flipAngle
        angle=fliplr(angle);
    end
    PX(2:9)=X(1)+(W/2)*cos(angle);
    PY(2:9)=Y(1)+(W/2)*sin(angle);

    PX(10)=X(1)+dX/2;
    PY(10)=a*PX(10)+b1;

    b2=-a*X(2)+Y(2);

    PX(11)=X(2)+dX/2;
    PY(11)=a*PX(11)+b2;

    angle=(atan(a)+(pi/9)):(pi/9):(atan(a)+pi-pi/9);
    if flipAngle
        angle=fliplr(angle);
    end
    PX(12:19)=X(2)+(W/2)*cos(angle);
    PY(12:19)=Y(2)+(W/2)*sin(angle);

    PX(20)=X(2)-dX/2;
    PY(20)=a*PX(20)+b2;

end


h=patch(PX,PY,C,'EdgeColor','none','FaceAlpha',transparency,varargin{:});


%line([XL(1) XL(2)],[i*aY+bY i*aY+bY],'LineWidth',2,'Color',[.1 .5 .7]);