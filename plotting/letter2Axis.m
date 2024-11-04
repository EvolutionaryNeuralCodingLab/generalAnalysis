function h=letter2Axis(ha,letter,shiftXY,corner,varargin)

if strcmp(corner,'TL')
    dx=0;
    dy=1;
elseif strcmp(corner,'TR')
    
elseif strcmp(corner,'BL')
    
elseif strcmp(corner,'BR')
    
end


xShift=(shiftXY(1)/ha.Position(3)+dx);%*diff(xlim);
yShift=(shiftXY(2)/ha.Position(4)+dy);%*diff(ylim);

h=text(xShift,yShift,letter,'Units','normalized',varargin{:});

