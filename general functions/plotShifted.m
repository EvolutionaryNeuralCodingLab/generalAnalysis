function varargout=plotShifted(varargin)
%plotShifted   Linear shifted plot. 
%   plotShifted(X,Y) plots matrix Y versus vector X, where each column of Y is one linear plot and these plots are shifted (on the y axis) by a constant.
%   plot shifted can get the same inputs as Matlab PLOT function   
%
%   plotShifted(Y) plots the columns of Y versus their index.
%   If Y is complex, plotShifted(Y) is equivalent to plotShifted(real(Y),imag(Y)).
%   In all other uses of plotShifted, the imaginary part is ignored.
%
%   Various line types, plot symbols and colors may be obtained with
%   plotShifted(X,Y,S) where S is a character string made from one element
%   from any or all the following 3 columns:
%
%          b     blue          .     point              -     solid
%          g     green         o     circle             :     dotted
%          r     red           x     x-mark             -.    dashdot 
%          c     cyan          +     plus               --    dashed   
%          m     magenta       *     star             (none)  no line
%          y     yellow        s     square
%          k     black         d     diamond
%          w     white         v     triangle (down)
%                              ^     triangle (up)
%                              <     triangle (left)
%                              >     triangle (right)
%                              p     pentagram
%                              h     hexagram
%                         
%   For example, plotShifted(X,Y,'c+:') plots a cyan dotted line with a plus 
%   at each data point; plotShifted(X,Y,'bd') plots blue diamond at each data 
%   point but does not draw any line.
%
%   The plotShifted command, if no color is specified, makes automatic use of
%   the colors specified by the axes ColorOrder property.  By default,
%   plotShifted cycles through the colors in the ColorOrder property.  For
%   monochrome systems, plotShifted cycles over the axes LineStyleOrder property.
%
%   Note that RGB colors in the ColorOrder property may differ from
%   similarly-named colors in the (X,Y,S) triples.  For example, the 
%   second axes ColorOrder property is medium green with RGB [0 .5 0],
%   while plotShifted(X,Y,'g') plots a green line with RGB [0 1 0].
%
%   If you do not specify a marker type, plotShifted uses no marker. 
%   If you do not specify a line style, plotShifted uses a solid line.
%
%   plotShifted(AX,...) plots into the axes with handle AX.
%
%   plotShifted returns a column vector of handles to lineseries objects, one
%   handle per plotted line. 
%
%   The X,Y pairs, or X,Y,S triples, can be followed by 
%   parameter/value pairs to specify additional properties 
%   of the lines. For example, plotShifted(X,Y,'LineWidth',2,'Color',[.6 0 0]) 
%   will create a plot with a dark red line width of 2 points.
%
%   Example
%      x = 1:100
%      y = randn(100,1000);
%      plotShifted(x,y,'--rs','verticalShift',30,...
%                      'MarkerEdgeColor','k',...
%                      'MarkerFaceColor','g',...
%                      'MarkerSize',10)
%
%   See also PLOTTOOLS, SEMILOGX, SEMILOGY, LOGLOG, PLOTYY, PLOT3, GRID,
%   TITLE, XLABEL, YLABEL, AXIS, AXES, HOLD, LEGEND, SUBPLOT, SCATTER.

%   If the NextPlot axes property is "replace" (HOLD is off), plotShifted resets 
%   all axes properties, except Position, to their default values,
%   deletes all axes children (line, patch, text, surface, and
%   image objects), and sets the View property to [0 90].

%   Up dated by Mark Shein-Idelson 
%   $Revision: 5.19.4.10 $  $Date: 2014/05/14 $

%find the data input variable
pNumeric=cellfun(@(x) isnumeric(x),varargin);
pFirst=find(pNumeric,1,'first');
if nargin>1 && pNumeric(pFirst+1) %check if format is time,data or data
    pData=pFirst+1;%for time,data
else
    pData=pFirst;% for data
end

%look for inputs into varargin which are within the local list, place them into Par and remove them from varargin
localPropertyList={'verticalShift'};
for i=1:numel(localPropertyList)
    pProp=find(strcmp(localPropertyList{i},varargin));
    if ~isempty(pProp) %assign input value from varargin
        Par.(localPropertyList{i})=varargin{pProp+1};
        varargin(pProp:pProp+1)=[];
    else %assign empty value
        Par.(localPropertyList{i})=[]; 
    end
end

%determine the vertical shifts
[nSamples,nCh]=size(varargin{pData});
if isempty(Par.verticalShift); %automatically detection
    Par.verticalShift=cumsum(nanstd(varargin{pData})*3);
elseif numel(Par.verticalShift)==1 %manual detection
    Par.verticalShift=0:Par.verticalShift:(Par.verticalShift*(nCh-1));
end

varargin{pData}=bsxfun(@plus,varargin{pData},Par.verticalShift);
[varargout]=plot(varargin{:});
varargout={varargout(1:nargout)};
