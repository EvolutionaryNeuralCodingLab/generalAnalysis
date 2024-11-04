function horizontalLegend(hLegObj,varargin)
%Default variables
height=0.5; %height of text and symbols (0.5 - in the middle) 
width=0.2;%The width of every legend segments
gapMarker2Text=0.05; %Spacing between segments
gapText2NextMarker=0.1;
fontSize=9;

% Output list of default variables
% print out default arguments and values if no inputs are given
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

nLeg=numel(hLegObj)/3;
cumPos=0.05;
for i=1:nLeg
    pHLine=nLeg+(i-1)*2+1;
    %for line
    hLegObj(pHLine).XData=[cumPos cumPos+width];
    hLegObj(pHLine).YData=[height height];
    %for marker
    hLegObj(pHLine+1).XData=cumPos+width/2;
    hLegObj(pHLine+1).YData=height;
    
    %for text
    cumPos=cumPos+width+gapMarker2Text;
    hLegObj(i).Position=[cumPos height 0];
    hLegObj(i).HorizontalAlignment='left';
    hLegObj(i).VerticalAlignment='middle';
    hLegObj(i).FontSize=fontSize;
    cumPos=cumPos+hLegObj(i).Extent(3)+gapText2NextMarker;
end

%set everything back
shift=cumPos-1;
for i=1:nLeg
    pHLine=nLeg+(i-1)*2+1;
    %for line
    hLegObj(pHLine).XData=hLegObj(pHLine).XData-shift;
    %for marker
    hLegObj(pHLine+1).XData=hLegObj(pHLine+1).XData-shift;
    %for text
    hLegObj(i).Position=[hLegObj(i).Position(1)-shift hLegObj(i).Position(2:3)];
end
