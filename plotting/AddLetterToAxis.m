function [Txth]=AddLetterToAxis(h,Letter,LetColor,LetSize,LetWeight,Pos,fixed,PosOutside)
% [Txth]=AddLetterToAxis(h,Letter,LetColor,LetSize,LetWeight,Pos)
% Function purpose : Adds a letter on the top left corner of an axis
%
% Recives : h - axes handle
%                       Letter - the letter of the axis (string)
%                       LetColor - the color of the letter, e.g. 'k','w','r',[1 1 0] (default='k') 
%                       LetSize - the size (default=14)
%                       LetWeight - the weight of the letter (default='Bold')
%                       Pos - 'Left'/'Right'/'LeftTop'/'RightTop' (default='Right')
%                       PosOutside - [0/1]
%                       fixed [0/1]- forces fixed distance from axes (default ==0)
%                       boundries - good for size varying multiple axes
%                       plots...
%                                                                                                                                                        
% Function give back : Txth - a handle for the letter
% Recomended usage  : AddLetterToAxis(h,'A')
% Last updated : 10/01/11
if nargin<6
    PosOutside=false;
end
if nargin<5
    LetWeight='Bold';
end
if nargin<4
    LetSize=14;
end
if nargin<3
    LetColor='k';
end
if nargin<6
    Pos='Left';
end
if nargin<7
    fixed=0;
end
axes(h);
if strcmpi(Pos,'Left') || strcmpi(Pos,'LeftTop')
    if ~fixed
        if strcmpi(Pos,'Left')
            Txth=text('Position',[0.03 0.99],'string',Letter,'FontSize',LetSize,'FontWeight',LetWeight,'color',LetColor,'Units','normalized','VerticalAlignment','top');
        else
            Txth=text('Position',[0.03 1.14],'string',Letter,'FontSize',LetSize,'FontWeight',LetWeight,'color',LetColor,'Units','normalized','VerticalAlignment','top');
        end
    else
        axUnits=get(h,'Units');
        set(h,'Units','Pixels');
        Txth=annotation('textbox',[0 0 0.01 0.01],'string',Letter,'FontSize',LetSize,'FontWeight',LetWeight,'color',LetColor,'Units','Pixels','VerticalAlignment','top',...
            'LineStyle','None');
        set(Txth,'FitBoxToText','On')
        TxPos=get(Txth,'Pos');
        axPos=get(h,'Pos');
        set(Txth,'Pos',[axPos(1)+5 axPos(2)+axPos(4)-TxPos(4) TxPos(3) TxPos(4)])
        set(h,'Units',axUnits)
        set(Txth,'Units','Normalized')
    end
else
    if ~fixed
        if strcmpi(Pos,'Right')
            Txth=text('Position',[0.93 0.99],'string',Letter,'FontSize',LetSize,'FontWeight',LetWeight,'color',LetColor,'Units','normalized','VerticalAlignment','top');
        else
            Txth=text('Position',[0.93 1.14],'string',Letter,'FontSize',LetSize,'FontWeight',LetWeight,'color',LetColor,'Units','normalized','VerticalAlignment','top');
        end
        
    else
        axUnits=get(h,'Units');
        set(h,'Units','Pixels');
        Txth=annotation('textbox',[0 0 0.01 0.01],'string',Letter,'FontSize',LetSize,'FontWeight',LetWeight,'color',LetColor,'Units','Pixels','VerticalAlignment','Top',...
                                                        'HorizontalAlignment','Right','LineStyle','None');
        set(Txth,'FitBoxToText','On')
        TxPos=get(Txth,'Pos');
        axPos=get(h,'Pos');
        set(Txth,'Pos',[axPos(1)+axPos(3)-TxPos(3) axPos(2)+axPos(4)-TxPos(4) TxPos(3) TxPos(4)])
        set(h,'Units',axUnits)
        set(Txth,'Units','Normalized')

    end
end
