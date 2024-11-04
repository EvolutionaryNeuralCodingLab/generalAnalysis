% [Fig,ax,Pos]=FigureAxisGenerator(FigWidth,FigHight,m,n,p,WidthStrech,HightStrech,WidthOffset,HightOffset,DefaultFontSize);
% Function purpose : Makes a figure with axis
%
% Function recives :    FigWidth - width of figure [cm]
%                       FigHight - hight of figure [cm]
%                       m,n,p (as in subplot) - breaks the figure window into an m-by-n matrix of small axes.
%                       p is a cell vector, specifies an axes position that covers all the subplot positions 
%                       listed in P, including those spanned by P for every axis
%                       WidthStrech [0-1] - the streach in axis width (if not entered == 0)
%                       HightStrech [0-1] - the streach in axis hight (if not entered == 0)
%                       WidthOffset [0-1] - the x offset of all axis (if not entered == 0)
%                       HightOffset [0-1] - the y offset of all axis (if not entered == 0)
%                       DefaultFontSize [Pt] - the default font size for all axes (if not entered == 10)
%
% Function give back :  Fig - the figure handle
%                       ax - the axis handles
%                       Pos - the position of the axis [left bottom width height]
%
% Last updated : 28/05/12
function [Fig,ax,Pos]=FigureAxisGenerator(FigWidth,FigHight,m,n,p,WidthStrech,HightStrech,WidthOffset,HightOffset,DefaultFontSize)

if ~exist('WidthStrech','var')
    WidthStrech=0;
end
if ~exist('HightStrech','var')
    HightStrech=0;
end
if ~exist('DefaultFontSize','var')
    DefaultFontSize=10;
end

%Max value for width is 30 and for hight is 15
FigRect = [1,1,FigWidth,FigHight];
Fig=figure;
set(Fig,'units','centimeters','position',FigRect,'resize','off','Color','w');

NAxes=length(p);
for i=1:NAxes
    ax(i)=subplot(m,n,p{i},'replace');
    Pos(i,:)=get(ax(i),'Position');
    OuterPos(i,:)=get(ax(i),'OuterPosition');
    %Ti(i,:)=get(ax(i),'TightInset'); %[left bottom right top]
    
    Pos(i,1)=Pos(i,1)-(Pos(i,1)-OuterPos(i,1))*WidthStrech;
    Pos(i,3)=Pos(i,3)-2*(Pos(i,3)-OuterPos(i,3))*WidthStrech;
    Pos(i,2)=Pos(i,2)-(Pos(i,2)-OuterPos(i,2))*HightStrech;
    Pos(i,4)=Pos(i,4)-2*(Pos(i,4)-OuterPos(i,4))*HightStrech;
    
    if exist('WidthOffset','var')
        Pos(i,1)=Pos(i,1)+WidthOffset;
    end
    if exist('HightOffset','var')
        Pos(i,2)=Pos(i,2)+HightOffset;
    end
    
    set(gca,'Position',Pos(i,:));
    set(gca,'FontSize',DefaultFontSize);
end
set(gcf,'paperunits','centimeters','papertype','A4');
set(gcf,'paperpositionmode','auto');

