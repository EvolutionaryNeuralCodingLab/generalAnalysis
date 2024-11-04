function [fps]=InterpActivityMovieMaker(MovieName,ActM,Channels,En,TimeSmoothing,TimeFactor,TimeBin,SF,OL,RImage)
% [fps]=InterpActivityMovieMaker(MovieName,ActM,Channels,ChannelMap,TimeSmoothing,TimeFactor,TimeBin,SF,OL);
% Function purpose : Makes smoothed (in time and space) avi movies from 2D activity arrays.
%
% Recives :    MovieName - the name of the movie
%              ActM [NxM] - activity
%              Channels - the channels names of every row in ActM 
%              En - the spatial map of the electrodes 
%              TimeSmoothingStd - the temporal width of the gaussian for time smoothing (in units of ActM)          
%              TimeFactor - the time factor used for slowing or increasing the movie speed, 
%                           if 0.1 is chosen the movie will be 10 times slower than in reality, if 5 is chosen the movie will be 5 times faster than in reality
%              TimeBin - the real time of every colomn in ActM [sec]
%              SF - spatial factor - how spatially big (how many pixels) is the signal of one channel  
%              OL - the overlap between the activity of adjucent channels
%               RImage - the real network image, image should be adjusted to fit the electrode grid                                                                                                                                                        
%
% Function give back :  
% Recomended usage  : 
% Last updated : 

maxFramesPerFile=3000;
backgroundColor='w';

% Apply a low pass filter over frames.
[NChannels Nframes]=size(ActM);
TimeKern=fspecial('gaussian',TimeSmoothing*4, TimeSmoothing);
for i=1:NChannels
    ActM(i,:)=convn(ActM(i,:),TimeKern,'same');
end

%creating translation between electrodes and locations
translation=[];
En=fliplr(En);
for i=1:NChannels
    [n,m]=find(En==Channels(i));
    ChannelMap(1:2,i)=[n;m];
    XLocations(i)=n*SF-ceil(SF/2);
    YLocations(i)=m*SF-ceil(SF/2);
end

%Calculating spatial spread
SpaceKern = fspecial('gaussian', (OL*SF), (OL*SF)/6);
Mx=max(ChannelMap(1,:))*SF;My=max(ChannelMap(2,:))*SF;
M=zeros([Mx My]);
M_XYInd = sub2ind(size(M),XLocations,YLocations);
M(M_XYInd)=ActM(:,1);

%Setting figure properties and calculating first frame
O1 = imfilter(M,SpaceKern);
f=figure('Position',[20 70 910 730]);
h=axes;
if exist('RImage','var'); %Add Real Image
    hS=surf(O1,RImage,'FaceColor','texturemap','EdgeColor','none', 'CDataMapping','direct');
    %view(-17,52);%view(116,59);
    %set(gca,'CameraPosition',1e3*[1.225 -0.0422 0.262]);
    %set(gca,'CameraPosition',1e3*[1.175 -0.180 0.360]);
else
    hS=surf(h,O1,'EdgeColor','flat');
end

axis tight;axis off;
set(h,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
set(h,'XTick',[],'YTick',[],'ZTick',[]);
whitebg(f,backgroundColor);
set(f,'color',backgroundColor)

ZTop=max(imfilter(max(ActM(:)),SpaceKern));
ZBottom=min(imfilter(min(ActM(:)),SpaceKern));

zlim([ZBottom ZTop]);
set(h,'CLim',[ZBottom ZTop]/2,'CLimMode','Manual');
set(h,'nextplot','replacechildren');
set(h,'position',get(h,'OuterPosition'));
set(f,'Renderer','zbuffer');

view(-20,70);
%view(0,90);
fps=round(1/TimeBin*TimeFactor);

writerObj = VideoWriter([MovieName '.avi']);
writerObj.FrameRate=fps;
open(writerObj);

%aviobj = avifile([MovieName '.avi'],'fps',fps,'compression','none','keyframe',ceil(fps/30),'quality',75);
delete(hS);
% Record the movie
partCounter=1;
for i = 1:Nframes
    M=zeros([Mx My]);
    M(M_XYInd)=ActM(:,i);
    O = imfilter(M,SpaceKern);
    if exist('RImage','var'); %Add Real Image
        hS=surface(h,O,RImage,'FaceColor','texturemap','EdgeColor','none', 'CDataMapping','direct');
    else
        hS=surf(h,O,'EdgeColor','flat','EdgeColor','none');
        view(-20,70);
        %view(0,90);
        shading interp;
    end
    set(h,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
	%surf(O,O1,'EdgeColor','flat');
    
    F = getframe;
    writeVideo(writerObj,F);
    
    if round(i/10)==i/10
        disp(['Frame: ' num2str(i) ' / ' num2str(Nframes)]);
    end
    delete(hS);
end
close(writerObj);

fprintf('\nFinished\n');

