function varargout = mriViewer(varargin)
% MRIVIEWER MATLAB code for mriViewer.fig
%      MRIVIEWER, by itself, creates a new MRIVIEWER or raises the existing
%      singleton*.
%
%      H = MRIVIEWER returns the handle to a new MRIVIEWER or the handle to
%      the existing singleton*.
%
%      MRIVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MRIVIEWER.M with the given input arguments.
%
%      MRIVIEWER('Property','Value',...) creates a new MRIVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mriViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mriViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mriViewer

% Last Modified by GUIDE v2.5 02-Jan-2013 18:07:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mriViewer_OpeningFcn, ...
    'gui_OutputFcn',  @mriViewer_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mriViewer is made visible.
function mriViewer_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
if nargin<4
    addpath('/media/sil1/Data/MRI/');
    %load('/media/sil1/Data/MRI/turtle03_5_ReducedOnlyBrain_mriViewerData.mat');
    load('/media/sil1/Data/MRI/turtle03_5_mriViewerData.mat');
    handles.data=turtleMRI;
    handles.surfPatches=surfPatches;
    handles.gui.pixelSize=0.039; %[mm]
    %{
    [FileName,PathName,FilterIndex] = uigetfile('*.mat','Please choose a matlab file with 3D data');
    S=load([PathName '\' FileName]);
    varName=fieldnames(S);
    handles.data=getfield(S,varName{1});
    clear S;
    handles.gui.pixelSize=1;
    %}
elseif nargin==6
    handles.data=double(varargin{1});
    handles.surfPatches=varargin{2}; %struct with faces and vertices
    handles.gui.pixelSize=varargin{3}; %uM units
elseif nargin==5
    handles.data=double(varargin{1});
    handles.gui.pixelSize=varargin{2}; %uM units
elseif nargin==4
    handles.data=double(varargin{1});
    handles.surfPatches.faces=[];
    handles.surfPatches.vertices=[];
else
    error('incorrect number of inputs! input should be: mriViewer(X), X beind a 3D matrix');
end
handles.gui.dataSiz=size(handles.data);

handles.gui.maxGridLength=ceil(sqrt(sum(handles.gui.dataSiz.^2))); %find maximal axis in image in pixels
handles.gui.section=1;
handles.gui.frameNumber=1;
set(handles.directionEdit,'String','[0 1 0]');
handles.gui.direction=[0 1 0];
handles.gui.angle=str2num(get(handles.angleEdit,'string'));
handles.gui.interpolationMethod=get(handles.interpolationMenu,'value');

set(handles.totalFrameNumber,'string',[' / ' num2str(handles.gui.dataSiz(handles.gui.section))]);
set(handles.frame,'String',num2str(handles.gui.frameNumber));

%plot Basic rectangle
%handles.PlaneAxis=copyobj(handles.rotationPlot,handles.figure1);

axes(handles.rotationPlot);
handles.gui.box=patch('Faces',[1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8],...
    'Vertices',bsxfun(@times,[0 0 0;0 1 0;1 1 0;1 0 0;0 0 1;0 1 1;1 1 1;1 0 1],handles.gui.dataSiz),...
    'FaceColor','r','FaceAlpha',0.05,'EdgeAlpha',0.1);

handles.surfacePatch=patch('Faces',handles.surfPatches.faces,'Vertices',handles.surfPatches.vertices);
set(handles.surfacePatch,'FaceColor',[1 0 0],'EdgeColor','none','FaceLighting','phong');

text(1,handles.gui.dataSiz(2)/2,handles.gui.dataSiz(3)/2,'C');
text(handles.gui.dataSiz(1),handles.gui.dataSiz(2)/2,handles.gui.dataSiz(3)/2,'R');
xlabel('X - rosto-caudal');
ylabel('Y - lateral');
zlabel('Z - ventro-dorsal');

light('Position',[1 1 0],'Style','infinite');
view(50,30);
axis equal tight;hold on;

%plot initial schematic plane
%axes(handles.PlaneAxis);
%set(handles.PlaneAxis,'Visible','off','NextPlot','add');

tmpVertValues=[1 0 0;1 1 0;1 1 1;1 0 1];
tmpPlaneValues=[handles.gui.frameNumber handles.gui.dataSiz(2) handles.gui.dataSiz(3)];
tmpVertices=bsxfun(@times,tmpVertValues,tmpPlaneValues);
handles.planePatch = patch('Faces',[1 2 3 4],'Vertices',tmpVertices,'FaceColor','b','FaceAlpha',0.3,'EdgeAlpha',0.1);

%view(50,30);
%axis equal tight vis3d;
%set(handles.PlaneAxis,'xLim',get(handles.rotationPlot,'xLim'));
%set(handles.PlaneAxis,'yLim',get(handles.rotationPlot,'yLim'));
%set(handles.PlaneAxis,'zLim',get(handles.rotationPlot,'zLim'));
%handles.rotation=rotate3d(handles.PlaneAxis);
%set(handles.rotation,'Enable','on');

handles=plotData(handles);
colormap(gray);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mriViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = mriViewer_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function handles=plotData(handles)
%update gui handles
histogramEqualization=true;

if get(handles.manualAngle,'value')
    handles.gui.pValid=find(handles.gui.Xp>=1 & handles.gui.Xp<=handles.gui.dataSiz(1)...
        & handles.gui.Yp>=1 & handles.gui.Yp<=handles.gui.dataSiz(2)...
        & handles.gui.Zp>=1 & handles.gui.Zp<=handles.gui.dataSiz(3));
    
    [I,J] = ind2sub(size(handles.gui.Xp),handles.gui.pValid);
    handles.gui.minValidI=min(I);
    handles.gui.minValidJ=min(J);
    handles.gui.maxValidI=max(I);
    handles.gui.maxValidJ=max(J);
    
    sliceImg=0./zeros(size(handles.gui.Xp));
    
    switch handles.gui.interpolationMethod
        case 1
            linearIdx=sub2ind(handles.gui.dataSiz,...
                round(handles.gui.Xp(handles.gui.pValid)),round(handles.gui.Yp(handles.gui.pValid)),round(handles.gui.Zp(handles.gui.pValid)));
            sliceImg(handles.gui.pValid) = handles.data(linearIdx);
        case 2
            sliceImg(handles.gui.pValid) = interp3(handles.data,handles.gui.Yp(handles.gui.pValid),handles.gui.Xp(handles.gui.pValid),handles.gui.Zp(handles.gui.pValid),'linear');
            %sliceImg(handles.gui.pValid) = interp3(handles.data,handles.gui.Xp(handles.gui.pValid),handles.gui.Yp(handles.gui.pValid),handles.gui.Zp(handles.gui.pValid),'linear');
        case 3
            sliceImg(handles.gui.pValid) = interp3(handles.data,handles.gui.Yp(handles.gui.pValid),handles.gui.Xp(handles.gui.pValid),handles.gui.Zp(handles.gui.pValid),'cubic');
    end
    %
    %sliceImg(pValid) = interp3(handles.data,Y(pValid),X(pValid),Z(pValid),'linear');

    sliceImg=sliceImg(handles.gui.minValidI:handles.gui.maxValidI,handles.gui.minValidJ:handles.gui.maxValidJ);
    
    %update mri slice plot
    axes(handles.imagePlot);
    imagesc(flipud(sliceImg));
    if histogramEqualization
        I=adapthisteq(I,'clipLimit',0.03);
    end

    axis image;
    disp('Finished calculating interpolation');
else
    axes(handles.imagePlot);
    switch handles.gui.section
        case 1
            net = denoisingNetwork('DnCNN');
            I=denoiseImage(rot90(squeeze(handles.data(handles.gui.frameNumber,:,:))),net);
            if histogramEqualization
                I=adapthisteq(I,'clipLimit',0.03);
            end
            handles.gui.p1=imagesc(I);
            %handles.gui.p1=imagesc(rot90(squeeze(handles.data(handles.gui.frameNumber,:,:))),[0 10000]);
            %colormap(gca,'jet');
            %colorbar;
        case 2
            handles.gui.p1=imagesc(rot90(squeeze(handles.data(:,handles.gui.frameNumber,:))));
        case 3
            handles.gui.p1=imagesc(squeeze(handles.data(:,:,handles.gui.frameNumber))');
    end
    axis image;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Callback functions %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function slider1_Callback(hObject, eventdata, handles)
tmpFrameNumber=1+ceil(get(handles.slider1,'Value')*(handles.gui.dataSiz(handles.gui.section)-1));
tmpFrameDiff=tmpFrameNumber-handles.gui.frameNumber;
handles.gui.frameNumber=tmpFrameNumber;

tmp=zeros(1,3);
tmp(handles.gui.section)=1;
tmpVertices=get(handles.planePatch,'Vertices')+tmpFrameDiff*ones(4,1)*tmp;
set(handles.planePatch,'Vertices',tmpVertices);

set(handles.frame,'String',num2str(handles.gui.frameNumber));
set(handles.slider1,'Value',(handles.gui.frameNumber-1)/(handles.gui.dataSiz(handles.gui.section)-1));

if get(handles.manualAngle,'value')
    switch handles.gui.section
        case 1
            handles.gui.Xp=handles.gui.Xp+tmpFrameDiff;
        case 2
            handles.gui.Yp=handles.gui.Yp+tmpFrameDiff;
        case 3
            handles.gui.Zp=handles.gui.Zp+tmpFrameDiff;
    end
end

handles=plotData(handles);
guidata(hObject, handles);

function Section_Callback(hObject, eventdata, handles)
handles.gui.section=get(handles.Section,'value');
set(handles.totalFrameNumber,'string',[' / ' num2str(handles.gui.dataSiz(handles.gui.section))]);
if get(handles.manualAngle,'value')
    handles.gui.frameNumber=round(handles.gui.origin(handles.gui.section));
    set(handles.frame,'String',num2str(handles.gui.frameNumber));
    set(handles.slider1,'Value',(handles.gui.frameNumber-1)/(handles.gui.dataSiz(handles.gui.section)-1));
    %******* For accuracy the planes (schematic+Xp,Yp,Zp) should be updated due to the slight movement due to frame rounding
else
    handles.gui.frameNumber=1;
    set(handles.frame,'String',num2str(handles.gui.frameNumber));
    set(handles.slider1,'Value',(handles.gui.frameNumber-1)/(handles.gui.dataSiz(handles.gui.section)-1));
    
    switch handles.gui.section
        case 1
            tmpVertValues=[1 0 0;1 1 0;1 1 1;1 0 1];
            tmpPlaneValues=[handles.gui.frameNumber handles.gui.dataSiz(2) handles.gui.dataSiz(3)];
            tmpVertices=bsxfun(@times,tmpVertValues,tmpPlaneValues);
            set(handles.planePatch,'Vertices',tmpVertices);
        case 2
            tmpVertValues=[0 1 0;1 1 0;1 1 1;0 1 1];
            tmpPlaneValues=[handles.gui.dataSiz(1) handles.gui.frameNumber handles.gui.dataSiz(3)];
            tmpVertices=bsxfun(@times,tmpVertValues,tmpPlaneValues);
            set(handles.planePatch,'Vertices',tmpVertices);
        case 3
            tmpVertValues=[0 0 1;1 0 1;1 1 1;0 1 1];
            tmpPlaneValues=[handles.gui.dataSiz(1) handles.gui.dataSiz(2) handles.gui.frameNumber];
            tmpVertices=bsxfun(@times,tmpVertValues,tmpPlaneValues);
            set(handles.planePatch,'Vertices',tmpVertices);
    end
end
handles=plotData(handles);
guidata(hObject, handles);

function backward_Callback(hObject, eventdata, handles)
if handles.gui.frameNumber<=1
    disp('First slice');
else
    handles.gui.frameNumber=handles.gui.frameNumber-1;
    set(handles.slider1,'Value',(handles.gui.frameNumber-1)/(handles.gui.dataSiz(handles.gui.section)-1));
    set(handles.frame,'string',num2str(handles.gui.frameNumber));
    
    %update schematic plot
    tmp=zeros(1,3);
    tmp(handles.gui.section)=1;
    tmpVertices=get(handles.planePatch,'Vertices')-ones(4,1)*tmp;
    set(handles.planePatch,'Vertices',tmpVertices);
    
    if get(handles.manualAngle,'value')
        switch handles.gui.section
            case 1
                handles.gui.Xp=handles.gui.Xp-1;
            case 2
                handles.gui.Yp=handles.gui.Yp-1;
            case 3
                handles.gui.Zp=handles.gui.Zp-1;
        end
    end
    
    handles=plotData(handles);
    guidata(hObject, handles);
end

function forward_Callback(hObject, eventdata, handles)
if handles.gui.frameNumber>=handles.gui.dataSiz(handles.gui.section)
    disp('Last slice');
else
    handles.gui.frameNumber=handles.gui.frameNumber+1;
    set(handles.slider1,'Value',(handles.gui.frameNumber-1)/(handles.gui.dataSiz(handles.gui.section)-1));
    set(handles.frame,'string',num2str(handles.gui.frameNumber));
    
    tmp=zeros(1,3);
    tmp(handles.gui.section)=1;
    tmpVertices=get(handles.planePatch,'Vertices')+ones(4,1)*tmp;
    set(handles.planePatch,'Vertices',tmpVertices);
    
    if get(handles.manualAngle,'value')
        switch handles.gui.section
            case 1
                handles.gui.Xp=handles.gui.Xp+1;
            case 2
                handles.gui.Yp=handles.gui.Yp+1;
            case 3
                handles.gui.Zp=handles.gui.Zp+1;
        end
    end
    
    handles=plotData(handles);
    guidata(hObject, handles);
end

function frame_Callback(hObject, eventdata, handles)
inputFrame=str2num(get(handles.frame,'String'));
if inputFrame>handles.gui.dataSiz(handles.gui.section)
    disp('The chosen slide is larger than stack size');
    set(handles.frame,'string',num2str(handles.gui.frameNumber));
elseif inputFrame<1
    disp('The smallest slide number is 1');
    set(handles.frame,'string',num2str(handles.gui.frameNumber));
else
    tmpFrameDiff=inputFrame-handles.gui.frameNumber;
    handles.gui.frameNumber=inputFrame;
    set(handles.slider1,'Value',(handles.gui.frameNumber-1)/(handles.gui.dataSiz(handles.gui.section)-1));
    
    tmp=zeros(1,3);
    tmp(handles.gui.section)=1;
    tmpVertices=get(handles.planePatch,'Vertices')+tmpFrameDiff*ones(4,1)*tmp;
    set(handles.planePatch,'Vertices',tmpVertices);
    
    if get(handles.manualAngle,'value')
        switch handles.gui.section
            case 1
                handles.gui.Xp=handles.gui.Xp+tmpFrameDiff;
            case 2
                handles.gui.Yp=handles.gui.Yp+tmpFrameDiff;
            case 3
                handles.gui.Zp=handles.gui.Zp+tmpFrameDiff;
        end
    end
    
    handles=plotData(handles);
    guidata(hObject, handles);
end

function Send2Fig_Callback(hObject, eventdata, handles)
figureHandle=figure;
copyobj(handles.imagePlot,figureHandle);
colormap(gray);
assignin('base','figureHandleMRI',figureHandle);

function planeFigure_Callback(hObject, eventdata, handles)
figureHandle=figure;
copyobj(handles.rotationPlot,figureHandle);
assignin('base','figureHandlePlane',figureHandle);

function angleEdit_Callback(hObject, eventdata, handles)
handles.gui.angle=str2num(get(handles.angleEdit,'string'));
guidata(hObject, handles);

function directionEdit_Callback(hObject, eventdata, handles)
eval(['handles.gui.direction=' get(handles.directionEdit,'string') ';']);
if numel(handles.gui.direction)~=3
    disp('Incorrect input to direction vector ([x y z]');
    return;
end
guidata(hObject, handles);

function manualAngle_Callback(hObject, eventdata, handles)
if get(handles.manualAngle,'value')
    %prepare pixels in initial plane
    switch handles.gui.section
        case 1 %X
            handles.gui.origin=[handles.gui.frameNumber handles.gui.dataSiz(2)/2 handles.gui.dataSiz(3)/2];
            [handles.gui.Yp,handles.gui.Zp]=meshgrid(handles.gui.origin(2)+(-handles.gui.maxGridLength:handles.gui.maxGridLength),handles.gui.origin(3)+(-handles.gui.maxGridLength:handles.gui.maxGridLength)); %prepare X,Y grid
            handles.gui.Xp=handles.gui.frameNumber*ones(size(handles.gui.Yp)); %prepare X grid
        case 2 %Y
            handles.gui.origin=[handles.gui.dataSiz(1)/2 handles.gui.frameNumber handles.gui.dataSiz(3)/2];
            [handles.gui.Xp,handles.gui.Zp]=meshgrid(-handles.gui.maxGridLength+handles.gui.origin(1):handles.gui.maxGridLength+handles.gui.origin(1),-handles.gui.maxGridLength+handles.gui.origin(3):handles.gui.maxGridLength+handles.gui.origin(3)); %prepare X,Y grid
            handles.gui.Yp=handles.gui.frameNumber*ones(size(handles.gui.Xp)); %prepare Y grid
        case 3 %Z
            handles.gui.origin=[handles.gui.frameNumber handles.gui.dataSiz(2)/2 handles.gui.dataSiz(3)/2];
            [handles.gui.Xp,handles.gui.Yp]=meshgrid(-handles.gui.maxGridLength+handles.gui.origin(1):handles.gui.maxGridLength+handles.gui.origin(1),-handles.gui.maxGridLength+handles.gui.origin(2):handles.gui.maxGridLength+handles.gui.origin(2)); %prepare X,Y grid
            handles.gui.Zp=handles.gui.frameNumber*ones(size(handles.gui.Xp)); %prepare Z grid
    end
else
    %************* Insert this into a subfunction and
    switch handles.gui.section
        case 1
            tmpVertValues=[1 0 0;1 1 0;1 1 1;1 0 1];
            tmpPlaneValues=[handles.gui.frameNumber handles.gui.dataSiz(2) handles.gui.dataSiz(3)];
            tmpVertices=bsxfun(@times,tmpVertValues,tmpPlaneValues);
            set(handles.planePatch,'Vertices',tmpVertices);
        case 2
            tmpVertValues=[0 1 0;1 1 0;1 1 1;0 1 1];
            tmpPlaneValues=[handles.gui.dataSiz(1) handles.gui.frameNumber handles.gui.dataSiz(3)];
            tmpVertices=bsxfun(@times,tmpVertValues,tmpPlaneValues);
            set(handles.planePatch,'Vertices',tmpVertices);
        case 3
            tmpVertValues=[0 0 1;1 0 1;1 1 1;0 1 1];
            tmpPlaneValues=[handles.gui.dataSiz(1) handles.gui.dataSiz(2) handles.gui.frameNumber];
            tmpVertices=bsxfun(@times,tmpVertValues,tmpPlaneValues);
            set(handles.planePatch,'Vertices',tmpVertices);
    end
    handles=plotData(handles);
end
guidata(hObject, handles);

function rotateImg_Callback(hObject, eventdata, handles)
if ~get(handles.manualAngle,'value')
    disp('Rotation is not enabled when manual angle is set to off');
else
    %rotate schematics
    rotate(handles.planePatch,handles.gui.direction,handles.gui.angle,handles.gui.origin);
    
    %calculate rotation matrix
    R = rotationmat3D(handles.gui.angle/180*pi,handles.gui.direction); %calculate the 3D rotation matrix
    
    %rotate image
    tmp=R*([handles.gui.Xp(:)-handles.gui.origin(1),handles.gui.Yp(:)-handles.gui.origin(2),handles.gui.Zp(:)-handles.gui.origin(3)]'); %calculate rotation on 3 coordinate vectors
    handles.gui.Xp(:)=tmp(1,:)+handles.gui.origin(1); %return to matrix presentation
    handles.gui.Yp(:)=tmp(2,:)+handles.gui.origin(2); %return to matrix presentation
    handles.gui.Zp(:)=tmp(3,:)+handles.gui.origin(3); %return to matrix presentation
end

handles=plotData(handles);
guidata(hObject, handles);

function interpolationMenu_Callback(hObject, eventdata, handles)
handles.gui.interpolationMethod=get(handles.interpolationMenu,'value');
handles=plotData(handles);
guidata(hObject, handles);
%1=nearest
%2=Linear
%3=Cubic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Create functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function directionEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function angleEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Section_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function frame_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function interpolationMenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
