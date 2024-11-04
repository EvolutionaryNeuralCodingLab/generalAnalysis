classdef imageRegistrationGUI < handle
    
    properties (SetAccess=public)
        h
        par
    end
    
    properties (SetObservable, AbortSet = true, SetAccess=public)
    end
    
    properties (Constant)
    end
    
    methods
        function obj=imageRegistrationGUI() %class constructor
            %define initial conditions
            obj.h.img={[],[]};
            obj.par.img={[],[]};
            obj.h.layoutPlot={cell(1,3),cell(1,3)}';
            obj.par.allignMarks{1}=[];
            obj.par.allignMarks{2}=[];
            obj.h.regAxes=[];
            obj.par.imageType={[],[]};
            obj.par.pointCounter = [0 0];
            obj.par.multiImageFile = [false false];
            
            obj.par.trasformNames={'similarity','nonreflective similarity','affine','projective','polynomial','piecewise linear','lwm'};
            obj.par.trasformMinPoints=[3 2 4 4 10 4 6];
            obj.par.transformStrings={'similarity (3)','nonref. similarity (2)','affine (4)','projective (4)','3rd deg. polynomial (10)','piecewise linear (4)','lwm (6,better=12)'};
            obj.par.overlayMethodNames={'blend','falsecolor','sourceDominance','transDominance'};
            obj.par.overlayMethodString={'blend','falsecolor','sourceDominance','transDominance'};
            obj.par.transformsInList=[];
            obj.createGUI;
            
            obj.par.currentTransformName=obj.par.trasformNames{1};
            obj.par.currentTransformMinPoints=obj.par.trasformMinPoints(1);
        end
        
        function createGUI(obj)
            
            obj.par.screenPositionsMatlab=get(0,'MonitorPositions');
            obj.par.GUIPosition=[obj.par.screenPositionsMatlab(1,[1 2])+100 obj.par.screenPositionsMatlab(1,[4 4])-200];
            
            obj.h.mainFigure=figure('Position',obj.par.GUIPosition,'Name','image registration GUI',...
                'NumberTitle','off',  'HandleVisibility','off','Toolbar','figure','MenuBar','none','Visible','off');%,'CloseRequestFcn',@closeMainGUIFigure);
            
            %keep only a few toolbar objects
            figObj=findall(obj.h.mainFigure);
            for i=1:numel(figObj)
                if strcmp(figObj(i).Parent.Tag,'FigureToolBar') & ~any(strcmp(figObj(i).TooltipString,{'Zoom Out','Zoom In','Pan'}))
                    figObj(i).Visible='off';
                    figObj(i).Separator='off';
                end
                if strcmp(figObj(i).Tag,'FigureToolBar')
                    obj.h.toolbar=figObj(i);
                end
            end
            
            %prepare icons for moving between images
            load forwardBackwardsIcons;
            obj.h.prevFrameToolbar = uipushtool('Parent',obj.h.toolbar,'cdata',cdataBack, 'tooltip','previous frame', 'ClickedCallback',@obj.multiImagePrevFrameCallback);
            obj.h.nextFrameToolbar = uipushtool('Parent',obj.h.toolbar,'cdata',cdataForward, 'tooltip','next frame', 'ClickedCallback',@obj.multiImageNextFrameCallback);
            
            % define zoom options
            %obj.h.mainFigureZoom = zoom(obj.h.mainFigure); %get zoom handle
            %set(obj.h.mainFigureZoom,'Enable','on','Motion','Both','RightClickAction','PostContextMenu');
            
            % set file menus
            obj.h.fileMenu=uimenu(obj.h.mainFigure, 'Label', 'File' );
            obj.h.Transform2MultiTiffMenu=uimenu(obj.h.fileMenu,'Label','Transform stack 2 multi Tiff','Callback', {@obj.transformStack2MultiTiffCallback});
            obj.h.maxProjectionMenu=uimenu(obj.h.fileMenu,'Label','Calculate max projection','Callback', {@obj.maxProjectionImgCallback});
            
            % Arrange the main interface windows
            obj.h.mainHBox = uix.HBox('Parent',obj.h.mainFigure, 'Spacing',8);
            obj.h.controlResultsVBox = uix.VBox('Parent',obj.h.mainHBox, 'Spacing',4, 'Padding',4);
            obj.h.imagesVBox = uix.VBox('Parent',obj.h.mainHBox, 'Spacing',4, 'Padding',4);
            obj.h.mainHBox.Widths=[-1 -1];
            
            % arrange input image display
            obj.h.imagePanel(1) = uix.Panel('Parent', obj.h.imagesVBox, 'Padding', 2);
            obj.h.imageAxes(1) = axes('Parent',obj.h.imagePanel(1));
            obj.h.imagePanel(2) = uix.Panel('Parent', obj.h.imagesVBox, 'Padding', 2);
            obj.h.imageAxes(2) = axes('Parent',obj.h.imagePanel(2));
            obj.h.imageAxes(1).Position=obj.h.imageAxes(1).OuterPosition;
            obj.h.imageAxes(2).Position=obj.h.imageAxes(2).OuterPosition;
            axis([obj.h.imageAxes(1) obj.h.imageAxes(2)],'image');
            obj.h.imagesVBox.Heights=[-1 -1];
            
            % arrange control - results display
            obj.h.controlVBox=uix.VBox('Parent', obj.h.controlResultsVBox, 'Padding', 2);
            obj.h.regAxesPanel = uix.Panel('Parent', obj.h.controlResultsVBox, 'Padding', 2);
            obj.h.regAxes = axes('Parent',obj.h.regAxesPanel);
            obj.h.regAxes.Position=obj.h.regAxes.OuterPosition;
            axis(obj.h.regAxes,'image');
            obj.h.controlResultsVBox.Heights=[-1 -1];
            
            obj.h.loadImagePanel=uix.Panel('Parent', obj.h.controlVBox, 'Padding', 2);
            obj.h.transformPanel=uix.Panel('Parent', obj.h.controlVBox, 'Padding', 2);
            obj.h.controlVBox.Heights=[-1 -2];
            
            % load images
            obj.h.markPointsVBox = uix.VBox('Parent',obj.h.loadImagePanel);
            
            obj.h.pointImportExportHBox = uix.HBox('Parent',obj.h.markPointsVBox);
            obj.h.openPointSelectionGUI=uicontrol('Parent', obj.h.pointImportExportHBox, 'Style','push', 'String','open CPS','Callback',{@obj.OpenPointSelectionGUICallback});
            obj.h.getPointsFromPointSelectionGUI=uicontrol('Parent', obj.h.pointImportExportHBox, 'Style','push', 'String','import points','Callback',{@obj.getPointsFromPointSelectionGUICallback});
            obj.h.savePoint2File=uicontrol('Parent', obj.h.pointImportExportHBox, 'Style','push', 'String','save points','Callback',{@obj.saveMarkers2FileCallback});
            obj.h.loadPointfromFile=uicontrol('Parent', obj.h.pointImportExportHBox, 'Style','push', 'String','load points','Callback',{@obj.loadMarkersFromFileCallback});

            obj.h.loadImagesGridBox = uix.Grid('Parent',obj.h.markPointsVBox);
            obj.h.loadImage(1)=uicontrol('Parent', obj.h.loadImagesGridBox, 'Style','push', 'String','load image 1','Callback',{@obj.LoadImagePushCallback,1});
            obj.h.loadImage(2)=uicontrol('Parent', obj.h.loadImagesGridBox, 'Style','push', 'String','load image 2','Callback',{@obj.LoadImagePushCallback,2});
            
            obj.h.loadLayout(1)=uicontrol('Parent', obj.h.loadImagesGridBox, 'Style','push', 'String','load layout','Callback',{@obj.LoadLayoutPushCallback,1});
            obj.h.loadLayout(2)=uicontrol('Parent', obj.h.loadImagesGridBox, 'Style','push', 'String','load layout','Callback',{@obj.LoadLayoutPushCallback,2});
            
            obj.h.markPos(1)=uicontrol('Parent', obj.h.loadImagesGridBox, 'Style','push', 'String','mark','Callback',{@obj.MarkPosPushCallback,1});
            obj.h.markPos(2)=uicontrol('Parent', obj.h.loadImagesGridBox, 'Style','push', 'String','mark','Callback',{@obj.MarkPosPushCallback,2});
            
            obj.h.clearPos(1)=uicontrol('Parent', obj.h.loadImagesGridBox, 'Style','push', 'String','clear','Callback',{@obj.ClearPosPushCallback,1});
            obj.h.clearPos(2)=uicontrol('Parent', obj.h.loadImagesGridBox, 'Style','push', 'String','clear','Callback',{@obj.ClearPosPushCallback,2});
            
            set(obj.h.loadImagesGridBox,'Heights',[-1 -1]);

            % transform
            obj.h.transformHBox = uix.HBox('Parent',obj.h.transformPanel);
            obj.h.transformControlVBox = uix.VBox('Parent',obj.h.transformHBox);
            obj.h.transformList=uicontrol('Parent', obj.h.transformHBox, 'Style','listbox','Callback',@obj.TransformListCallback);
            obj.h.transformHBox.Widths=[-1 -1];
            
            obj.h.calculateTransform=uicontrol('Parent', obj.h.transformControlVBox, 'Style','push', 'String','Calculate','Callback',@obj.CalculateTransformPushCallback);
            obj.h.calculateListTransform=uicontrol('Parent', obj.h.transformControlVBox, 'Style','push', 'String','Calculate all list trans.','Callback',@obj.CalculateListTransformPushCallback,'TooltipString','Converts lower image (2) to upper image (1) going through transform list from top to bottom');
            obj.h.calculateSelectedListTransform=uicontrol('Parent', obj.h.transformControlVBox, 'Style','push', 'String','Calculate selected list trans.','Callback',@obj.CalculateSelectedListTranformCallback,'TooltipString','Converts lower image (2) to upper image (1)');
            obj.h.methodTransform=uicontrol('Parent', obj.h.transformControlVBox, 'Style','popupmenu','Callback',@obj.MethodTransformMenuCallback,'String',obj.par.transformStrings);
            obj.h.methodOverlay=uicontrol('Parent', obj.h.transformControlVBox, 'Style','popupmenu','String',obj.par.overlayMethodString);
            obj.h.saveTransform=uicontrol('Parent', obj.h.transformControlVBox, 'Style','push', 'String','Save to file','Callback',@obj.SaveTransformPushCallback);
            obj.h.loadTransform=uicontrol('Parent', obj.h.transformControlVBox, 'Style','push', 'String','Load to list','Callback',@obj.LoadTransformPushCallback);
            obj.h.addTransform=uicontrol('Parent', obj.h.transformControlVBox, 'Style','push', 'String','Add to list','Callback',@obj.AddTransformPushCallback);
            obj.h.deleteTransform=uicontrol('Parent', obj.h.transformControlVBox, 'Style','push', 'String','Delete from list','Callback',@obj.DeleteTransformPushCallback);
            obj.h.clearListTransform=uicontrol('Parent', obj.h.transformControlVBox, 'Style','push', 'String','Clear list','Callback',@obj.ClearTransformListPushCallback);
            obj.h.flipTransform=uicontrol('Parent', obj.h.transformControlVBox, 'Style','push', 'String','Flip','Callback',@obj.FlipTransformPushCallback);
            obj.h.moveTrasUp=uicontrol('Parent', obj.h.transformControlVBox, 'Style','push', 'String','move Up','Callback',{@obj.moveTrasUpDnCallback,1});
            obj.h.moveTrasDn=uicontrol('Parent', obj.h.transformControlVBox, 'Style','push', 'String','move Down','Callback',{@obj.moveTrasUpDnCallback,-1});
            
            obj.h.saveTransformedImages=uicontrol('Parent', obj.h.transformControlVBox, 'Style','push', 'String','save trans. image','Callback',@obj.SaveTransformedImagesCallback);
            obj.h.saveOverlayImage=uicontrol('Parent', obj.h.transformControlVBox, 'Style','push', 'String','save overlay image','Callback',@obj.SaveOverlayImageCallback);
            
            obj.h.mainFigure.Visible='on';
        end
    end
    
    methods (Access=protected)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%% callback functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function OpenPointSelectionGUICallback(obj,hObj,event)
            pos=getXYfromImPoint(obj);
            disp('Opening cpselect, please wait (it may take some time depending on the number of points you have');
            if ~isempty(obj.par.img{1}) & ~isempty(obj.par.img{2})
                if isempty(pos)
                    obj.h.cpselectHandle=cpselect(obj.par.img{1},obj.par.img{2});
                else
                    obj.h.cpselectHandle=cpselect(obj.par.img{1},obj.par.img{2},pos{1},pos{2});
                end
            else
                msgbox('Two images should be selected before opening point selection tool');
            end
        end
        
        function pos=getXYfromImPoint(obj)
            pos={};
            for img=[1 2]
                for i=1:obj.par.pointCounter(img)
                    pos{img}(i,:) = getPosition(obj.par.allignMarks{img}{i});
                end
            end
        end
        
        function saveMarkers2FileCallback(obj,hObj,event)
            [outFileName,outPathName] = uiputfile('*.mat');
            pos=getXYfromImPoint(obj);
            fileNames=obj.par.fullImageName;
            images=obj.par.img;
            save([outPathName outFileName],'pos','fileNames','images');
        end
        
        function loadMarkersFromFileCallback(obj,hObj,event)
            [FileName,PathName] = uigetfile('*.mat');
            load([PathName FileName],'pos');
            for img=[1 2]
                nPoints=size(pos{img},1);
                for i=1:nPoints
                    MarkPosPushCallback(obj,[],[],img,pos{img}(i,:));
                end
            end
        end
        
        function getPointsFromPointSelectionGUICallback(obj,hObj,event)
            disp('Exporting points...');
            pointObj=findobj(obj.h.cpselectHandle,'Type','hggroup','Tag','impoint');
            %point can be selected more specifically by checking the exact figure they belong to, first get the handle for the relevant parent and than check for this parent
            
            allPoints=1:numel(pointObj)/2; %all points for the lower images in cpselect
            img2Points=allPoints(1:numel(allPoints)/2); %the first markers are from image 2
            img1Points=allPoints(numel(allPoints)/2+1:numel(allPoints));
            
            for i=img1Points
                MarkPosPushCallback(obj,[],[],1,[pointObj(i).Children(5).XData pointObj(i).Children(5).YData])
                %obj.par.pointCoord{1}(obj.par.pointCounter(imageNum)+i,:)=[pointObj(i).Children(5).XData pointObj(i).Children(5).YData];
            end
            for i=img2Points
                MarkPosPushCallback(obj,[],[],2,[pointObj(i).Children(5).XData pointObj(i).Children(5).YData])
                %obj.par.pointCoord{2}(obj.par.pointCounter(imageNum)+i,:)=[pointObj(i).Children(5).XData pointObj(i).Children(5).YData];
            end
            
            delete(obj.h.cpselectHandle);
            disp('Done');
        end
        
        function maxProjectionImgCallback(obj,hObj,event)
            [FileName,PathName] = uigetfile({'*.*'},'Select files for max projection (for one file, takes all withs With BaseName_XXX)','MultiSelect','on');
            
            if isempty(FileName), return, end;
            if ischar(FileName)
                [~,name,ext] = fileparts(FileName);
                strDelimiters = strsplit(name,'_');
                choice = questdlg(['All files with ending of the form :' strDelimiters{end} 'will be selected']);
                if strcmp(choice,'Cancel'), return, end;
                FileName=dir([PathName strjoin(strDelimiters(1:end-1),'_') '_*' ext]);
                FileName={FileName.name};
            end
            [outFileName,outPathName] = uiputfile('*.tiff');
            maxProjectionFileName=[outPathName outFileName];
            nImages=numel(FileName);
            %Iinfo=imfinfo([PathName FileName{1}]);
            
            hWaitBar=waitbar(0,'Calculate max projection...');
            IMax=imread([PathName FileName{1}]);
            for i=2:nImages
                waitbar(i/nImages);
                Itmp=imread([PathName FileName{i}]);
                IMax=IMax+Itmp;
            end
            IMax = imadjust(IMax,stretchlim(IMax),[0 1]);
            imwrite(IMax,maxProjectionFileName);
            close(hWaitBar);
        end
        
        function transformStack2MultiTiffCallback(obj,hObj,event)
            [FileName,PathName] = uigetfile({'*.*'},'Select files to convert (for one file, takes all withs With BaseName_XXX)','MultiSelect','on');
            if isempty(FileName), return, end;
            if ischar(FileName)
                [~,name,ext] = fileparts(FileName);
                strDelimiters = strsplit(name,'_');
                choice = questdlg(['All files with ending of the form :' strDelimiters{end} 'will be selected']);
                if strcmp(choice,'Cancel'), return, end;
                FileName=dir([PathName strjoin(strDelimiters(1:end-1),'_') '_*' ext]);
                FileName={FileName.name};
            end
            [outFileName,outPathName] = uiputfile('*.tiff');
            multiTiffFileName=[outPathName outFileName];
            nImages=numel(FileName);
            %Iinfo=imfinfo([PathName FileName{1}]);
            
            options.append = true;
            options.big = false; % Use BigTIFF format
            options.comp = 'lzw';
            
            hWaitBar=waitbar(0,'Converting to multiTiff file...');
            for i=1:nImages
                waitbar(i/nImages);
                Itmp=imread([PathName FileName{i}]);
                %Itmp = imadjust(Itmp,stretchlim(Itmp),[0 1]);
                %Itmp=loadtiff([PathName FileName{i}]);
                saveastiff(Itmp, multiTiffFileName, options);
                %imwrite(Itmp,multiTiffFileName,'WriteMode','append');
            end
            close(hWaitBar);
        end
        
        function clearAxis(obj,imageNum)
            delete(obj.h.img{imageNum});
            for i=1:numel(obj.par.allignMarks{imageNum})
                delete(obj.par.allignMarks{imageNum}{i});
            end
            for i=1:numel(obj.h.layoutPlot{imageNum})
                delete(obj.h.layoutPlot{imageNum}{i});
            end
            obj.par.pointCounter(imageNum)=0;
            obj.par.multiImageIdx=1;
            childrenList=get(obj.h.imageAxes(imageNum),'Children');
            delete(childrenList);
            set(obj.h.imageAxes(imageNum),'Visible','on','TickDir','in','DataAspectRatioMode','auto','Box','off','YDir','normal');
        end
        
        function LoadImagePushCallback(obj,hObj,event,imageNum)
            obj.clearAxis(imageNum);
            [obj.par.imgFileName{imageNum},obj.par.imgPathName{imageNum}] = uigetfile('*.*','loadImage',cd);
            if obj.par.imgFileName{imageNum}==0
                disp('No figure selected');
                return;
            end
            obj.par.fullImageName{imageNum}=[obj.par.imgPathName{imageNum} obj.par.imgFileName{imageNum}];
            obj.par.imgInfo{imageNum}=imfinfo(obj.par.fullImageName{imageNum});
            if numel(obj.par.imgInfo{imageNum})>1
                obj.par.multiImageFile(imageNum)=true;
                obj.par.multiImageIdx=1;
                obj.par.multiImages=numel(obj.par.imgInfo{imageNum});
            else
                obj.par.multiImageFile(imageNum)=false;
            end
            obj.par.img{imageNum}=imread(obj.par.fullImageName{imageNum});
            if size(obj.par.img{imageNum},3)==4
                obj.par.img{imageNum}(:,:,4)=[];
                disp('Image had a transparency channel. The information in this channel was discarded');
            end
            obj.h.img{imageNum}=imshow(obj.par.img{imageNum},'Parent',obj.h.imageAxes(imageNum));
            obj.par.isElectrodeLayout{imageNum}=false;
        end
        
        function multiImageNextFrameCallback(obj,hObj,event)
            pMultiFile=find(obj.par.multiImageFile);
            if isempty(pMultiFile)
                disp('Frame change is only supported for multi-image files')
            else
                if (obj.par.multiImageIdx+1)<=obj.par.multiImages
                    delete(obj.h.img{pMultiFile});
                    obj.par.img{pMultiFile}=imread(obj.par.fullImageName{pMultiFile},obj.par.multiImageIdx+1);
                    obj.h.img{pMultiFile}=imshow(obj.par.img{pMultiFile},'Parent',obj.h.imageAxes(pMultiFile));
                    obj.par.multiImageIdx=obj.par.multiImageIdx+1;
                else
                    disp('Last image in multi-image file');
                end
            end
        end
        
        function multiImagePrevFrameCallback(obj,hObj,event)
            pMultiFile=find(obj.par.multiImageFile);
            if isempty(pMultiFile)
                disp('Frame change is only supported for multi-image files')
            else
                if (obj.par.multiImageIdx-1)>=1
                    delete(obj.h.img{pMultiFile});
                    obj.par.img{pMultiFile}=imread(obj.par.fullImageName{pMultiFile},obj.par.multiImageIdx-1);
                    obj.h.img{pMultiFile}=imshow(obj.par.img{pMultiFile},'Parent',obj.h.imageAxes(pMultiFile));
                    obj.par.multiImageIdx=obj.par.multiImageIdx-1;
                else
                    disp('Last image in multi-image file');
                end
            end
        end
        
        function MarkPosPushCallback(obj,hObj,event,imageNum,point)
            %[x, y] = getpts(ax)
            obj.par.pointCounter(imageNum)=obj.par.pointCounter(imageNum)+1;
            posConstFcn = makeConstrainToRectFcn('impoint',obj.h.imageAxes(imageNum).XLim,obj.h.imageAxes(imageNum).YLim); % this line should be moved to image initialization
            if nargin~=5
                obj.par.allignMarks{imageNum}{obj.par.pointCounter(imageNum)} = impoint(obj.h.imageAxes(imageNum),'PositionConstraintFcn',posConstFcn);
            else
                obj.par.allignMarks{imageNum}{obj.par.pointCounter(imageNum)} = impoint(obj.h.imageAxes(imageNum),point(1),point(2),'PositionConstraintFcn',posConstFcn);
            end
            
            obj.par.pointCoord{imageNum}(obj.par.pointCounter(imageNum),:)=getPosition(obj.par.allignMarks{imageNum}{obj.par.pointCounter(imageNum)});
            
            if obj.par.isElectrodeLayout{imageNum}
                [~,p]=min((obj.par.pointCoord{imageNum}(obj.par.pointCounter(imageNum),1)-obj.par.Xc{imageNum}).^2+(obj.par.pointCoord{imageNum}(obj.par.pointCounter(imageNum),2)-obj.par.Yc{imageNum}).^2);
                setPosition(obj.par.allignMarks{imageNum}{obj.par.pointCounter(imageNum)},[obj.par.Xc{imageNum}(p) obj.par.Yc{imageNum}(p)]);
                obj.par.pointCoord{imageNum}(obj.par.pointCounter(imageNum),:)=getPosition(obj.par.allignMarks{imageNum}{obj.par.pointCounter(imageNum)});
            end
            
            setString(obj.par.allignMarks{imageNum}{obj.par.pointCounter(imageNum)},num2str(obj.par.pointCounter(imageNum)));
        end
        
        function ClearPosPushCallback(obj,hObj,event,imageNum)
            for i=1:numel(obj.par.allignMarks{imageNum})
                delete(obj.par.allignMarks{imageNum}{i});
            end
            obj.par.pointCoord{imageNum}=[];
            obj.par.pointCounter(imageNum)=0;
        end
        
        function LoadLayoutPushCallback(obj,hObj,event,imageNum)
            obj.clearAxis(imageNum);
            obj.par.isElectrodeLayout{imageNum}=true;
            gridColor=[0.8 0.8 0.8];
            elecNumFontSize=12;
            layoutDir=fileparts(which('layout_100_12x12.mat'));
            [obj.par.layoutName] = uigetfile('layout_*.mat','loadImage',layoutDir);
            if obj.par.layoutName==0
                disp('No layout selected');
                return;
            end
            obj.par.layout=load(obj.par.layoutName);
            
            %obj.par.layout=flipud(obj.par.layout.En);
            obj.par.layout=obj.par.layout.En;
            obj.par.electrodePitch=textscan(obj.par.layoutName(1:end-4), '%s', 'delimiter', '_');
            obj.par.electrodePitch=str2num(obj.par.electrodePitch{1}{2});
            obj.par.channels=unique(obj.par.layout(~isnan(obj.par.layout)));
            
            [meshY,meshX]=meshgrid(1:size(obj.par.layout,1),1:size(obj.par.layout,2));
            meshX=meshX';meshY=meshY';
            obj.par.Xc{imageNum}(obj.par.layout(~isnan(obj.par.layout)))=meshX(~isnan(obj.par.layout))*obj.par.electrodePitch;
            obj.par.Yc{imageNum}(obj.par.layout(~isnan(obj.par.layout)))=meshY(~isnan(obj.par.layout))*obj.par.electrodePitch;
            
            %create electrode grid with numbers for plotting
            [max_grid_y,max_grid_x]=size(obj.par.layout);
            XL=[0.5 max_grid_x+0.5]*obj.par.electrodePitch;
            YL=[0.5 max_grid_y+0.5]*obj.par.electrodePitch;
            xM=(1.5:(max_grid_y-0.5))'*obj.par.electrodePitch;
            yM=(1.5:(max_grid_x-0.5))'*obj.par.electrodePitch;
            
            
            obj.h.layoutPlot{imageNum}{1}=line([XL(1) XL(2)],[xM xM],'LineWidth',1,'Color',gridColor,'Parent',obj.h.imageAxes(imageNum));
            obj.h.layoutPlot{imageNum}{2}=line([yM yM],[YL(1) YL(2)],'LineWidth',1,'Color',gridColor,'Parent',obj.h.imageAxes(imageNum));
            obj.h.layoutPlot{imageNum}{3}=text(obj.par.Xc{imageNum},obj.par.Yc{imageNum},num2str(obj.par.channels),'fontsize',elecNumFontSize,'horizontalAlignment','center','Parent',obj.h.imageAxes(imageNum));
            xlim(obj.h.imageAxes(imageNum),[XL(1) XL(2)]);
            ylim(obj.h.imageAxes(imageNum),[YL(1) YL(2)]);
        end
        
        function CalculateTransformPushCallback(obj,hObj,event)
            %check if the number of points is equal in both plots
            if ~isfield(obj.par,'pointCoord')
                disp('Cant calculate transformation since no anchor points were chosen. To calculate transform from the list use "Calculate from list"');
                return;
            end
            if any(size(obj.par.pointCoord{1})~=size(obj.par.pointCoord{2}))
                disp('Cant calculate transformation due to unequal number of transform points in the two plots');
                return;
            end
            
            % Create a transformation structure for transformation (from image 2 to image 1)
            obj.par.currentTrans = cp2tform(obj.par.pointCoord{2},obj.par.pointCoord{1},obj.par.currentTransformName); %for two point transformation
            obj.applyCurrentTransform;
        end
        
        function applyCurrentTransform(obj)
            
            disp('Calculating transform...');
            childrenList=get(obj.h.regAxes,'Children');
            delete(childrenList);
            overlaymethod=obj.par.overlayMethodNames{obj.h.methodOverlay.Value};
            
            if obj.par.isElectrodeLayout{2} && ~obj.par.isElectrodeLayout{1}
                [obj.par.transPointsX, obj.par.transPointsY] = tformfwd(obj.par.currentTrans, obj.par.Xc{2}, obj.par.Yc{2});
                imshow(obj.par.img{1},'Parent',obj.h.regAxes);hold(obj.h.regAxes,'on');
                cMap=jet(numel(obj.par.channels));
                scatter(obj.h.regAxes,obj.par.transPointsX,obj.par.transPointsY,10,cMap(obj.par.channels,:));hold(obj.h.regAxes,'off');
                [obj.par.mergedImg] = frame2im(getframe(obj.h.regAxes));
            elseif ~obj.par.isElectrodeLayout{2} && obj.par.isElectrodeLayout{1}
                obj.par.transImg = imtransform(obj.par.img{2},obj.par.currentTrans,'FillValues',0,'XData',[0 max(obj.par.Xc{1})+100], 'YData',[0 max(obj.par.Yc{1})+100],'XYScale',1);
                imshow(obj.par.transImg,'Parent',obj.h.regAxes);hold(obj.h.regAxes,'on');
                cMap=jet(numel(obj.par.channels));
                scatter(obj.h.regAxes,obj.par.Xc{1},obj.par.Yc{1},10,cMap(obj.par.channels,:));hold(obj.h.regAxes,'off');
                set(obj.h.regAxes,'YDir','normal');
                [obj.par.mergedImg] = frame2im(getframe(obj.h.regAxes));
            elseif ~obj.par.isElectrodeLayout{2} && ~obj.par.isElectrodeLayout{1}
                if strcmp(overlaymethod,'sourceDominance') || strcmp(overlaymethod,'transDominance')
                    maxSample(1)=2^(obj.par.imgInfo{1}(1).BitDepth/size(obj.par.img{1},3))-1;
                    maxSample(2)=2^(obj.par.imgInfo{2}(1).BitDepth/size(obj.par.img{2},3))-1;
                    obj.par.transImg = imtransform(obj.par.img{2},obj.par.currentTrans,'FillValues',maxSample(2),'XData',[1 obj.par.imgInfo{1}(1).Width], 'YData',[1 obj.par.imgInfo{1}(1).Height],'XYScale',1);
                    if strcmp(overlaymethod,'sourceDominance')
                        obj.par.mergedImg=obj.par.transImg;
                        obj.par.mergedImg(obj.par.img{1}~=maxSample(1))=obj.par.img{1}(obj.par.img{1}~=maxSample(1));
                    else
                        if size(obj.par.transImg,3)==3 && size(obj.par.img{1},3)==1
                            obj.par.mergedImg=cast(zeros([obj.par.imgInfo{1}(1).Height,obj.par.imgInfo{1}(1).Width,3]),'like',obj.par.img{1});
                            obj.par.mergedImg(:,:,3)=obj.par.img{1};
                            obj.par.mergedImg(:,:,2)=obj.par.img{1};
                            obj.par.mergedImg(:,:,1)=obj.par.img{1};
                        else
                            obj.par.mergedImg=obj.par.img{1};
                        end
                        p=obj.par.transImg(:,:,1)~=maxSample(2) & obj.par.transImg(:,:,2)~=maxSample(2) & obj.par.transImg(:,:,3)~=maxSample(2);
                        p=repmat(p,[1 1 3]);
                        obj.par.mergedImg(p)=cast(obj.par.transImg(p),'like',obj.par.img{1}).*((maxSample(1)+1)/(maxSample(2)+1));
                    end
                elseif strcmp(overlaymethod,'blend')
                    obj.par.transImg = imtransform(obj.par.img{2},obj.par.currentTrans,'FillValues',0,'XData',[1 obj.par.imgInfo{1}(1).Width], 'YData',[1 obj.par.imgInfo{1}(1).Height],'XYScale',1);
                    obj.par.mergedImg = imfuse(obj.par.img{1},obj.par.transImg,overlaymethod);
                elseif strcmp(overlaymethod,'falsecolor')
                    obj.par.transImg = imtransform(obj.par.img{2},obj.par.currentTrans,'FillValues',0,'XData',[1 obj.par.imgInfo{1}(1).Width], 'YData',[1 obj.par.imgInfo{1}(1).Height],'XYScale',1);
                    obj.par.mergedImg = imfuse(obj.par.img{1},obj.par.transImg,overlaymethod,'ColorChannels','green-magenta');
                end
                imshow(obj.par.mergedImg,'Parent',obj.h.regAxes);
            elseif obj.par.isElectrodeLayout{2} && obj.par.isElectrodeLayout{1}
                disp('Transforms from electrode layout to electrode layout are not supported');
                return;
            end
            disp('Done!');
        end
        
        function CalculateListTransformPushCallback(obj,hObj,event)
            transforms=obj.par.transformsInList;
            obj.par.currentTrans=maketform('composite',flipud(cell2mat(transforms)));
            obj.applyCurrentTransform;
            disp('Done!');
        end
        
        function MethodTransformMenuCallback(obj,hObj,event)
            obj.par.currentTransformName=obj.par.trasformNames{hObj.Value};
            obj.par.currentTransformMinPoints=obj.par.trasformMinPoints(hObj.Value);
        end
        
        function LoadTransformPushCallback(obj,hObj,event)
            [transformFilename,transformPath] = uigetfile('*.mat','Load transformes to list',cd,'MultiSelect','on');
            if isnumeric(transformFilename) %meaning no files were load (opperation cancelled
                disp('Transform not selected');
                return;
            elseif ischar(transformFilename)
                transformFilename={transformFilename};
            end
            for i=1:numel(transformFilename)
                tmp=load([transformPath transformFilename{i}]);
                nTransforms=numel(obj.h.transformList.String);
                obj.h.transformList.Value=max(nTransforms,1);
                obj.h.transformList.String=[obj.h.transformList.String;tmp.transformNames];
                obj.par.transformsInList=[obj.par.transformsInList;tmp.transforms];
                
            end
            disp('Transforms loaded!');
        end
        
        function SaveTransformPushCallback(obj,hObj,event)
            if isempty(obj.h.transformList.String)
                disp('No transforms calculated, cant save anything');
            else
                strPartsFirst=strsplit(obj.h.transformList.String{1},'->');
                strPartsLast=strsplit(obj.h.transformList.String{end},'->');
                [imgFileName,imgPathName] = uiputfile('*.mat','Save transformes in list',[cd filesep 'T_' strPartsFirst{1} '-' strPartsLast{2}]);
                if imgFileName==0
                    disp('Transform save cancelled');
                    return;
                end
                transformNames=obj.h.transformList.String;
                transforms=obj.par.transformsInList;
                save([imgPathName imgFileName],'transformNames','transforms');
                nTransforms=numel(obj.h.transformList.String);
                obj.h.transformList.Value=max(nTransforms,1);
                disp('Transform saved!');
            end
        end
        
        function AddTransformPushCallback(obj,hObj,event)
            defaultVal = {'Image 1','Image 2'};
            inputNames = inputdlg({'Source figure','warped figure'},'Annotate transformation',1,defaultVal);
            if numel(inputNames)==0
                disp('Adding a new transform to list cancelled!!!');
                return;
            end
            nTransforms=numel(obj.h.transformList.String);
            obj.h.transformList.String{nTransforms+1}=[inputNames{2} '->' inputNames{1}];
            obj.par.transformsInList{nTransforms+1,1}=obj.par.currentTrans;
            obj.h.transformList.Value=nTransforms+1;
        end
        
        function CalculateSelectedListTranformCallback(obj,hObj,event)
            selectedTransform=obj.h.transformList.Value;
            obj.par.currentTrans=obj.par.transformsInList{selectedTransform};
            obj.applyCurrentTransform;
        end
        
        function ClearTransformListPushCallback(obj,hObj,event)
            obj.h.transformList.String='';
            obj.par.transformsInList=[];
            obj.h.transformList.Value=1;
        end
        
        function DeleteTransformPushCallback(obj,hObj,event)
            selectedTransform=obj.h.transformList.Value;
            obj.h.transformList.Value=max(selectedTransform-1,1);
            obj.h.transformList.String(selectedTransform)=[];
            obj.par.transformsInList(selectedTransform)=[];
        end
        
        function TransformListCallback(obj,hObj,event)
        end
        
        function FlipTransformPushCallback(obj,hObj,event)
            selectedTransform=obj.h.transformList.Value;
            sourceWarpedImageNames=textscan(obj.h.transformList.String{selectedTransform}, '%s', 'delimiter', '->','MultipleDelimsAsOne',true);
            obj.par.transformsInList{selectedTransform}=fliptform(obj.par.transformsInList{selectedTransform}); %inverse the transform direction
            obj.h.transformList.String{selectedTransform}=[sourceWarpedImageNames{1}{2} '->' sourceWarpedImageNames{1}{1}];
        end
        
        function moveTrasUpDnCallback(obj,hObj,event,moveWhere)
            selectedTransform=obj.h.transformList.Value;
            
            if moveWhere==1 %move up
                if selectedTransform==1
                    disp('Cant move up, since I am already on the top');
                    return;
                end
            elseif moveWhere==-1 %move up
                nTrasforms=numel(obj.h.transformList.String);
                if selectedTransform==nTrasforms
                    disp('Cant move down, since I am already on the botom');
                    return;
                end
            end
            tmpStr=obj.h.transformList.String{selectedTransform};
            tmpTrans=obj.par.transformsInList{selectedTransform};
            obj.h.transformList.String{selectedTransform}=obj.h.transformList.String{selectedTransform-moveWhere};
            obj.par.transformsInList{selectedTransform}=obj.par.transformsInList{selectedTransform-moveWhere};
            obj.h.transformList.String{selectedTransform-moveWhere}=tmpStr;
            obj.par.transformsInList{selectedTransform-moveWhere}=tmpTrans;
            obj.h.transformList.Value=obj.h.transformList.Value-moveWhere;
        end
        
        function SaveTransformedImagesCallback(obj,hObj,event)
            [imgFileName,imgPathName] = uiputfile('*.*','Save transformed image',cd);
            if imgFileName==0
                disp('Save image cancelled');
                return;
            end
            transformedImage=obj.par.transImg;
            save([imgPathName imgFileName],'transformedImage');
            imwrite(obj.par.transImg,[imgPathName imgFileName '.tiff'],'tiff');
            disp('Transformed image saved!');
        end
        
        function SaveOverlayImageCallback(obj,hObj,event)
            export2separateFig=false;
            saveOverlayMatfile=false;
            if obj.par.isElectrodeLayout{1} || obj.par.isElectrodeLayout{2}
                export2separateFig=true;
            end
            
            pMultiImage=find(obj.par.multiImageFile);
            pref_value='no';
            if ~isempty(pMultiImage)
                pref_value = uigetpref('tmpGroup','tmpPref','save options','Save overlay using last calculated transform for all images in the multi-image file?',{'yes','no'});
            end
            if strcmp(pref_value,'no')
                [imgFileName,imgPathName] = uiputfile('*.*','Save overlay image',[cd filesep 'overlay_']);
                if imgFileName==0
                    disp('Save image cancelled');
                    return;
                end
                overlayImage=obj.par.mergedImg;
                imwrite(obj.par.mergedImg,[imgPathName imgFileName],'tiff');
                %export_fig(obj.h.regAxes,[imgPathName imgFileName],'-tif','-r500');
                if saveOverlayMatfile
                    save([imgPathName imgFileName],'overlayImage');
                end
                
                if export2separateFig
                    %cmap=get(obj.h.regAxes,'colormap');
                    hNewFig=figure('Position',[50 100 500 500]);
                    copyobj(obj.h.regAxes,hNewFig);
                    set(hNewFig,'PaperPositionMode','auto');
                    assignin('base','exportedPlot',hNewFig);
                end
                disp('Overlayed image saved!');
            else
                [imgFileName,imgPathName] = uiputfile('*.*','Save overlay image',[cd filesep 'overlay_']);
                if imgFileName==0
                    disp('Save image cancelled');
                    return;
                end
                obj.par.multiImageIdx=0;
                for i=1:obj.par.multiImages
                    disp(['Calculating transform for image ' num2str(i) ' of stack']); 
                    multiImageNextFrameCallback(obj)
                    applyCurrentTransform(obj);
                    imwrite(obj.par.mergedImg,[imgPathName imgFileName '.tif'],'WriteMode','append');
                end
                disp('Done!');
            end
        end
        
        function closeMainGUIFigure(obj,hObj,event)
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%% end callback functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
end