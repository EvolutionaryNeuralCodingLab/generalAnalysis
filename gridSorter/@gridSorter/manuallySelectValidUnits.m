function [obj]=manuallySelectValidUnits(obj)
%manually select valid spike units after sorting

obj=obj.findSortingFiles;

if ~obj.sortingFileNames.fittingExist || ~all(obj.sortingFileNames.postProcessingAnalysisExist)
    error('Cant run manual selection since spike sorting on this data set was not completed');
end

load(obj.sortingFileNames.fittingFile,'t','ic');
plotNames=dir([obj.sortingDir filesep '*spikeShape.jpg']);
plotNames=cellfun(@(x) [obj.sortingDir filesep x],{plotNames.name},'UniformOutput',0)';

f=figure('Position',[100 50 1100 900]);
uipAxes = uipanel('Parent',f,'Title','Spike plot','FontSize',12,'Position',[.005 .1 0.99 .9]);
uipControl = uipanel('Parent',f,'Title','Spike selection','FontSize',12,'Position',[.005 .005 0.99 .09]);
pushAccept = uicontrol('Parent',uipControl,'String','Accept','Position',[10 10 72 36],'Callback',{@acceptReject_Callback,true});
pushReject = uicontrol('Parent',uipControl,'String','Reject','Position',[100 10 72 36],'Callback',{@acceptReject_Callback,false});

hAxesSpike=axes('Parent',uipAxes,'Position',[.005 .005 0.6 0.9]);
hAxesHistHighRes=axes('Parent',uipAxes,'Position',[.62 .005 .35 .4]);
hAxesHist=axes('Parent',uipAxes,'Position',[.62 .5 .35 .4]);

%set(hAxesSpike,'Position',get(hAxesSpike,'OuterPosition'));

edgesH=[0:10:500];
edgesL=[0:100:10000];

nNeu=size(ic,2);
for i=1:size(ic,2)
    printFile=[obj.sortingDir filesep 'neuron' 'Neu ' num2str(ic(1,i)) '-' num2str(ic(2,i)) '-spikeShape.jpg'];
    if exist(printFile,'file')
        
        
        histHigh=histc(diff(t(ic(3,i):ic(4,i))),edgesH);
        histLow=histc(diff(t(ic(3,i):ic(4,i))),edgesL);
        bar(hAxesHistHighRes,(edgesH(1:end-1)+edgesH(2:end))/2,histHigh(1:end-1));
        bar(hAxesHist,(edgesL(1:end-1)+edgesL(2:end))/2,histLow(1:end-1));
        
        imshow(printFile,'Parent',hAxesSpike);
        set(pushAccept,'Selected','off');
        set(pushReject,'Selected','off');
        uiwait;
        selectedTrue=get(pushAccept,'Selected');
        RejectedTrue=get(pushReject,'Selected');
        if strcmp(selectedTrue,'on') & strcmp(RejectedTrue,'off')
            neurons2keep(i)=true;
        else
            neurons2keep(i)=false;
        end
    else
        disp(['No plot exists for unit: ' num2str(ic(1:2,i)')]);
    end
end

[t,ic]=RemainNeurons(t,ic,ic(1:2,neurons2keep));
save([obj.sortingDir filesep 'manuallySelectedSpikeSorting'],'t','ic','neurons2keep');

delete(f);

function acceptReject_Callback(hObject, eventdata, handle, isValid)
set(hObject,'Selected','on');
uiresume(gcbf);