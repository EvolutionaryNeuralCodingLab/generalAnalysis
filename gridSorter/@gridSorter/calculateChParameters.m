function obj=calculateChParameters(obj)
%Populate the chPar structure with morphological channel parameters to be later used in soring
%test that a valid channel layout is available
if isempty(obj.dataRecordingObj.chLayoutNumbers)
    disp('No channel layout found in recording object,exiting...');
    return;
end

if obj.localGridSize/2==round(obj.localGridSize/2)
    error('Error constructing channel parameter: localGridSize should be an even integer');
end
%determine channels
if isempty(obj.selectedChannelSubset)
    obj.selectedChannelSubset=obj.dataRecordingObj.channelNumbers; %take all channel in the recording
    disp('No specific channels selected: Sorting all channels in recording');
end
obj.nCh=numel(obj.selectedChannelSubset);
obj.chPar.s2r=obj.selectedChannelSubset; %transformation between the serial channel number and the real ch number
obj.chPar.r2s(obj.selectedChannelSubset)=1:obj.nCh; %transformation between the real channel number and the serial ch number

%get channel layout
obj.chPar.En=obj.dataRecordingObj.chLayoutNumbers;

%limit layout only to the serial channels selected by user
obj.chPar.rEn=obj.chPar.En;
obj.chPar.rEn(~isnan(obj.chPar.En))=obj.chPar.r2s(obj.chPar.En(~isnan(obj.chPar.En)));
obj.chPar.rEn(obj.chPar.rEn==0)=NaN;
[nRowsTmp,nColsTmp]=size(obj.chPar.En);

localGridSizeExt=obj.localGridSize+obj.localGridExt*2;

%initiate arrays
obj.chPar.surChExt=cell(1,obj.nCh);obj.chPar.pValidSurChExt=cell(1,obj.nCh);obj.chPar.surChExtVec=cell(1,obj.nCh);obj.chPar.pCenterCh=zeros(1,obj.nCh);obj.chPar.nValidChExt=zeros(1,obj.nCh);obj.chPar.pSurCh=cell(1,obj.nCh);obj.chPar.pSurChOverlap=cell(1,obj.nCh);
if ~isempty(localGridSizeExt)
    obj.arrayExt=(obj.localGridSize-1)/2;
    overheadGridExtension=(localGridSizeExt-1)/2;
    EnExt=NaN(nRowsTmp+overheadGridExtension*2,nColsTmp+overheadGridExtension*2);
    EnExt(1+overheadGridExtension:end-overheadGridExtension,1+overheadGridExtension:end-overheadGridExtension)=obj.chPar.rEn;
    for i=1:obj.nCh
        [x,y]=find(EnExt==i);
        
        %find the surrounding channels on which feature extraction will be performed
        surCh=EnExt(x-obj.arrayExt:x+obj.arrayExt,y-obj.arrayExt:y+obj.arrayExt);
        pValidSurCh=find(~isnan(surCh)); %do not remove find
        
        %find the channels that are overhead step from the central channel - these are the channels who's waveforms should be checked for merging
        surChOverlap=EnExt(x-obj.localGridExt:x+obj.localGridExt,y-obj.localGridExt:y+obj.localGridExt);
        surChOverlap(obj.localGridExt+1,obj.localGridExt+1)=NaN;
        pValidSurChOverlap=find(~isnan(surChOverlap)); %do not remove find
        
        %find the extended channels for merging of the same neurons detected on nearby channels
        obj.chPar.surChExt{i}=EnExt(x-overheadGridExtension:x+overheadGridExtension,y-overheadGridExtension:y+overheadGridExtension);
        obj.chPar.pValidSurChExt{i}=find(~isnan(obj.chPar.surChExt{i})); %do not remove find
        obj.chPar.nValidChExt(i)=numel(obj.chPar.pValidSurChExt{i});
        
        obj.chPar.surChExtVec{i}=obj.chPar.surChExt{i}(obj.chPar.pValidSurChExt{i}(:))';
        obj.chPar.pCenterCh(i)=find(obj.chPar.surChExtVec{i}==i); %the position of the central channel in surChExtVec
        [~,obj.chPar.pSurCh{i}]=intersect(obj.chPar.surChExtVec{i},surCh(pValidSurCh(:)));
        [~,obj.chPar.pSurChOverlap{i}]=intersect(obj.chPar.surChExtVec{i},surChOverlap(pValidSurChOverlap(:)));
    end
end
%map channel intersections
for i=1:obj.nCh
    for j=obj.chPar.surChExtVec{i}(obj.chPar.pSurChOverlap{i}) %the trivial case : (i,i) is also included (can be remove in not required)
        [obj.chPar.sharedChNames{i}{j},obj.chPar.pSharedCh1{i}{j},obj.chPar.pSharedCh2{i}{j}]=intersect(obj.chPar.surChExtVec{i},obj.chPar.surChExtVec{j});
    end
end

