function obj=getSpikeTrigWF(obj)
%Detect spikes in raw data and save waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error('This method is an old version, please use MEAAnalysis.getSpikeTrigWF instead');

load(obj.sortingFileNames.fittingFile,'t','ic');
saveFileName='tmpPostProcessing';
neuronNames=ic(1:2,:);

%determine quantization
obj.detectionInt2uV=obj.detectionMaxSpikeAmp/2^(obj.detectionNQuantizationBits-1);

obj=obj.getHighpassFilter;

%start spike detection
fprintf('\nRunning waveform extraction on %s...',obj.dataRecordingObj.dataFileNames{1});

%determine the chunck size
if obj.detectionMaxChunkSize>obj.dataRecordingObj.recordingDuration_ms
    startTimes=0;
    endTimes=obj.dataRecordingObj.recordingDuration_ms;
else
    startTimes=(0:obj.detectionMaxChunkSize:obj.dataRecordingObj.recordingDuration_ms);
    endTimes=[startTimes(2:end) obj.dataRecordingObj.recordingDuration_ms];
end
nChunks=numel(startTimes);

obj.nCh=numel(obj.chPar.s2r);
nNeurons=size(ic,2);
nSpkTotal=ic(4,:)-ic(3,:)+1;

windowSamplesRaw=obj.postTotalRawWindow*obj.dataRecordingObj.samplingFrequency/1000;
windowSamplesFiltered=obj.postTotalFilteredWindow*obj.dataRecordingObj.samplingFrequency/1000;
HPstartIdx=((obj.postTotalFilteredWindow-obj.postPreFilteredWindow)*obj.dataRecordingObj.samplingFrequency/1000+1);
HPIdx=HPstartIdx:(HPstartIdx+windowSamplesFiltered-1);

keepAllInMemory=1;
if ispc
    userview = memory;
    doubleBytes=8;
    if userview.MemAvailableAllArrays<(nNeurons*obj.nCh*windowSamplesRaw*doubleBytes*4) & keepAllInMemory
        disp('Arrays too big for memory, moving to saving on disk (for using memory change chunk size');
        keepAllInMemory=0;
    end
end

par.binTRaw=1; %ms
par.binTHP=0.1; %ms
par.binV=0.5;
par.maxV=75;

nBinsRaw=ceil(obj.postTotalRawWindow/par.binTRaw);
nBinsHP=ceil(obj.postTotalFilteredWindow/par.binTHP);
nBinsV=2*par.maxV/par.binV;

if keepAllInMemory
    avgRawWF=zeros(nNeurons,obj.nCh,windowSamplesRaw);
    stdRawWF=zeros(nNeurons,obj.nCh,windowSamplesRaw);
    avgHPWF=zeros(nNeurons,obj.nCh,windowSamplesFiltered);
    stdHPWF=zeros(nNeurons,obj.nCh,windowSamplesFiltered);
    histRawWF=zeros(nNeurons,obj.nCh,nBinsRaw,nBinsV,'uint8');
    histHPWF=zeros(nNeurons,obj.nCh,nBinsHP,nBinsV,'uint8');
else
    matFileObj = matfile(obj.sortingFileNames.STWaveformFile,'Writable',true);
    if obj.postExtractRawLongWaveformsFromSpikeTimes
        matFileObj.avgRawWF=zeros(nNeurons,obj.nCh,windowSamplesRaw);
        matFileObj.stdRawWF=zeros(nNeurons,obj.nCh,windowSamplesRaw);
        matFileObj.histRawWF=zeros(nNeurons,obj.nCh,nBinsRaw,nBinsV,'uint8');
    else
        matFileObj.avgRawWF=[];matFileObj.stdRawWF=[];matFileObj.histRawWF=[];
    end
    
    if obj.postExtractFilteredWaveformsFromSpikeTimes
        matFileObj.avgHPWF=zeros(nNeurons,obj.nCh,windowSamplesFiltered);
        matFileObj.stdHPWF=zeros(nNeurons,obj.nCh,windowSamplesFiltered);
        matFileObj.histHPWF=zeros(nNeurons,obj.nCh,nBinsHP,nBinsV,'uint8');
    else
        matFileObj.avgHPWF=[];matFileObj.stdHPWF=[];matFileObj.histHPWF=[];
    end
end

%load files if matlab crashed during analysis
if exist([obj.sortingDir filesep saveFileName '_tmp.mat'],'file')
    load([obj.sortingDir filesep saveFileName '_tmp.mat']);
    startChunk=max(1,lastGoodChunck+1);
else
    startChunk=1;
end

tBinRaw=shiftdim(ceil(((1:windowSamplesRaw)/obj.dataRecordingObj.samplingFrequency*1000)/par.binTRaw),-1);
tBinHP=shiftdim(ceil(((1:windowSamplesFiltered)/obj.dataRecordingObj.samplingFrequency*1000)/par.binTHP),-1);
VBins=(-par.maxV+par.binV/2):par.binV:par.maxV;
pPreBaselineSamples=((obj.dataRecordingObj.samplingFrequency(1)/1000)*(obj.postPreRawWindow-obj.detectionPreSpikeWindow-1)):((obj.dataRecordingObj.samplingFrequency(1)/1000)*(obj.postPreRawWindow-obj.detectionPreSpikeWindow));
tBinsRaw=par.binTRaw/2:par.binTRaw:obj.postTotalRawWindow;
tBinsHP=par.binTHP/2:par.binTHP:obj.postTotalFilteredWindow;
try
    tic;
    %initiate arrays
    fprintf('\nExtracting spikes from chunks (total %d): ',nChunks);
    for i=startChunk:nChunks
        fprintf('%d',i);
        %get data
        MAll=obj.dataRecordingObj.getData(obj.chPar.s2r(1:obj.nCh),startTimes(i)-obj.postPreRawWindow,endTimes(i)-startTimes(i)+obj.postTotalRawWindow);
        MFAll=squeeze(obj.filterObj.getFilteredData(MAll))'; %filter class needs unsqueezed input
        MAll=squeeze(MAll)';
        nSamples=size(MAll,1);
        %tAll=((1:nSamples)/obj.dataRecordingObj.samplingFrequency*1000)-obj.postPreRawWindow+startTimes(i);
        for j=1:nNeurons
            tTmp=t(ic(3,j):ic(4,j));
            tTmp=tTmp(find(tTmp>=startTimes(i) & tTmp< (endTimes(i)-obj.postTotalRawWindow+obj.postPreRawWindow) ));
            nSpkTmp=numel(tTmp);
            if nSpkTmp>0
                startIdx=round((tTmp-startTimes(i))*obj.dataRecordingObj.samplingFrequency/1000);
                idx=bsxfun(@plus,startIdx,(0:windowSamplesRaw-1)');
                WF=MAll(idx,:);
                WF=permute(reshape(WF,[size(idx,1) size(idx,2) obj.nCh]),[2 3 1]);
                %WF=bsxfun(@minus,WF,mean(WF(:,:,pPreBaselineSamples),3));
                avgRawWF(j,:,:)=avgRawWF(j,:,:)+sum(WF,1);
                stdRawWF(j,:,:)=stdRawWF(j,:,:)+sum(WF.^2,1);
                
                idx=idx(HPIdx,:);
                WFH=MFAll(idx,:);
                WFH=permute(reshape(WFH,[size(idx,1) size(idx,2) obj.nCh]),[2 3 1]);
                avgHPWF(j,:,:)=avgHPWF(j,:,:)+sum(WFH,1);
                stdHPWF(j,:,:)=stdHPWF(j,:,:)+sum(WFH.^2,1);
                
                for k=1:obj.nCh
                    hVTmp=round((WF(:,k,:)+par.maxV)/(par.maxV*2/nBinsV)); %the reshape is important since otherwise when nSpkTmp==1 squeeze changes matrix from row to column
                    p=find(hVTmp>0 & hVTmp<nBinsV);
                    if numel(p)>1 %if there is only one point to include in the statistics we can reject this  trial.
                        htTmp=repmat(tBinRaw,[numel(tTmp) 1]);
                        if nSpkTmp>1
                            histRawWF(j,k,:,:)=histRawWF(j,k,:,:)+shiftdim(uint8(accumarray([htTmp(p) hVTmp(p)],1,[nBinsRaw nBinsV])),-2);
                        else
                            histRawWF(j,k,:,:)=histRawWF(j,k,:,:)+shiftdim(uint8(accumarray(squeeze([htTmp(p) hVTmp(p)])',1,[nBinsRaw nBinsV])),-2);
                        end
                        
                        hVTmp=round((WFH(:,k,:)+par.maxV)/(par.maxV*2/nBinsV));
                        p=find(hVTmp>0 & hVTmp<nBinsV);
                        htTmp=repmat(tBinHP,[numel(tTmp) 1]);
                        if nSpkTmp>1
                            histHPWF(j,k,:,:)=histHPWF(j,k,:,:)+shiftdim(uint8(accumarray([htTmp(p) hVTmp(p)],1,[nBinsHP nBinsV])),-2);
                        else
                            histHPWF(j,k,:,:)=histHPWF(j,k,:,:)+shiftdim(uint8(accumarray(squeeze([htTmp(p) hVTmp(p)])',1,[nBinsHP nBinsV])),-2);
                        end
                    end
                end
            end
        end
        tChunk(i)=toc;tic;
        fprintf('(%d) ',round(tChunk(i)));
    end
    
    for i=1:nNeurons
        tmpN=(ic(4,i)-ic(3,i)+1);
        stdRawWF(i,:,:)=sqrt((stdRawWF(i,:,:)-avgRawWF(i,:,:))/(tmpN-1));
        stdHPWF(i,:,:)=sqrt((stdHPWF(i,:,:)-avgHPWF(i,:,:))/(tmpN-1));
        
        avgRawWF(i,:,:)=avgRawWF(i,:,:)/tmpN;
        avgHPWF(i,:,:)=avgHPWF(i,:,:)/tmpN;
    end
    
    [mElecs,nElecs]=size(obj.chPar.En);
    if obj.postPlotRawLongWaveforms
        figurePosition=[100 50 min(1000,100*nElecs) min(900,100*mElecs)];
        for i=1:nNeurons
            f=figure('position',figurePosition);
            neuronString=['Neu ' num2str(ic(1,i)) '-' num2str(ic(2,i))];
            infoStr={['nSpk=' num2str(nSpkTotal(i))],neuronString};
            [h,hParent]=spikeDensityPlotPhysicalSpace(squeeze(histRawWF(i,:,:,:)),obj.dataRecordingObj.samplingFrequency(1),obj.chPar.s2r,obj.chPar.En,...
                'hParent',f,'avgSpikeWaveforms',permute(avgRawWF(i,:,:),[3 1 2]),'logColorScale',0,'spikeShapesIsAHist',1,'xBin',tBinsRaw,'yBin',unique(VBins));
            annotation('textbox',[0.01 0.89 0.1 0.1],'FitHeightToText','on','String',infoStr);
            
            printFile=[obj.sortingDir filesep 'neuron' neuronString '-spikeShapeRaw'];
            set(f,'PaperPositionMode','auto');
            print(printFile,'-djpeg','-r300');
            close(f);
        end
    end
    
    if obj.postPlotFilteredWaveforms	
        figurePosition=[100 50 min(1000,100*nElecs) min(900,100*mElecs)];
        for i=1:nNeurons
            f=figure('position',figurePosition);
            neuronString=['Neu ' num2str(ic(1,i)) '-' num2str(ic(2,i))];
            infoStr={['nSpk=' num2str(nSpkTotal(i))],neuronString};
            [h,hParent]=spikeDensityPlotPhysicalSpace(squeeze(histHPWF(i,:,:,:)),obj.dataRecordingObj.samplingFrequency(1),obj.chPar.s2r,obj.chPar.En,...
                'hParent',f,'avgSpikeWaveforms',permute(avgHPWF(i,:,:),[3 1 2]),'logColorScale',0,'spikeShapesIsAHist',1,'xBin',tBinsHP,'yBin',unique(VBins));
            annotation('textbox',[0.01 0.89 0.1 0.1],'FitHeightToText','on','String',infoStr);
            
            printFile=[obj.sortingDir filesep 'neuron' neuronString '-spikeShapeHP'];
            set(f,'PaperPositionMode','auto');
            print(printFile,'-djpeg','-r300');
            close(f);
        end
    end
    
catch errorMsg
    if i>1
        lastGoodChunck=i-1;
        save([obj.sortingDir filesep saveFileName '_tmp.mat'],'avgRawWF','avgHPWF','stdRawWF','stdHPWF','histRawWF','histHPWF','par','lastGoodChunck','-v7.3');
    end
    rethrow(errorMsg);
end
fprintf('\nTriggered waveform analysis took (%f) hours\nSaving data...\n',sum(tChunk)/60/60);

save(obj.sortingFileNames.STWaveformFile,'avgRawWF','avgHPWF','stdRawWF','stdHPWF','neuronNames','histRawWF','histHPWF','nSpkTotal','par','-v7.3');
fprintf('Deleting temporary files...');
delete([obj.sortingDir filesep saveFileName '_tmp.mat']);
obj=obj.findSortingFiles; %update sorted files
