function getSpikeTrigWF(obj,startEnd,printFolder)
%Detect spikes in raw data and save waveforms
%Todo:
%Add detection of multiple spikes
if nargin<3
    printFolder=obj.sortingDir;
end

%determine the chunck size

if nargin>1
    if ~isempty(startEnd)
        disp(['Plotting post processing on a subset of data ' num2str(round(startEnd(1)/1000)) ' - ' num2str(round(startEnd(2)/1000)) 's']);
        postFix=[num2str(round(startEnd(1)/1000)) '_' num2str(round(startEnd(2)/1000)) 's'];
    else
        startEnd=[0 obj.dataRecordingObj.recordingDuration_ms];
        postFix='';
    end
else
    startEnd=[0 obj.dataRecordingObj.recordingDuration_ms];
    postFix='';
end

load([obj.sortingFileNames.STWaveformFile(1:end-4) postFix],'avgRawWF','avgHPWF','neuronNames','histRawWF','histHPWF','nSpkTotal','par');
nNeurons=size(neuronNames,2);

windowSamplesRaw=obj.postTotalRawWindow*obj.dataRecordingObj.samplingFrequency(1)/1000;
windowSamplesFiltered=obj.postTotalFilteredWindow*obj.dataRecordingObj.samplingFrequency(1)/1000;
tBinRaw=shiftdim(ceil(((1:windowSamplesRaw)/obj.dataRecordingObj.samplingFrequency(1)*1000)/par.binTRaw),-1);
tBinHP=shiftdim(ceil(((1:windowSamplesFiltered)/obj.dataRecordingObj.samplingFrequency(1)*1000)/par.binTHP),-1);
VBins=(-par.maxV+par.binV/2):par.binV:par.maxV;
tBinsRaw=par.binTRaw/2:par.binTRaw:obj.postTotalRawWindow;
tBinsHP=par.binTHP/2:par.binTHP:obj.postTotalFilteredWindow;

[mElecs,nElecs]=size(obj.chPar.En);
if obj.postPlotRawLongWaveforms
    figurePosition=[100 50 min(1000,100*nElecs) min(900,100*mElecs)];
    for i=1:nNeurons
        f=figure('position',figurePosition);
        neuronString=['Neu' num2str(neuronNames(1,i)) '-' num2str(neuronNames(2,i))];
        infoStr={['nSpk=' num2str(nSpkTotal(i))],neuronString};
        [h,hParent]=spikeDensityPlotPhysicalSpace(squeeze(histRawWF(i,:,:,:)),obj.dataRecordingObj.samplingFrequency(1),obj.chPar.s2r,obj.chPar.En,...
            'hParent',f,'avgSpikeWaveforms',permute(avgRawWF(i,:,:),[3 1 2]),'logColorScale',1,'spikeShapesIsAHist',1,'xBin',tBinsRaw,'yBin',unique(VBins));
        annotation('textbox',[0.01 0.89 0.1 0.1],'FitHeightToText','on','String',infoStr);
        
        printFile=[printFolder filesep 'neuron' neuronString '-spikeShapeRaw' postFix];
        set(f,'PaperPositionMode','auto');
        print(printFile,'-djpeg','-r300');
        close(f);
    end
end

if obj.postPlotFilteredWaveforms
    figurePosition=[100 50 min(1000,100*nElecs) min(900,100*mElecs)];
    for i=1:nNeurons
        f=figure('position',figurePosition);
        neuronString=['Neu' num2str(neuronNames(1,i)) '-' num2str(neuronNames(2,i))];
        infoStr={['nSpk=' num2str(nSpkTotal(i))],neuronString};
        [h,hParent]=spikeDensityPlotPhysicalSpace(squeeze(histHPWF(i,:,:,:)),obj.dataRecordingObj.samplingFrequency(1),obj.chPar.s2r,obj.chPar.En,...
            'hParent',f,'avgSpikeWaveforms',permute(avgHPWF(i,:,:),[3 1 2]),'logColorScale',1,'spikeShapesIsAHist',1,'xBin',tBinsHP,'yBin',unique(VBins));
        annotation('textbox',[0.01 0.89 0.1 0.1],'FitHeightToText','on','String',infoStr);
        
        printFile=[printFolder filesep 'neuron' neuronString '-spikeShapeHP' postFix];
        set(f,'PaperPositionMode','auto');
        print(printFile,'-djpeg','-r300');
        close(f);
    end
end