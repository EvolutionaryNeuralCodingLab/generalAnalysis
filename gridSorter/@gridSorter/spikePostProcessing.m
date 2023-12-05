function obj=spikePostProcessing(obj,tStartEnd)
% Run post-processing steps on spike data. Calculate post-spike fields, and spike neuron identification reliability
% obj=spikePostProcessingNSK(obj,tStart,tEnd)
% Method recieves : gridSorter object
%                   tStartEnd - [1 x 2] the start and end of a desired segment of the recording one which post processing is calculated
%                               Default - uses the t,ic of the whole recording
%                 
% Function give back
% Last updated : 28/03/15
  
load(obj.sortingFileNames.fittingFile,'t','ic','md','allWaveforms','isNoiseAll');

if nargin==2
    disp(['Calculating post processing on a subset of data ' num2str(tStartEnd(1)/1000/60/60) ' - ' num2str(tStartEnd(2)/1000/60/60) 'h']);
    [t,ic]=CutSortChannel(t,ic,tStartEnd(1),tStartEnd(2),false);
    postFix=[num2str(tStartEnd(1)/1000/60/60,3) '_' num2str(tStartEnd(2)/1000/60/60,3)];
    obj.sortingFileNames.postProcessingAnalysisExist=0; %overwriting
else
    postFix=[];
end

if obj.postPlotRawLongWaveforms || obj.postPlotFilteredWaveforms %calculate figure aspect ration based on electrode layout
    [mElecs,nElecs]=size(obj.chPar.En);
    figurePosition=[100 50 min(1000,100*nElecs) min(900,100*mElecs)];
end
fprintf('\nPost processing spikes...');

[nSamples,nNeurons,obj.nCh]=size(allWaveforms);

if ~isempty(ic)
    neuronNames=ic(1:2,:);
else
    neuronNames=[];
end

if isempty(obj.filterObj)
    obj=getHighpassFilter(obj);
end

%for fast plotting do not pad filtered waveforms (this increases the speed considerably since padding requires about 20 sec!!!!!!!!!!!
obj.filterObj.padding=false;

nonExistingNeurons=find(~obj.sortingFileNames.postProcessingAnalysisExist | obj.overwritePostProcessingAnalysis);
nNonExistingNeurons=numel(nonExistingNeurons);

matFileObj = matfile([obj.sortingFileNames.postProcessingAnalysisFile(1:end-4) postFix '.mat'],'Writable',true);
if nNonExistingNeurons==1 %there is no matfile
    matFileObj.postProcessingAnalysisExist=zeros(1,nNeurons);
    nonExistingNeurons=1:nNeurons;
    nNonExistingNeurons=numel(nonExistingNeurons);
    matFileObj.neuronNames=neuronNames;
    if obj.postExtractRawLongWaveformsFromSpikeTimes
        matFileObj.avgLongWF=zeros(nNeurons,obj.nCh,round(obj.postTotalRawWindow*obj.dataRecordingObj.samplingFrequency(1)/1000));
        matFileObj.stdLongWF=zeros(nNeurons,obj.nCh,round(obj.postTotalRawWindow*obj.dataRecordingObj.samplingFrequency(1)/1000));
        matFileObj.PSDSNR=zeros(1,nNeurons);
    else
        matFileObj.avgLongWF=[];matFileObj.stdLongWF=[];matFileObj.PSDSNR=[];
    end
    
    if obj.postExtractFilteredWaveformsFromSpikeTimes
        matFileObj.avgFinalWF=zeros(nNeurons,obj.nCh,round(obj.postTotalFilteredWindow*obj.dataRecordingObj.samplingFrequency(1)/1000));
        matFileObj.stdFinalWF=zeros(nNeurons,obj.nCh,round(obj.postTotalFilteredWindow*obj.dataRecordingObj.samplingFrequency(1)/1000));
        matFileObj.spkSNR=zeros(1,nNeurons);
        matFileObj.nSpkTotal=zeros(1,nNeurons);
        matFileObj.spkMaxAmp=zeros(1,nNeurons);
    else
        matFileObj.avgFinalWF=[];matFileObj.stdFinalWF=[];matFileObj.spkSNR=[];matFileObj.nSpkTotal=[];matFileObj.spkMaxAmp=[];
    end
end

fprintf('Calculating SIFs (nNeurons=%d)...',numel(nonExistingNeurons));
tic;
for i=nonExistingNeurons
    tInter=toc;
    fprintf('(%d) %d-%d, ',round(tInter/60),ic(1,i),ic(2,i));
    tic;
    tTmp=t(ic(3,i):ic(4,i));
    mdTmp=md(ic(3,i):ic(4,i));
    nSpk=numel(tTmp);
    if obj.postExtractRawLongWaveformsFromSpikeTimes
        [~,pOrder]=sort(mdTmp);
        tSpk=sort(tTmp(pOrder(1:min(nSpk,obj.postMaxSpikes2Present)))); %the the best spikes for SIFs according to lowest minkowski distance
        
        [V_uV,T_ms]=obj.dataRecordingObj.getData(obj.chPar.s2r(1:obj.nCh),tSpk-obj.postPreRawWindow,obj.postTotalRawWindow);
        
        %standard deviation
        stdLongWF(1,:,:)=squeeze(std(V_uV,[],2));
        
        %average substruction mean
        avgLongWF(1,:,:)=squeeze(median(V_uV,2));
        tmpMean=mean(avgLongWF(1,:,T_ms<(obj.postPreRawWindow-1)),3); %calculate the baseline according to the average voltage before the spike (1ms before spike peak)
        V_uV=bsxfun(@minus,V_uV,tmpMean');
        avgLongWF(1,:,:)=bsxfun(@minus,avgLongWF(1,:,:),tmpMean);
        
        p=find(T_ms>obj.postRawSNRStartEnd(1) & T_ms<=obj.postRawSNRStartEnd(2)); %find relevant time points surrounding the spike
        tmpSNR=avgLongWF(1,obj.chPar.surChExtVec{obj.chPar.r2s(ic(1,i))},p)./stdLongWF(1,obj.chPar.surChExtVec{obj.chPar.r2s(ic(1,i))},p); %calculated at the position of the extended grid
        
        if obj.postPlotRawLongWaveforms
            
            f=figure('position',figurePosition);
            neuronString=['Neu ' num2str(ic(1,i)) '-' num2str(ic(2,i))];
            infoStr={['nSpk=' num2str(nSpk)],neuronString,['noise=' num2str(isNoiseAll(i))]};
            
            [h,hParent]=spikeDensityPlotPhysicalSpace(permute(V_uV,[3 2 1]),obj.dataRecordingObj.samplingFrequency(1),obj.chPar.s2r,obj.chPar.En,...
                'hParent',f,'avgSpikeWaveforms',permute(avgLongWF(1,:,:),[3 1 2]),'logDensity',true);
            
            annotation('textbox',[0.01 0.89 0.1 0.1],'FitHeightToText','on','String',infoStr);
            
            printFile=[obj.sortingDir filesep 'neuron' neuronString postFix '-spikeShapeRaw'];
            
            set(f,'PaperPositionMode','auto');
            print(printFile,'-djpeg','-r300');
            close(f);
        end
        
        matFileObj.PSDSNR(1,i)=mean(abs(tmpSNR(:)));
        matFileObj.stdLongWF(i,:,:)=stdLongWF;
        matFileObj.avgLongWF(i,:,:)=avgLongWF;
    end
    
    if obj.postExtractFilteredWaveformsFromSpikeTimes
        if obj.postExtractRawLongWaveformsFromSpikeTimes
            [V_uV]=obj.filterObj.getFilteredData(V_uV( : , : , :));
            V_uV=V_uV(:,:,T_ms>=(obj.postPreRawWindow-obj.postPreFilteredWindow) & T_ms<(obj.postPreRawWindow-obj.postPreFilteredWindow+obj.postTotalFilteredWindow));
            T_ms=T_ms(T_ms>=(obj.postPreRawWindow-obj.postPreFilteredWindow) & T_ms<(obj.postPreRawWindow-obj.postPreFilteredWindow+obj.postTotalFilteredWindow));
            %[V_uV,T_ms]=V_uV( : , : , T_ms>=(obj.postPreRawWindow-obj.postPreFilteredWindow) & T_ms<(obj.postPreRawWindow-obj.postPreFilteredWindow+obj.postTotalFilteredWindow));
        else
            %[V_uV,T_ms]=obj.filterObj.getFilteredData(obj.dataRecordingObj.getData(obj.chPar.s2r(1:obj.nCh),tTmp(1:min(nSpk,obj.postMaxSpikes2Present))-obj.postPreFilteredWindow,obj.postTotalFilteredWindow));
            [V_uV,T_ms]=obj.dataRecordingObj.getData(obj.chPar.s2r(1:obj.nCh),tTmp(1:min(nSpk,obj.postMaxSpikes2Present))-obj.postPreFilteredWindow,obj.postTotalFilteredWindow);
        end
        
        avgFinalWF(1,:,:)=squeeze(median(V_uV,2));
        stdFinalWF(1,:,:)=squeeze(std(V_uV,[],2));
        
        p=find(T_ms>obj.postFilteredSNRStartEnd(1) & T_ms<=obj.postFilteredSNRStartEnd(2)); %find relevant time points surrounding the spike
        tmpSNR=avgFinalWF(1,obj.chPar.surChExtVec{obj.chPar.r2s(ic(1,i))}(obj.chPar.pSurCh{obj.chPar.r2s(ic(1,i))}),p)./...
            stdFinalWF(1,obj.chPar.surChExtVec{obj.chPar.r2s(ic(1,i))}(obj.chPar.pSurCh{obj.chPar.r2s(ic(1,i))}),p); %calculated at the positions of the surrounding grid
        
        if obj.postPlotFilteredWaveforms
            
            f=figure('position',figurePosition);
            neuronString=['Neu ' num2str(ic(1,i)) '-' num2str(ic(2,i))];
            infoStr={['nSpk=' num2str(nSpk)],neuronString,['noise=' num2str(isNoiseAll(i))]};
            
            [h,hParent]=spikeDensityPlotPhysicalSpace(permute(V_uV,[3 2 1]),obj.dataRecordingObj.samplingFrequency(1),obj.chPar.s2r,obj.chPar.En,...
                'hParent',f,'avgSpikeWaveforms',permute(avgFinalWF(1,:,:),[3 1 2]),'logDensity',true);
            
            annotation('textbox',[0.01 0.89 0.1 0.1],'FitHeightToText','on','String',infoStr);
            
            %print to file
            printFile=[obj.sortingDir filesep 'neuron' neuronString postFix '-spikeShape'];
            set(f,'PaperPositionMode','auto');
            print(printFile,'-djpeg','-r300');
            close(f);
        end
        
        matFileObj.avgFinalWF(i,:,:)=avgFinalWF;
        matFileObj.stdFinalWF(i,:,:)=stdFinalWF;
        matFileObj.spkSNR(1,i)=mean(abs(tmpSNR(:)));
        matFileObj.spkMaxAmp(1,i)=max(max(abs(avgFinalWF))); %the extremal spike amplitude
        matFileObj.nSpkTotal(1,i)=ic(4,i)-ic(3,i)+1;
        
    end
    matFileObj.postProcessingAnalysisExist(1,i)=1;
end %loop over all neurons
fprintf('Done!\n');

if obj.postPlotAllAvgSpikeTemplates
    nPlotAxis=min(ceil(sqrt(nNeurons)),5);
    nPlotsPerFigure=nPlotAxis.^2;
    if nNeurons>0
        for i=1:(nNeurons+1)
            if mod((i-1),nPlotsPerFigure)==0 || i==(nNeurons+1)
                if i~=1
                    printFile=[obj.sortingDir filesep 'avgSpikeShapes' postFix num2str(ceil((i-1)/nPlotsPerFigure))];
                    if obj.fastPrinting
                        imwrite(frame2im(getframe(f)),[printFile '.jpeg'],'Quality',90);
                    else
                        set(f,'PaperPositionMode','auto');
                        print(printFile,'-djpeg','-r300');
                    end
                    close(f);
                end
                if i~=(nNeurons+1)
                    f=figure('position',figurePosition);
                end
            end
            if i<nNeurons
                neuronString=['Neu ' num2str(ic(1,i)) '-' num2str(ic(2,i)),',N' num2str(isNoiseAll(i))];
                h=subaxis(f,nPlotAxis,nPlotAxis,mod(i,nPlotsPerFigure),'S',0.02,'M',0.02);
                activityTracePhysicalSpacePlot(h,obj.chPar.s2r,0.03*squeeze(allWaveforms(:,i,:))',obj.chPar.En,'scaling','none');
                text(0,0,neuronString);
            elseif i==nNeurons
                neuronString=['Neu ' num2str(ic(1,i)) '-' num2str(ic(2,i)),',N' num2str(isNoiseAll(i))];
                h=subaxis(f,nPlotAxis,nPlotAxis,nNeurons,'S',0.02,'M',0.02);
                activityTracePhysicalSpacePlot(h,obj.chPar.s2r,0.03*squeeze(allWaveforms(:,i,:))',obj.chPar.En,'scaling','none');
                text(0,0,neuronString);
            end
        end
    end
end

if obj.postPlotSpikeReliability
    if nNeurons>0
        f=figure('position',[300 50 900 700]);
        mNSpkTotal=mean(matFileObj.nSpkTotal);
        %generate size legend
        plotSpikeSNR=[matFileObj.spkSNR (max(matFileObj.spkSNR)-(max(matFileObj.spkSNR)-min(matFileObj.spkSNR))*0.01)*ones(1,5)];
        plotPSDSNR=[matFileObj.PSDSNR [0.9 1 1.1 1.2 1.3]*((max(matFileObj.PSDSNR)+min(matFileObj.PSDSNR))/2)];
        legendSpikeNums=round([mNSpkTotal/6 mNSpkTotal/3 mNSpkTotal mNSpkTotal*3 mNSpkTotal*6]);
        plotnSpkTotal=[matFileObj.nSpkTotal legendSpikeNums];
        plotspkMaxAmp=[matFileObj.spkMaxAmp ones(1,5)*min(matFileObj.spkMaxAmp)];
        
        scatter(plotSpikeSNR,plotPSDSNR,(plotnSpkTotal/mNSpkTotal)*50+5,plotspkMaxAmp,'linewidth',2);
        text(plotSpikeSNR(end-4:end)*1.02,plotPSDSNR(end-4:end),num2str(legendSpikeNums'/(obj.dataRecordingObj.recordingDuration_ms/1000),3),'FontSize',8)
        xlabel('$$\sqrt{SNR_{spike}}$$','Interpreter','latex','FontSize',14);ylabel('$$\sqrt{SNR_{PSD}}$$','Interpreter','latex','FontSize',14);
        cb=colorbar('position',[0.8511    0.6857    0.0100    0.2100]);ylabel(cb,'Max spk. amp.');
        
        printFile=[obj.sortingDir filesep 'SNR_spikePSD' postFix];
        set(f,'PaperPositionMode','auto');
        print(printFile,'-djpeg','-r300');
    end
end

%[obj, er]=assessQuality(obj);

obj=obj.findSortingFiles; %update sorted files
