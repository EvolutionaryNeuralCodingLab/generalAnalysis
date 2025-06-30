function data=getMultiBandSpectralAnalysis(obj,varargin)
obj.checkFileRecording;

parseObj = inputParser;
addParameter(parseObj,'ch',obj.recTable.defaulLFPCh(obj.currentPRec),@isnumeric);
addParameter(parseObj,'avgOnCh',[],@isnumeric); %uses several averaged channels for d/b extraction
addParameter(parseObj,'movLongWin',1000*60*30,@isnumeric); %max freq. to examine
addParameter(parseObj,'movWin',10000,@isnumeric);
addParameter(parseObj,'movOLWin',9000,@isnumeric);
addParameter(parseObj,'segmentWelch',1000,@isnumeric);
addParameter(parseObj,'dftPointsWelch',2^10,@isnumeric);
addParameter(parseObj,'OLWelch',0.5);
addParameter(parseObj,'tStart',0,@isnumeric);
addParameter(parseObj,'win',0,@isnumeric); %if 0 uses the whole recording duration
addParameter(parseObj,'band1Low',5,@isnumeric);
addParameter(parseObj,'band1High',30,@isnumeric);
addParameter(parseObj,'band2Low',35,@isnumeric);
addParameter(parseObj,'band2High',60,@isnumeric);
addParameter(parseObj,'applyNotch',0,@isnumeric);
addParameter(parseObj,'saveSpectralProfiles',0,@isnumeric);
addParameter(parseObj,'maxVoltage',1000,@isnumeric);
addParameter(parseObj,'overwrite',0,@isnumeric);
addParameter(parseObj,'inputParams',false,@isnumeric);
parseObj.parse(varargin{:});
if parseObj.Results.inputParams
    disp(parseObj.Results);
    return;
end

%evaluate all input parameters in workspace
for i=1:numel(parseObj.Parameters)
    eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
end

%make parameter structure
par=parseObj.Results;

if isnan(ch)
    disp('Error: no reference channel for multiBand extraction');
    return;
end
%check if analysis was already done done
obj.files.multiBand=sprintf('%smultiBand_ch%d_1B%d_%d_2B%d_%d.mat',[obj.currentAnalysisFolder filesep],ch,band1Low,band1High,band2Low,band2High);
if exist(obj.files.multiBand,'file') & ~overwrite
    if nargout==1
        data=load(obj.files.multiBand);
    else
        disp('MultiBand analysis already exists for this recording');
    end
    return;
end

obj.getFilters;
movWinSamples=movWin/1000*obj.filt.FFs;%obj.filt.FFs in Hz, movWin in samples
movOLWinSamples=movOLWin/1000*obj.filt.FFs;
timeBin=(movWin-movOLWin); %ms

segmentWelchSamples = round(segmentWelch/1000*obj.filt.FFs);
samplesOLWelch = round(segmentWelchSamples*OLWelch);

%run welch once to get frequencies for every bin (f) determine frequency bands
[~,f] = pwelch(randn(1,movWinSamples),segmentWelchSamples,samplesOLWelch,dftPointsWelch,obj.filt.FFs);
pBand1=find(f>=band1Low & f<band1High);
pBand2=find(f>=band2Low & f<band2High);

%if obj.currentDataObj.recordingDuration_ms<movLongWin
%    movLongWin=obj.currentDataObj.recordingDuration_ms;
%end

if win==0
    win=obj.currentDataObj.recordingDuration_ms-tStart;
    endTime=obj.currentDataObj.recordingDuration_ms;
else
    endTime=min(win+tStart,obj.currentDataObj.recordingDuration_ms);
end
startTimes=tStart:(movLongWin-movOLWin):endTime;
nChunks=numel(startTimes);
band2to1RatioAll=cell(1,nChunks);
t_ms=cell(1,nChunks);
%band2to1RatioAllLow=cell(1,nChunks);;band2to1RatioAllHigh=cell(1,nChunks);

if saveSpectralProfiles
    FMLongB = buffer(true(1,movLongWin/1000*obj.filt.FFs),movWinSamples,movOLWinSamples,'nodelay');
    fftInBuffer=size(FMLongB,2);
    allFreqProfiles=zeros(ceil(dftPointsWelch/2)+1,nChunks*fftInBuffer);
else
    allFreqProfiles=[];
end
if applyNotch
    obj.filt.FN=filterData(obj.currentDataObj.samplingFrequency(1));
    obj.filt.FN.filterDesign='cheby1';
    obj.filt.FN.padding=true;
    obj.filt.FN=obj.filt.FN.designNotch;
end

loadCh=ch;
if ~isempty(avgOnCh)
    loadCh=avgOnCh; %this can not be moved to other positions
end

fprintf('\nDelta2Beta extraction (%d chunks)-',nChunks);
for i=1:nChunks
    fprintf('%d,',i);
    MLong=obj.currentDataObj.getData(loadCh,startTimes(i),movLongWin);
    if applyNotch
        MLong=obj.filt.FN.getFilteredData(MLong); %for 50Hz noise
    end
    FMLong=obj.filt.F.getFilteredData(MLong);
    if ~isempty(avgOnCh)
        FMLong=mean(FMLong,1);
    end

    FMLong(FMLong<-maxVoltage | FMLong>maxVoltage)=nan; %remove high voltage movement artifacts

    FMLongB = buffer(FMLong,movWinSamples,movOLWinSamples,"nodelay");
    pValid=all(~isnan(FMLongB));

    band2to1RatioAll{i}=nan(1,numel(pValid)); %changes from zeros to nan in these 3 lines (Mark)
    band1All{i}=nan(1,numel(pValid));
    band2All{i}=nan(1,numel(pValid));
    if any(pValid)
        [pxx,f] = pwelch(FMLongB(:,pValid),segmentWelchSamples,samplesOLWelch,dftPointsWelch,obj.filt.FFs);
        band2to1RatioAll{i}(pValid)=(mean(pxx(pBand2,:))./mean(pxx(pBand1,:)))';
        band1All{i}(pValid)=mean(pxx(pBand1,:))';
        band2All{i}(pValid)=mean(pxx(pBand2,:))';
    else
        pxx=zeros(dftPointsWelch/2+1,numel(pValid));
    end

    if saveSpectralProfiles
        allFreqProfiles(:,(fftInBuffer*(i-1)+find(pValid)))=pxx;
    end

    t_ms{i}=startTimes(i)+((movWin/2):timeBin:(movLongWin-movWin/2));
end
fprintf('\n');
band2to1RatioAll{end}(t_ms{end}>(endTime-movWin/2))=NaN;
band1All{end}(t_ms{end}>(endTime-movWin/2))=NaN;
band2All{end}(t_ms{end}>(endTime-movWin/2))=NaN;

band2to1Ratio=cell2mat(band2to1RatioAll);band2to1Ratio=band2to1Ratio(:);
band1=cell2mat(band1All);band1=band1(:);
band2=cell2mat(band2All);band2=band2(:);

t_ms=cell2mat(t_ms);

save(obj.files.multiBand,'t_ms','band2to1Ratio','par','band2','band1','allFreqProfiles');