function results = CalculateReceptiveFields(obj,params)

arguments (Input)
    obj
    params.overwrite logical = false
    params.analysisTime = datetime('now')
    params.inputParams = false
    params.preBase = 200
    params.bin = 10
    params.exNeurons = 0
    params.AllResponsiveNeurons = true
    params.fixedWindow = false
    params.noEyeMoves = false
    params.delay = 250
    params.nShuffle = 20 %Number of shuffles to generate shuffled receptive fields. 
    params.testConvolution = false
    params.reduceFactor = 20;
    params.duration = 300; %response window
    params.offsetR = 50; %Response after onset of stim
end


if params.inputParams,disp(params),return,end

filename = obj.getAnalysisFileName;


if isfile(obj.getAnalysisFileName) && ~params.overwrite
    if nargout==1
        fprintf('Loading saved results from file.\n');
        results=load(filename);
    else
        fprintf('Analysis already exists (use overwrite option to recalculate).\n');
    end

    return
end

NeuronResp = obj.ResponseWindow;
Stats = obj.ShufflingAnalysis;
goodU = NeuronResp.goodU;
p = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
phy_IDg = p.phy_ID(string(p.label') == 'good');

pvals = Stats.pvalsResponse;
C = NeuronResp.C;
stimDur = NeuronResp.stimDur;

if params.AllResponsiveNeurons
    respU = find(pvals<0.05);
    if isempty(respU)
        fprintf('No responsive neurons.\n')
        return
    end
end

if params.exNeurons >0
    respU = params.exNeurons;
end

seqMatrix = obj.VST.pos;
sizes = obj.VST.tilingRatios;
uSize = unique(sizes);
nSize = length(uSize);
uLums = unique(obj.VST.rectLuminosity(obj.VST.luminosities));
nLums = length(uLums); %%mAKE IT TO BE ABLE TO COMPARE TWO LUMINOSITIES. 

trialDiv  = length(seqMatrix)/length(unique(seqMatrix))/nSize/nLums;
directimesSorted = C(:,1)';

[Mr] = BuildBurstMatrix(goodU,round(p.t/params.bin),round((directimesSorted+params.offsetR)/params.bin),round(params.duration/params.bin));
[Mro] = BuildBurstMatrix(goodU,round(p.t/params.bin),round((directimesSorted+stimDur)/params.bin),round(params.duration/params.bin));
[Mb1] = BuildBurstMatrix(goodU,round(p.t/params.bin),round((directimesSorted-params.duration)/params.bin),round(params.duration/params.bin));

[nT,nN,NB] = size(Mr);
[nTo,nNo,NBo] = size(Mro);


%%%%%%%%%%%%%%%%%%%% Shuffle raster before point multiplication in order
%%%%%%%%%%%%%%%%%%%% to calculate tuning index
%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%

nShuffle =params.nShuffle;

Raster = BuildBurstMatrix(goodU,round(p.t/params.bin),round((directimesSorted)/params.bin),round((stimDur)/params.bin));

shuffledData = zeros(size(Raster,1), size(Raster,2), size(Raster,3), nShuffle);

for i =1:nShuffle

    % Shuffle along the first dimension
    idx1 = randperm(size(Raster,1));

    % Shuffle along the third dimension
    idx3 = randperm(size(Raster,3));

    shuffledData(:,:,:,i) = Raster(idx1, :, idx3);

end

if params.noEyeMoves
    % EyePositionAnalysis
    % Create spike Sums with NaNs when the eye is not present.
    % 
    % 
    % file = dir (NP.recordingDir);
    % filenames = {file.name};
    % files= filenames(contains(filenames,"timeSnipsNoMov-31"));
    % cd(NP.recordingDir)
    % %Run eyePosition Analysis to find no movement timeSnips
    % timeSnips = load(files{1}).timeSnips;
    % timeSnipsMode = timeSnips(:,timeSnips(3,:) == mode(timeSnips(3,:)));
    % 
    % selecInd = [];
    % for i = 1:size(timeSnipsMode,2)
    % 
    %     %Find stimOns and offs that are between each timeSnip
    %     selecInd = [selecInd find(directimesSorted>=timeSnipsMode(1,i) & directimesSorted<(timeSnipsMode(2,i)-stimDur))];
    % end
    % 
    % %MrC = nan(round(nT/trialDiv),nN, NB+NBo);
    % 
    % MrC = nan(round(nT/trialDiv),nN, NB);
    % 
    % 
    % %%Create summary of identical trials
    % 
    % MrMean = nan(round(nT/trialDiv),nN);
    % 
    % for u = 1:length(goodU)
    %     j=1;
    % 
    %     for i = 1:trialDiv:nT
    % 
    %         indexVal = selecInd(selecInd>=i & selecInd<=i+trialDiv-1);
    % 
    %         if ~isempty(indexVal)
    % 
    % 
    %             meanRon =  reshape(mean(Mr(indexVal,u,:),1),[1,size(Mr,3)]);
    % 
    %             meanRoff =  reshape(mean(Mro(indexVal,u,:),1),[1,size(Mro,3)]);
    % 
    %             %meanBase =  reshape(mean(Mb1(indexVal,u,:),1),[1,size(Mb1,3)]);
    % 
    %             %MrC(j,u,:) = [meanRon-meanBase meanRoff-meanBase]; %Combine on and off response and substract to each the mean baseline
    % 
    %             MrC(j,u,:) = [meanRon];
    %             MrMean(j,u) = mean(MrC(j,u,:),3);%-Nbase;
    % 
    %         else
    %             2+2
    %         end
    % 
    % 
    %         j = j+1;
    % 
    %     end
    % end
    % 
    % 2+2

else

    MrC = zeros(2,nLums,nSize,round(nT/trialDiv),nN, NB);

    MRtotal = zeros(2,size(Mr,1),size(Mr,2),size(Mr,3)); %includes on and off response

    MRtotal(1,:,:,:) =  Mr;

    MRtotal(2,:,:,:) =  Mro;

    %%Create summary of identical trials

    for u = 1:length(goodU)
       

        for o = 1:2
             j=1;
            for i = 1:trialDiv:nT

                meanR = mean(squeeze(MRtotal(o,i:i+trialDiv-1,u,:)));

                MrC(o,uLums == C(i,4),uSize == C(i,3),j,u,:) =meanR; %Combine on and off response

                j = j+1;
            end
        end
    end

    MrMean = mean(MrC,6);%-Nbase;

end



screenSide = obj.VST.rect; %Same as moving ball

screenRed = screenSide(4)/params.reduceFactor;
[x, y] = meshgrid(1:screenRed, 1:screenRed);

pxyScreen = zeros(screenRed,screenRed);

VideoScreen = zeros(screenRed,screenRed,size(C,1)/trialDiv);

rectData = obj.VST.rectData;

j=1;

for i = 1:trialDiv:length(C)

    xyScreen = zeros(screenRed,screenRed)'; %%Make calculations if sizes>1 and if experiment is new and the shape is a circle.

    % string(obj.VST.shape) == "circle"  %%%Asumes that shape is circle

    Xc = round((rectData.X2{1,C(i,3)}(C(i,2))-rectData.X1{1,C(i,3)}(C(i,2)))/2)+rectData.X1{1,C(i,3)}(C(i,2));%...
    Xc = Xc/params.reduceFactor;

    Yc = round((rectData.Y4{1,C(i,3)}(C(i,2))-rectData.Y1{1,C(i,3)}(C(i,2)))/2)+rectData.Y1{1,C(i,3)}(C(i,2));%...
    Yc = Yc/params.reduceFactor;

    r = round((rectData.X2{1,C(i,3)}(C(i,2))-rectData.X1{1,C(i,3)}(C(i,2)))/2);
    r= r/params.reduceFactor;

    % Calculate the distance of each point from the center
    distances = sqrt((x - Xc).^2 + (y - Yc).^2);

    % Set the values inside the circle to 1 (or any other value you prefer)
    xyScreen(distances <= r) = 1;

    %
    %                 hold on; plot(rect.VSMetaData.allPropVal{21,1}.X1{1,C(i,3)}(C(i,2)),rect.VSMetaData.allPropVal{21,1}.Y1{1,C(i,3)}(C(i,2)),points{C(i,3)},MarkerSize=10)%,...
    %                 hold on; plot(rect.VSMetaData.allPropVal{21,1}.X2{1,C(i,3)}(C(i,2)),rect.VSMetaData.allPropVal{21,1}.Y2{1,C(i,3)}(C(i,2)),points{C(i,3)},MarkerSize=10)
    %                 hold on; plot(rect.VSMetaData.allPropVal{21,1}.X3{1,C(i,3)}(C(i,2)),rect.VSMetaData.allPropVal{21,1}.Y3{1,C(i,3)}(C(i,2)),points{C(i,3)},MarkerSize=10)%,...
    %                 hold on; plot(rect.VSMetaData.allPropVal{21,1}.X4{1,C(i,3)}(C(i,2)),rect.VSMetaData.allPropVal{21,1}.Y4{1,C(i,3)}(C(i,2)),points{C(i,3)},MarkerSize=10)%,...
    %                 hold on; plot(Xc,Yc,points{C(i,3)},MarkerSize=10);

    %figure;imagesc(xyScreen')

    VideoScreen(:,:,j) = xyScreen';

    pxyScreen = pxyScreen+xyScreen;

    j = j+1;

end

%  M = MrMean(:,u)'./Nbase(u);

VD = reshape(VideoScreen,[1 1 1 size(VideoScreen,1) size(VideoScreen,1) size(VideoScreen,3)]);
VD = repmat(VD,[1,1,1,1,1,1,nN]);


NanPos = isnan(MrMean);

MrMean(NanPos) = 0;

Res = reshape(MrMean,[size(MrMean,1),size(MrMean,2),size(MrMean,3),1,1,size(MrMean,4),nN]).*1000;

RFu = reshape(sum(VD.*Res,6),[size(MrMean,1),size(MrMean,2),size(MrMean,3),size(VD,4),size(VD,4),nN]);

offsetN = sqrt(max(seqMatrix));

TwoDGaussian = fspecial('gaussian',floor(size(RFu,4)/(offsetN/2)),screenRed/offsetN);

RFuFilt = zeros(size(RFu));


for d = 1:size(RFu,1) %On off response
    for s = 1:size(RFu,2) %Lums
        for l = 1:size(RFu,3) %size
            for ui =1:size(RFu,6) %units
                slice = squeeze(RFu(d,s,l,:,:,ui));

                slicek = conv2(slice,TwoDGaussian,'same');

                RFuFilt(d,s,l,:,:,ui) =slicek;
            end
        end
    end
end

% figure;imagesc(squeeze(RFu(2,:,:,:,:,83)));
S.RFu = RFu;

S.RFuFilt = RFuFilt;

S.params = params;

save(filename,'-struct','S');

end

