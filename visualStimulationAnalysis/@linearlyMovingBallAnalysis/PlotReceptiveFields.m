function [colorbarLims] = PlotReceptiveFields(obj,params)

arguments (Input)
    obj
    params.overwrite logical = false
    params.analysisTime = datetime('now')
    params.inputParams = false
    params.exNeurons = nan;
    params.AllSomaticNeurons = false;
    params.AllResponsiveNeurons = true;
    params.fixedWindow = false;
    params.speed = 1; %min =1, max = 2;
    params.noEyeMoves = false
    params.reduceFactor = 20
    params.allCombined = false
    params.eye_to_monitor_distance = 21.5 % Distance from eye to monitor in cm
    params.pixel_size = 33
    params.resolution = 1080
    params.meanAllNeurons = false %get mean of receptive fields
    params.PaperFig logical = false
    params.OneDirection string = "all"
    params.OneLuminosity string = "all"
    params.OneSize string = "all"
    params.colorbarLims = []

end

if params.inputParams,disp(params),return,end

Stats = obj.ShufflingAnalysis;
RFs = obj.CalculateReceptiveFields('speed',params.speed);
%Parameters
%check receptive field neurons first
fieldName = sprintf('Speed%d', params.speed);
pvals = Stats.(fieldName).pvalsResponse;

responses = obj.ResponseWindow;
uDir = unique(responses.(fieldName).C(:,2));
uSize = unique(responses.(fieldName).C(:,4));
uLum = unique(responses.(fieldName).C(:,6));

%%% switch cases to plot one or all cases of one parameter
if params.OneDirection ~= "all"
    switch params.OneDirection
        case "up"
            dirIDX = find(uDir==0);
        case "left"
            dirIDX = find(uDir==1.57);
        case "down"
            dirIDX = find(uDir==3.14);
        case "right"
            dirIDX = find(uDir==1.57);
        otherwise
            error("Unknown inputPa value: %s", params.OneDirection)
    end
    DirectionSelected = uDir(dirIDX);
else
     DirectionSelected = uDir;
end

if params.OneLuminosity ~= "all"
    switch params.OneLuminosity
        case "black"
            lumIDX = find(uLum==1);
        case "white"
            lumIDX = find(uLum==255);
        otherwise
            error("Unknown inputPa value: %s", params.OneLuminosity)
    end
    LuminositySelected = uLum(lumIDX);
else
     LuminositySelected = uLum;
end

if params.OneSize ~= "all"
    switch params.OneSize
        case "small"
            sizeIDX = 1;
        case "middle"
            sizeIDX = 2;
        case "big"
            sizeIDX = 3;
        otherwise
            error("Unknown inputPa value: %s", params.OneLuminosity)
    end
    SizeSelected = uSize(sizeIDX);
else
     SizesSelected = uSize;
end


if isnan(params.exNeurons)
    if params.AllSomaticNeurons
        eNeuron = 1:numel(pvals);
        pvals = [eNeuron;pvals(eNeuron)];
    elseif params.AllResponsiveNeurons
        eNeuron = find(pvals<0.05);
        pvals = [eNeuron;pvals(eNeuron)];% Select all good neurons if not specified
        if isempty(eNeuron)
            fprintf('No responsive neurons.\n')
            return
        end
    end
else
    eNeuron = params.exNeurons;
    pvals = [eNeuron;pvals(eNeuron)];
end


coorRect = obj.VST.rect';
reduceFactor = min([params.reduceFactor min(obj.VST.ballSizes)]); %has to be bigger than the smallest ball size
redCoorX = round(coorRect(3)/reduceFactor);
redCoorY = round(coorRect(4)/reduceFactor);

pixel_size = params.pixel_size/(params.resolution/reduceFactor); % Size of one pixel in cm (e.g., 25 micrometers)
monitor_resolution = [redCoorX, redCoorY]; % Width and height in pixels
[theta_x,theta_y] = pixels2eyeDegrees(params.eye_to_monitor_distance,pixel_size,monitor_resolution);

theta_x = theta_x(:,1+(redCoorX-redCoorY)/2:(redCoorX-redCoorY)/2+redCoorY);

if params.noEyeMoves %%%mode quadrant
    RFu = squeeze(load(sprintf('NEM-RFuST-Q1-Div-X-%s',NP.recordingName)).RFuST);
else
    RFu = RFs.RFuST; %Sum of RFUs

    RFuDirSizeLum = RFs.RFuDirSizeLumFilt; %Size and dir and lum


end

offsetN = numel(unique(obj.VST.offsets));
TwoDGaussian = fspecial('gaussian',floor(size(RFu,2)/(offsetN/2)),redCoorY/offsetN); %increase size of gaussian by 100%.

for u = eNeuron


    ru = find(eNeuron == u);

    if params.allCombined

        % %%%Filter with gaussian:

        figRF=figure;
        imagesc((squeeze(conv2(RFu(:,:,ru),TwoDGaussian,'same'))));

        c = colorbar;
        title(c,'spk/s')

        colormap('turbo')
        title(sprintf('u-%d',u))

        xt = xticks;
        xt = xt((1:2:numel(xt)));
        xticks(xt);
        xticklabels(round(theta_x(1,xt)))

        yt = yticks;
        yt = yt(1:2:numel(yt));
        yticks(yt);
        yticklabels(round(theta_y(yt,1)))

        axis equal tight

        figRF.Position = [ 680   577   156   139];
        if  params.noEyeMoves
            %print(figRF, sprintf('%s-NEM-MovBall-ReceptiveField-eNeuron-%d.pdf',NP.recordingName,u), '-dpdf', '-r300', '-vector');
            if params.PaperFig
                if params.overwrite,obj.printFig(figRF,sprintf('%s-NEM-MovBall-ReceptiveField-eNeuron-%d.pdf',obj.dataObj.recordingName,u), PaperFig = params.PaperFig),end
            else
                if params.overwrite,obj.printFig(figRF,sprintf('%s-NEM-MovBall-ReceptiveField-eNeuron-%d.pdf',obj.dataObj.recordingName,u)),end
            end

        else
            if params.PaperFig
                if params.overwrite,obj.printFig(figRF,sprintf('%s-MovBall-ReceptiveField-eNeuron-%d',obj.dataObj.recordingName,u), PaperFig = params.PaperFig),end
            else
                if params.overwrite,obj.printFig(figRF,sprintf('%s-MovBall-ReceptiveField-eNeuron-%d',obj.dataObj.recordingName,u)),end
            end
        end
        
    end

    %%%% Plot receptive field per direction
    %%%% find max and min of colorbar limits


    if params.meanAllNeurons
        RFuRed =reshape(mean(RFuDirSizeLum,6),[size(RFuDirSizeLum,1),size(RFuDirSizeLum,2),...
            size(RFuDirSizeLum,3),size(RFuDirSizeLum,4)...
            ,size(RFuDirSizeLum,5)]); %%Takes mean across all neurons
        
        for i = 1:numel(hasNotString) %Take mean of elements that are not going to be compared (like luminosities, or directions, etc)
            RFuRed = mean(RFuRed,hasNotString(i));
            size(RFuRed)
        end
    else
        RFuRed =reshape(RFuDirSizeLum(:,:,:,:,:,ru),[size(RFuDirSizeLum,1),size(RFuDirSizeLum,2),size(RFuDirSizeLum,3),size(RFuDirSizeLum,4)...
            ,size(RFuDirSizeLum,5)]);

        if params.OneSize ~= "all"
            RFuRed = RFuRed(:,sizeIDX,:,:,:);
        end

        if params.OneLuminosity ~=  "all"
             RFuRed = RFuRed(:,:,lumIDX,:,:);
        end
        
        if params.OneDirection ~=  "all"
             RFuRed = RFuRed(dirIDX,:,:,:,:);
        end
    end
  
    cMax = max(RFuRed,[],'all');
    cMin = min(RFuRed,[],'all');

    tilesSize = prod(size(RFuRed,[1 2 3]));

    if numel(tilesSize) ==1 %%Create tile grid for RF ploting 
        if tilesSize<4
            tilesSize = [1 tilesSize];
        else
            tilesSize = [floor(tilesSize/2) ceil(tilesSize/2)];
        end
    end

    if numel(tilesSize) ==3 %%Create tile grid for RF ploting
        tilesSize = [tilesSize(1) tilesSize(2)*tilesSize(3)];
    end


    %%%%%%%%%%%%%%% Create tiled plot showcasing different luminosities and
    %%%%%%%%%%%%%%% directions
    figRF = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Full screen figure;
    NeuronLayout = tiledlayout(tilesSize(1),tilesSize(2),"TileSpacing","tight","Padding","compact");

    j=0;

    for d = 1:size(RFuRed,1)
        for s = 1:size(RFuRed,2)
            for l = 1:size(RFuRed,3)

                ax = nexttile;
                imagesc((squeeze(RFuRed(d,s,l,:,:))));

                caxis([cMin cMax]);

                axi = gca;
                axi.YAxis.FontSize = 8;
                axi.YAxis.FontName = 'helvetica';

                xlabel('Degrees','FontSize',10,'FontName','helvetica')
                ylabel('Degrees','FontSize',10,'FontName','helvetica')

                axi = gca;
                axi.XAxis.FontSize = 8;
                axi.XAxis.FontName = 'helvetica';

                colormap('turbo')
                title(sprintf('Dir-%s-Size-%s-Lum-%s',string(uDir(d)),string(uSize(s)),string(uLum(l))),'FontSize',4)

                %xlim([(redCoorX-redCoorY)/2 (redCoorX-redCoorY)/2+redCoorY])
                xt = xticks;%(linspace((redCoorX-redCoorY)/2,(redCoorX-redCoorY)/2+redCoorY,offsetN*2));
                xt = xt((1:2:numel(xt)));
                xticks(xt);
                xticklabels(round(theta_x(1,xt)))

                yt = yticks;
                yt = yt(1:2:numel(yt));
                yticks(yt);
                yticklabels(round(theta_y(yt,1)))

                j = j+1;
                
                if j ==size(RFuRed,1)*size(RFuRed,2)*size(RFuRed,3)
                    c = colorbar;
                    title(c,'spk/s','FontSize',8,'FontName','helvetica')

                    if ~isempty(params.colorbarLims)
                        clim([params.colorbarLims]);
                    end
                    [colorbarLims] = c.Limits;
                end
                
                axis(ax, 'equal'); 
                pbaspect(ax, [1 1 1]);

            end
        end
    end

    
    %title(NeuronLayout, sprintf('Unit-%d',u),'FontSize',4);
    %figRF.Position = [  0.2328125                     0.315                0.23515625                   0.38125];

    Sdir= strjoin(string(DirectionSelected),"-");
    Ssize= strjoin(string(SizesSelected),"-");
    Slum= strjoin(string(LuminositySelected),"-");
    
    if params.meanAllNeurons
        title(NeuronLayout,'MeanAllUnits');
        if params.overwrite,obj.printFig(figRF,sprintf('%s-%s-MovBall-RF-sep-%s-Mean',...
                obj.dataObj.recordingName,fieldName, sprintf('Dir-%s-Size-%s-Lum-%s',Sdir,Ssize,Slum))),end
        return
    end

    if tilesSize(1) == 1
        set(figRF, 'Units', 'centimeters');
        set(figRF, 'Position', [2 2 4 4]);
    end


    if params.noEyeMoves

    else
        if params.PaperFig
            if params.overwrite,obj.printFig(figRF,sprintf('%s-%s-MovBall-RF-sep-%s-eNeuron-%d',...
                    obj.dataObj.recordingName,fieldName, sprintf('Dir-%s-Size-%s-Lum-%s',Sdir,Ssize,Slum),u),PaperFig = params.PaperFig),end
        else
            if params.overwrite,obj.printFig(figRF,sprintf('%s-%s-MovBall-RF-sep-%s-eNeuron-%d',...
                    obj.dataObj.recordingName,fieldName, sprintf('Dir-%s-Size-%s-Lum-%s',Sdir,Ssize,Slum),u)),end
        end
    end

    if u ~= eNeuron(end)
        close
    end


end %%%End onDir

end
