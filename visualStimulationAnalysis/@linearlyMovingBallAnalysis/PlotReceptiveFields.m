function PlotReceptiveFields(obj,params)

arguments (Input)
    obj
    params.overwrite logical = false
    params.analysisTime = datetime('now')
    params.inputParams = false
    params.exNeurons = 1;
    params.AllSomaticNeurons = false;
    params.AllResponsiveNeurons = true;
    params.fixedWindow = false;
    params.speed = 1; %min =1, max = 2;
    params.noEyeMoves = false
    params.reduceFactor = 20
    params.allDir = false
    params.eye_to_monitor_distance = 21.5 % Distance from eye to monitor in cm
    params.pixel_size = 33
    params.resolution = 1080
end

if params.inputParams,disp(params),return,end

Stats = obj.ShufflingAnalysis;
RFs = obj.CalculateReceptiveFields('speed',params.speed);
%Parameters
%check receptive field neurons first
fieldName = sprintf('Speed%d', params.speed);
pvals = Stats.(fieldName).pvalsResponse;


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
    RFu = RFs.RFuST;
    q=1;
    rf{1} = RFs.RFuSTDirSizeFilt; %Size and dir
    RF{1} = rf;
    names{1} ="";
end

offsetN = numel(unique(obj.VST.offsets));
TwoDGaussian = fspecial('gaussian',floor(size(RFu,2)/(offsetN/2)),redCoorY/offsetN); %increase size of gaussian by 100%.

for u = eNeuron

    ru = find(eNeuron == u);

    if params.allDir

        % %%%Filter with gaussian:

        figRF=figure;
        imagesc((squeeze(conv2(RFu(:,:,ru),TwoDGaussian,'same'))));

        % figRF=figure;
        % imagesc((squeeze(conv2(RFu(:,:),TwoDGaussian,'same'))));

        caxis

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
        %                     xlabel('X degrees')
        %                     ylabel('Y degrees')

        axis equal tight

        figRF.Position = [ 680   577   156   139];
        if  params.noEyeMoves
            %print(figRF, sprintf('%s-NEM-MovBall-ReceptiveField-eNeuron-%d.pdf',NP.recordingName,u), '-dpdf', '-r300', '-vector');
            if params.overwrite,obj.printFig(figRF,sprintf('%s-NEM-MovBall-ReceptiveField-eNeuron-%d.pdf',obj.dataObj.recordingName,u)),end
        else
            %print(figRF, sprintf('%s-MovBall-ReceptiveField-eNeuron-%d.pdf',NP.recordingName,u), '-dpdf', '-r300', '-vector');
            if params.overwrite,obj.printFig(figRF,sprintf('%s-MovBall-ReceptiveField-eNeuron-%d',obj.dataObj.recordingName,u)),end
        end
        close

    else %%%% Plot receptive field per direction

        if   params.noEyeMoves

            if string(DivisionType) == "XY"

                file = dir (NP.recordingDir);
                filenames = {file.name};
                cd(NP.recordingDir)

                names = {'X','Y'};

                files = filenames(contains(filenames,"NEM-RFuSTDirSizeFilt")& contains(filenames,"Div"));


                for n = 1:2

                    filesD = filenames(contains(files,names{n}));

                    for f = 1:numel(filesD)
                        rf{1,f} = load(sprintf('NEM-RFuSTDirSizeFilt-Q%d-Div-%s-%s',f,names{n},NP.recordingName)).RFuSTDirSizeFilt;

                        rf{2,f} = sprintf('NEM-RFuSTDirSizeFilt-Q%d-Div-%s-%s',f,names{n},NP.recordingName);
                    end

                    RF{1,n} = rf;

                    RF{2,n} = sprintf('Division along %s',names{n});
                end

                %RF contains two cells, one for each type of eye movement
                %division (along X or along Y)

                %                 q=1;
                %                 RF{1} = load(sprintf('NEM-RFuSTDirSizeFilt-Q%d-%s',q,NP.recordingName)).RFuSTDirSizeFilt; %Size and dir
                %                 q=2;
                %                 RF{2} = load(sprintf('NEM-RFuSTDirSizeFilt-Q%d-%s',q,NP.recordingName)).RFuSTDirSizeFilt; %Size and dir

            else

                rf{1} = load(sprintf('NEM-RFuSTDirSizeFilt-Q1-Div-X-%s',NP.recordingName)).RFuSTDirSizeFilt; %Size and dir

                RF{1} = rf;

                names{1} ="";
            end
        else
            q=1;
            rf{1} = load(sprintf('RFuSTDirSizeFilt-%s',NP.recordingName)).RFuSTDirSizeFilt; %Size and dir

            RF{1} = rf;

            names{1} ="";
        end

        %%%% find max and min of colorbar limits

        cMax = -inf;
        cMin = inf;
        for n = 1:size(RF,2)

            rf = RF{1,n};
            for q = 1:size(rf,2)
                RFuSTDirSizeFilt = rf{1,q};

                if size(RFuSTDirSizeFilt,2) >1 %If there is more than on size %%Solve for 2 quadrants
                    %%Select size closest to 120
                    [minS indx] = min(abs(unique(sizes)-120));
                    RFuSTDirFilt = squeeze(RFuSTDirSizeFilt(:,indx,:,:,:));
                else
                    RFuSTDirFilt = squeeze(RFuSTDirSizeFilt);
                end

                for d = 1:direcN

                    if cMax < max(RFuSTDirFilt(:,:,:,ru),[],'all')

                        cMax = max(RFuSTDirFilt(:,:,:,ru),[],'all');

                    end

                    if cMin > min(RFuSTDirFilt(:,:,:,ru),[],'all')

                        cMin = min(RFuSTDirFilt(:,:,:,ru),[],'all');

                    end

                end
            end

        end

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figRF = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Full screen figure;
        NeuronLayout = tiledlayout(size(RF,2),size(RF,2),"TileSpacing","tight","Padding","tight");

        tl = 1;
        for n = 1:size(RF,2)

            rf = RF{1,n};

            for q = 1:size(rf,2)

                if names{n} == 'X'
                    qNames = {'Right','Left'};

                elseif names{n} == 'Y'
                    qNames = {'Down','Up'};
                else
                    qNames = {"-"};
                end

                RFuSTDirSizeFilt = rf{1,q};

                if size(RFuSTDirSizeFilt,2) >1 %If there is more than on size %%Solve for 2 quadrants
                    %%Select size closest to 120
                    [minS indx] = min(abs(unique(sizes)-120));
                    RFuSTDirFilt = squeeze(RFuSTDirSizeFilt(:,indx,:,:,:));
                else
                    RFuSTDirFilt = squeeze(RFuSTDirSizeFilt);
                end


                xqLayout=tiledlayout(NeuronLayout,direcN/2,direcN/2,"TileSpacing","compact");
                xqLayout.Layout.Tile = tl;

                %DirLayout=tiledlayout(direcN/2,direcN/2,"TileSpacing","tight","Padding","none");


                for d = 1:direcN

                    nexttile(xqLayout);

                    if d ==1 || d==3
                        %imagesc(rot90(squeeze(RFuSTDirFilt(d,:,1+(redCoorX-redCoorY)/2:(redCoorX-redCoorY)/2+redCoorY,ru)),2));

                        imagesc((squeeze(RFuSTDirFilt(d,:,:,ru))));

                        %                         c = colorbar;
                        2+2
                    else
                        imagesc((squeeze(RFuSTDirFilt(d,:,:,ru))));

                    end

                    if d==2
                        c = colorbar;
                        title(c,'Spike Rate (Hz)')
                    end



                    caxis([cMin cMax]);

                    colormap('turbo')
                    title(string(uDir(d)))

                    %xlim([(redCoorX-redCoorY)/2 (redCoorX-redCoorY)/2+redCoorY])
                    xt = xticks;%(linspace((redCoorX-redCoorY)/2,(redCoorX-redCoorY)/2+redCoorY,offsetN*2));
                    xt = xt((1:2:numel(xt)));
                    xticks(xt);
                    xticklabels(round(theta_x(1,xt)))

                    yt = yticks;
                    yt = yt(1:2:numel(yt));
                    yticks(yt);
                    yticklabels(round(theta_y(yt,1)))
                    %                     xlabel('X degrees')
                    %                     ylabel('Y degrees')
                    if d == 1 || d==2

                        xticks([])

                    end

                    if d==2 || d==4

                        yticks([])

                    end

                    %%%%Use same colorbar scale

                    axis equal tight


                end

                ti = title(xqLayout,sprintf('nUnit-%d-Div-%s-Q-%s\n\',eNeuron(ru),names{n},qNames{q}),'FontSize',10);

                %ti.Position(2) = ti.Position(2) + 0.05;

                tl = tl+1;

                %         %%%%%%%%% plot rasters used in noEyeMoves RFs
                %
                %         spikeSum = spikeSums{q};
                %
                %         eNspkS = squeeze(spikeSum(:,ru,:));
                %         figure;imagesc(eNspkS);colormap(flipud(gray(64)));
                %         %     rowsWithNaN = find(any(isnan(eNspkS), 2));
                %         %     yline(rowsWithNaN,'g')
                %         yline(trialDivision*sizeN:trialDivision*sizeN:size(spikeSums{q},1));
                %         yline(trialDivision*offsetN*sizeN:trialDivision*offsetN*sizeN:size(spikeSums{q},1)-1,'r','LineWidth',5)
                %         title(string(q))

            end %%%%%%%%%%%%%%%%
        end

        if savePlot
            cd(saveDir)
            figRF.Position = [  0.2328125                     0.315                0.23515625                   0.38125];
            if noEyeMoves
                print(gcf, sprintf('%s-NEM-MovBall-ReceptiveField-Div-%s-Q-%d-eNeuron-%d.pdf',NP.recordingName,names{n},q,u), '-dpdf', '-r300', '-vector');
            else
                print(gcf, sprintf('%s-MovBall-ReceptiveField-Div-%s-Q-%d-eNeuron-%d.pdf',NP.recordingName,names{n},q,u), '-dpdf', '-r300', '-vector');
            end
        end

        if params.overwrite,obj.printFig(f1,'On_SpkVsPosition'),end

    end %%%End onDir

end

end