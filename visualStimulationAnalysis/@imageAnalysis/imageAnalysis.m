classdef imageAnalysis < VStimAnalysis

    properties
        Session
    end

    properties (Constant)
        trialType = 'imageTrials'
    end

    methods (Hidden)
        %class constructor - name of class should be identical to the visual stimulation with the addition of Analysis
        function [obj] = imageAnalysis(dataObj,params)
            arguments (Input) %ResponseWindow.mat
                dataObj
                params.Session = 1;
            end
            if nargin==0
                dataObj=[];
            end
            % Call superclass constructor
            obj@VStimAnalysis(dataObj,'Session',params.Session);
            obj.Session = params.Session;
        end
    end

    %static methods that are not directly related to the class are defined here
    methods (Static)
        activityTracePhysicalSpacePlot
    end

    methods

        plotSpatialTuningLFP(obj,params)

        function result = getCorrSpikePattern(obj,varargin)
            %Plots the correlation matrix between the responses of all pairs of stimlui
            T = obj.getSyncedDiodeTriggers;

            %order trials accoridng ot direction of movement and offset for each direction.
            [~,pOrdered]=sort(obj.VST.pos);

            trialCat="X="+ num2str(obj.VST.pos2X(obj.VST.pos(pOrdered)'),2)+",Y=" + num2str(obj.VST.pos2Y(obj.VST.pos(pOrdered)'),2);
            %Adds to the window 600 ms of the off response
            result = getCorrSpikePattern@VStimAnalysis(obj,T.stimOnFlipTimes(pOrdered),trialCat,'win',obj.VST.stimDuration*1000+600,varargin{:});

        end

        function results = ResponseWindow(obj,params)
            arguments (Input)
                obj
                params.overwrite logical = false
                params.analysisTime = datetime('now')
                params.inputParams = false
                params.binRaster = 1
                params.durationWindow = 100
                params.preBase = 200
                params.GaussianLength =5;
            end
            if params.inputParams,disp(params),return,end

            if isfile(obj.getAnalysisFileName) && ~params.overwrite
                if nargout==1
                    fprintf('Loading saved results from file.\n');
                    results=load(obj.getAnalysisFileName);
                else
                    fprintf('Analysis already exists (use overwrite option to recalculate).\n');
                end

                return
            end

            ST = obj.getSessionTime;

            DiodeCrossings = obj.getSyncedDiodeTriggers;
            
            stimOn = DiodeCrossings.stimOnFlipTimes;
            stimOff = DiodeCrossings.stimOffFlipTimes;
            stimInter = obj.VST.interTrialDelay*1000;

            stimDur = mean(-stimOn+stimOff);

            imgOrder = obj.VST.imgSequence(obj.VST.order)';

            A = [stimOn' imgOrder];
            [C indexS] = sortrows(A,[2]);

            B = [stimOff' imgOrder];
            [Coff indexSo] = sortrows(B,[2]);

            directimesSorted = C(:,1)';
            bin = params.binRaster;

            trialDivision = numel(directimesSorted)/numel(obj.VST.imgNames);

            % Load Kilosort and phy results
            p = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
            label = string(p.label');
            goodU = p.ic(:,label == 'good');

            %%Create window to scan rasters and get the maximum response
            duration = params.durationWindow; %window in ms, same as in MB
            Mr = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted)/bin),round((stimDur+duration)/bin)); %response matrix
            [MrC]=ConvBurstMatrix(Mr,fspecial('gaussian',[1 params.GaussianLength],3),'same');

            [nT,nN,nB] = size(MrC);

            window_size = [1, round(params.durationWindow /bin)];

            %5. Initialize the maximum mean value and its position

            max_position = zeros(nN,2);
            max_mean_value = zeros(1,nN);
            max_mean_valueB = zeros(1,nN);


            NeuronVals = zeros(nN,nT/(trialDivision),5); %Each neuron, which has a matrix where the first column is maxVal of bins, 2nd and 3rd position of window in matrix...
            % 4th Z-score.
            % responsive compared to the baseline.

            %[1, 2, 3, 7, 8, 9, 10, 11, 12, 13, 17, 20, 22, 23, 24, 28, 32, 34, 36]

            mergeTrials = round(trialDivision);

            %Baseline = size window
            preBase = params.preBase;

            [Mbd] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round(preBase/bin)); %Baseline matrix plus

            %Merge trials:

            Bd = reshape(Mbd, [mergeTrials, size(Mbd, 1)/mergeTrials, size(Mbd, 2), size(Mbd,3)]);

            Mbd2 = squeeze(mean(Bd, 1));

            B = reshape(MrC, [mergeTrials, size(Mr, 1)/mergeTrials, size(Mr, 2), size(Mr,3)]);

            % Take the mean across the first dimension (rows)
            Mr2 = squeeze(mean(B, 1));

            [nT,nN,nB] = size(Mr2);

            for u =1:nN
                % Slide the window over the matrix
                %unit matrix
                max_mean_value(u) = -Inf; %General max? not needed
                max_mean_valueB(u)=-Inf;
                NeuronRespProfile = zeros(nT,5); %4 columns plus: ofsett, dir, size, speed, frec.

                k =1;
                max_position_Trial = zeros(nT,2); %Create 2 matrices, for mean inside max window, and for position of window, for each trial category
                max_mean_value_Trial = zeros(1,nT);
                max_mean_value_TrialB = zeros(1,nT);

                for i = 1:nT %Iterate across trials
                    uM = squeeze(Mr2(i,u,:))';%Create matrix per unique trial conditions
                    uMB = squeeze(Mbd2(i,u,:))';%Create matrix per unique trial conditions

                    max_mean_value_Trial(k) = -Inf;
                    max_mean_value_TrialB(k) = -Inf;

                    for j = 1:2:size(uM, 2) - window_size(2) + 1 %Iterate across bins
                        % Extract the sub-matrix
                        sub_matrix = uM(j:min(j+window_size(2)-1,end)); %Create matrix of size window per bin
                        sub_matrixB = uMB(j:min(j+window_size(2)-1,end));

                        % Compute the mean value
                        mean_value = mean(sub_matrix(:)); %Compute mean of each window
                        mean_valueB = mean(sub_matrixB(:)); %Compute mean of each window

                        % Update the maximum mean value and its position (a
                        % window is selected across each trial division
                        if mean_value >  max_mean_value_Trial(k)
                            max_mean_value_Trial(k) = mean_value;
                            max_position_Trial(k,:) = [i j];
                        end

                        if mean_valueB >  max_mean_value_TrialB(k)
                            max_mean_value_TrialB(k) = mean_valueB;
                        end

                    end
                    %Save across each trial (in a row) the max bin window
                    %1%. Response
                    NeuronRespProfile(k,1) = max_mean_value_Trial(k);

                    %2%. WindowTrial
                    NeuronRespProfile(k,2) = max_position_Trial(k,1);

                    %3%. WindowBin
                    NeuronRespProfile(k,3) = max_position_Trial(k,2);

                    %4%. Resp - Baseline
                    % NeuronRespProfile(i,4) = (max_mean_value_Trial(i) - (spkRateBM(u)+max_mean_value_Trial(i))/2)/denom(u); %Zscore
                    %NeuronRespProfile(k,4) = (max_mean_value_Trial(k) - spkRateBM(u))/denom(u); %Zscore
                    NeuronRespProfile(k,4) = max_mean_value_Trial(k)-max_mean_value_TrialB(k);
                    k = k+1;

                end

                %figure;imagesc(uM);xline(max_position_Trial(i,2));xline(max_position_Trial(i,2)+window_size(2))

                NeuronVals(u,:,:) = NeuronRespProfile;
            end            

            colNames = {'Resp','MaxWinTrial','MaxWinBin','RespSubsBaseline','Position','Size','Luminosities'};

            S.params = params;
            S.colNames = {colNames};
            S.goodU = goodU;
            S.stimInter = stimInter;
            S.stimDur = stimDur;
            S.C = C;
            S.NeuronVals = NeuronVals;

            %save results in the right file
            fprintf('Saving results to file.\n');
            save(obj.getAnalysisFileName,'-struct', 'S');
            results = S;

        end

        function results = ShufflingAnalysis(obj,params)
            arguments (Input) %ResponseWindow.mat
                obj
                params.overwrite logical = false
                params.analysisTime = datetime('now')
                params.inputParams = false
                params.trialThreshold = 0.6
                params.N_bootstrap = 2000;
                
            end
            if params.inputParams,disp(params),return,end

            if isfile(obj.getAnalysisFileName) && ~params.overwrite
                if nargout==1
                    fprintf('Loading saved results from file.\n');
                    results=load(obj.getAnalysisFileName);
                else
                    fprintf('Analysis already exists (use overwrite option to recalculate).\n');
                end

                return
            end

            ResponseWindowPath = [fileparts(obj.getAnalysisFileName) filesep 'ResponseWindow.mat'];
            
            if isfile(ResponseWindowPath)
                NeuronResp = load(ResponseWindowPath);
            else
                fprintf('Run ResponseWindow method first!\n');
                return
            end

            directimesSorted = NeuronResp.C(:,1)';
            goodU = NeuronResp.goodU;
            p = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);

            trialDivision = numel(directimesSorted)/numel(obj.VST.imgNames);

            %Load same parameters used in response window calculation
            preBase = NeuronResp.params.preBase;
            duration =  NeuronResp.params.durationWindow;
            bin =  NeuronResp.params.binRaster;

            baseline = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round((preBase)/bin));
            baseline = single(baseline);
            [nT,nN,nB] = size(baseline);

            % Bootstrapping settings
            boot_means = zeros(params.N_bootstrap, nN,'single');
            resampled_indicesTr = single(randi(nT, [trialDivision, params.N_bootstrap]));% To store bootstrapped means
            resampled_indicesTi = single(randi(nB, [duration, params.N_bootstrap]));
            kernel = ones(trialDivision, duration) / (trialDivision * duration); % Normalize for mean


            tic
            for i = 1:params.N_bootstrap
                % Resample trials with replacement
                resampled_trials = baseline(resampled_indicesTr(:, i), :,resampled_indicesTi(:, i));
                for ui = 1:nN
                    % Extract the slice for the current unit (t x b matrix)
                    slice = resampled_trials(:, ui, :);
                    slice = squeeze(slice); % Result is t x b

                    % Compute the mean using 2D convolution
                    means = conv2(slice, kernel, 'valid'); % 'valid' ensures the window fits within bounds

                    % Find the maximum mean in this slice
                    boot_means(i, ui) = max(means(:));
                end
            end
            toc

            %[bootstats] = get_bootstrapped_equalsamples(data,nruns,num_trials,param)

            [respVal,respVali]= max(NeuronResp.NeuronVals(:,:,1),[],2);

            Mr = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted)/bin),round((NeuronResp.stimDur+duration)/bin)); %response matrix

            %%% Calculate p-value & Filter out neurons in which max response window is empty for more than
            %%% 60% of trials

            pvalsResponse = zeros(1,nN);
            ZScoreU = zeros(1,nN);

            for u = 1:nN

                posTr = NeuronResp.NeuronVals(u,respVali(u),2);
                posBin = NeuronResp.NeuronVals(u,respVali(u),3);

                maxWindow = squeeze(Mr(posTr*trialDivision-trialDivision+1:posTr*trialDivision,u,posBin:posBin+duration-1));

                emptyRows = sum(all(maxWindow == 0, 2));

                pvalsResponse(u) = mean(boot_means(:,u)>respVal(u));
                ZScoreU(u) = (respVal(u)-mean(boot_means(:,u)))/(std(boot_means(:,u))+1/(params.N_bootstrap*trialDivision));


                if emptyRows/trialDivision >= params.trialThreshold %%%Check which unit is 34 i respU
                    pvalsResponse(u) = 1;
                end

            end
            S.pvalsResponse = pvalsResponse;
            S.ZScoreU = ZScoreU;
            S.boot_means = boot_means;

            S.params = params;

            %save results in the right file
            fprintf('Saving results to file.\n');
            save(obj.getAnalysisFileName,'-struct','S');

            results = S;

        end
    end
end