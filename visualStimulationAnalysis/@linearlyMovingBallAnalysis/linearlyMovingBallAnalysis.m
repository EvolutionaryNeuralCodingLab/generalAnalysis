classdef linearlyMovingBallAnalysis < VStimAnalysis

    properties

    end

    properties (Constant)
        trialType = 'videoTrials'
    end

    methods (Hidden)
        %class constructor - name of class should be identical to the visual stimulation with the addition of Analysis
        function [obj] = linearlyMovingBallAnalysis(dataObj)
            if nargin==0
                dataObj=[];
            end
            % Call superclass constructor
            obj@VStimAnalysis(dataObj);
        end
    end

    methods

        plotSpatialTuningLFP(obj,params)
        plotSpatialTuningSpikes(obj,params)

        function result = getCorrSpikePattern(obj,varargin)
            %Plots the correlation matrix between the responses of all pairs of stimlui
            T = obj.getSyncedDiodeTriggers;

            %order trials accoridng ot direction of movement and offset for each direction.
            [~,pOrdered2]=sort(obj.VST.offsets);
            [originalOrder,pOrdered1]=sort(obj.VST.directions(pOrdered2));
            pOrdered=pOrdered2(pOrdered1);
            trialCat="Dir="+ num2str(obj.VST.directions(pOrdered)',2)+",offset=" + num2str(obj.VST.offsets(pOrdered)',4);

            result = getCorrSpikePattern@VStimAnalysis(obj,T.stimOnFlipTimes(pOrdered),trialCat,'win',obj.VST.stimDuration*1000,varargin{:});

        end

        % function results = setUpAnalysis(obj, params)
        function results = ResponseWindow(obj, params)
            arguments (Input)
                obj
                params.oneSpeed logical = true
                params.overwrite logical = false
                params.analysisTime = datetime('now')
                params.inputParams = false
                params.binRaster = 1
                params.durationWindow = 100
                params.preBase = 200
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
            try
                DiodeCrossings = obj.getSyncedDiodeTriggers;
            catch
                obj.getDiodeTriggers("extractionMethod",'digitalTriggerDiode','overwrite',true);
                DiodeCrossings = obj.getSyncedDiodeTriggers;
            end


            stimOn = DiodeCrossings.stimOnFlipTimes;
            stimOff = DiodeCrossings.stimOffFlipTimes;

            A = [stimOn' obj.VST.directions' obj.VST.offsets' obj.VST.ballSizes' obj.VST.speeds'];
            [C indexS] = sortrows(A,[2 3 4 5]);

            B = [stimOff' obj.VST.directions' obj.VST.offsets' obj.VST.ballSizes' obj.VST.speeds'];
            [Coff indexSo] = sortrows(B,[2 3 4 5]);

            stimInter= mean(stimOn(2:end)-stimOff(1:end-1));

            if params.oneSpeed

                %Select maximum speed
                A =  A(A(:,5) == max(A(:,5)),:);

                B =  B(B(:,5) == max(B(:,5)),:);

                [C indexS] = sortrows(A,[2 3 4 5]);

                [Coff indexSo] = sortrows(B,[2 3 4 5]);

                stimOn = A(:,1);
                stimOff = B(:,1);
                stimInter = obj.VST.interTrialDelay*1000;
                stimDur = mean(-stimOn+stimOff);

            end

            preBase = round(stimInter-params.preBase);
            % Load Kilosort and phy results
            p = obj.dataObj.convertPhySorting2tIc(obj.dataObj.recordingDir);
            label = string(p.label');
            goodU = p.ic(:,label == 'good');

            %4. Sort directions:
            directimesSorted = C(:,1)';
            Mr = BuildBurstMatrix(goodU,round(p.t/params.binRaster),round(directimesSorted/params.binRaster),round((stimDur+ params.durationWindow)/params.binRaster)); %response matrix
            [MrC]=ConvBurstMatrix(Mr,fspecial('gaussian',[1 5],3),'same');

            [nT,nN,nB] = size(Mr);

            window_size = [1, round(params.durationWindow/params.binRaster)];
            trialDivision = nT/(numel(unique(obj.VST.directions))*numel(unique(obj.VST.offsets))*numel(unique(obj.VST.ballSizes))*numel(unique(obj.VST.speeds)));

            %5. Initialize the maximum mean value and its position
            max_position = zeros(nN,2);
            max_mean_value = zeros(1,nN);
            max_mean_valueB = zeros(1,nN);
            NeuronVals = zeros(nN,nT/trialDivision,8);

            %Baseline = size window

            [Mbd] = BuildBurstMatrix(goodU,round(p.t/params.binRaster),round((directimesSorted-preBase)/params.binRaster),round(preBase/params.binRaster)); %Baseline matrix plus

            %Merge trials:

            mergeTrials = trialDivision;

            Bd = reshape(Mbd, [mergeTrials, size(Mbd, 1)/mergeTrials, size(Mbd, 2), size(Mbd,3)]);

            Mbd2 = squeeze(mean(Bd, 1));

            B = reshape(MrC, [mergeTrials, size(Mr, 1)/mergeTrials, size(Mr, 2), size(Mr,3)]);

            % Take the mean across the first dimension (rows)
            Mr2 = squeeze(mean(B, 1));
            [nT,nN,nB] = size(Mr2);
            %Real data:
            for u =1:nN
                % Slide the window over the matrix
                %unit matrix
                max_mean_value(u) = -Inf; %General max? not needed
                max_mean_valueB(u)=-Inf;
                NeuronRespProfile = zeros(nT,8); %4 columns plus: ofsett, dir, size, speed, frec.

                k =1;
                max_position_Trial = zeros(nT,2); %Create 2 matrices, for mean inside max window, and for position of window, for each trial category
                max_mean_value_Trial = zeros(1,nT);
                max_mean_value_TrialB = zeros(1,nT);
                for i = 1:nT %Iterate across trials
                    uM = squeeze(Mr2(i,u,:))';%Create matrix per unique trial conditions

                    uMB = squeeze(Mbd2(i,u,:))';%Create matrix per unique trial conditions
                    %uMb = Mbd2(i,u);

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
                    %Assign visual stats to last columns of NeuronRespProfile. Select
                    %according  to trial (d) the appropiate parameter > directions'offsets' sizes' speeds' freq'
                    NeuronRespProfile(k,5) = C(i*mergeTrials,2);
                    NeuronRespProfile(k,6) = C(i*mergeTrials,3);
                    NeuronRespProfile(k,7) = C(i*mergeTrials,4);
                    NeuronRespProfile(k,8) = C(i*mergeTrials,5);
                    k = k+1;

                end

                %figure;imagesc(uM);xline(max_position_Trial(i,2));xline(max_position_Trial(i,2)+window_size(2))
                NeuronVals(u,:,:) = NeuronRespProfile;
            end

            colNames = {'Resp','MaxWinTrial','MaxWinBin','RespSubsBaseline','Directions','Offsets','Sizes','Speeds'};

            %save results in the right file
            fprintf('Saving results to file.\n');
            save(obj.getAnalysisFileName,'params','NeuronVals','C','Coff','goodU','stimDur','stimInter','colNames');
        end

        function results = ShufflingAnalysis(obj, params)

            arguments (Input) %ResponseWindow.mat
                obj
                params.overwrite logical = false
                params.analysisTime = datetime('now')
                params.inputParams = false
                params.trialThreshold = 0.6
                params.N_bootstrap = 2000;
                % params.preBase = 200;
                % params.durationWindow = 100
                % params.bin = 1;
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
            p = obj.dataObj.convertPhySorting2tIc(obj.dataObj.recordingDir);

            trialDivision = numel(directimesSorted)/numel(unique(NeuronResp.C(:,2)))/numel(unique(NeuronResp.C(:,3)))/numel(unique(NeuronResp.C(:,4)))...
                /numel(unique(NeuronResp.C(:,5)));

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
            
            %save results in the right file
            fprintf('Saving results to file.\n');
            save(obj.getAnalysisFileName,'params','pvalsResponse','ZScoreU','boot_means');
        end
    
    
    end %end Methods
end %end clas

        %%%1. Get Diode
        %%%2. Create A Matrix.
        %%%3. Load Kilos0rt and phy results.
        %%%4. Create response matrix.
%%%5. Create shuffling analysis
%%%6. Create receptive field analysis with eyes moving.
%%%7. Create receptive field analysis with eyes not moving. 