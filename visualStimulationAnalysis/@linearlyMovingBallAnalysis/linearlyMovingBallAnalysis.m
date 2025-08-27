classdef linearlyMovingBallAnalysis < VStimAnalysis

    properties
        Session
    end

    properties (Constant)
        trialType = 'videoTrials'
    end

    methods (Hidden)
        %class constructor - name of class should be identical to the visual stimulation with the addition of Analysis
        function [obj] = linearlyMovingBallAnalysis(dataObj,params)
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

    methods

        plotSpatialTuningLFP(obj,params)
        plotSpatialTuningSpikes(obj,params)
        result = getCorrSpikePattern(obj,varargin)

        % function results = setUpAnalysis(obj, params)
        function results = ResponseWindow(obj, params)
            arguments (Input)
                obj
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

            A = [stimOn' obj.VST.directions' obj.VST.offsets' obj.VST.ballSizes' obj.VST.speeds' obj.VST.Luminosities'];
            [C indexS] = sortrows(A,[2 3 4 5]);

            B = [stimOff' obj.VST.directions' obj.VST.offsets' obj.VST.ballSizes' obj.VST.speeds' obj.VST.Luminosities'];
            [Coff indexSo] = sortrows(B,[2 3 4 5]);

            stimInter = obj.VST.interTrialDelay*1000;

            speeds = sort(unique(A(:,5)));

            x = length(speeds);

            % Generate field names for each unique speed
            topFields = arrayfun(@(n) sprintf('Speed%d', n), 1:x, 'UniformOutput', false);

            % Initialize each as an empty struct with subfields
            subFields = {'NeuronVals','C','Coff','stimDur'};  % subfield names
            emptySubStruct = cell2struct(cell(1, numel(subFields)), subFields, 2);
            S = cell2struct(repmat({emptySubStruct}, 1, x), topFields, 2);

            for s = 1:x %iterate among unique speeds

                fieldName = sprintf('Speed%d', s);

                %Select maximum speed
                At =  A(A(:,5) ==speeds(s),:);

                Bt =  B(B(:,5) ==speeds(s),:);

                [C indexS] = sortrows(At,[2 3 4 5 6]);

                [Coff indexSo] = sortrows(Bt,[2 3 4 5 6]);

                stimOnT = At(:,1);
                stimOffT = Bt(:,1);
               
                stimDur = mean(-stimOnT+stimOffT);

                preBase = round(stimInter-params.preBase);
                % Load Kilosort and phy results
                p = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
                label = string(p.label');
                goodU = p.ic(:,label == 'good');

                %4. Sort directions:
                directimesSorted = C(:,1)';
                Mr = BuildBurstMatrix(goodU,round(p.t/params.binRaster),round(directimesSorted/params.binRaster),round((stimDur+ params.durationWindow)/params.binRaster)); %response matrix
                [MrC]=ConvBurstMatrix(Mr,fspecial('gaussian',[1 5],3),'same');

                [nT,nN,nB] = size(Mr);

                window_size = [1, round(params.durationWindow/params.binRaster)];
                trialDivision = numel(directimesSorted)/numel(unique(C(:,2)))/numel(unique(C(:,3)))/numel(unique(C(:,4)))...
                    /numel(unique(C(:,5)))/numel(unique(C(:,6)));

                %5. Initialize the maximum mean value and its position
                max_position = zeros(nN,2);
                max_mean_value = zeros(1,nN);
                max_mean_valueB = zeros(1,nN);
                NeuronVals = zeros(nN,nT/trialDivision,9);

                %Baseline = size window

                [Mbd] = BuildBurstMatrix(goodU,round(p.t/params.binRaster),round((directimesSorted-preBase)/params.binRaster),round(preBase/params.binRaster)); %Baseline matrix plus

                %Merge trials:

                mergeTrials = trialDivision;

                Bd = reshape(Mbd, [mergeTrials, size(Mbd, 1)/mergeTrials, size(Mbd, 2), size(Mbd,3)]);

                Mbd2 = squeeze(mean(Bd, 1));

                Ba = reshape(MrC, [mergeTrials, size(Mr, 1)/mergeTrials, size(Mr, 2), size(Mr,3)]);

                % Take the mean across the first dimension (rows)
                Mr2 = squeeze(mean(Ba, 1));
                [nT,nN,nB] = size(Mr2);
                %Real data:
                for u =1:nN
                    % Slide the window over the matrix
                    %unit matrix
                    max_mean_value(u) = -Inf; %General max? not needed
                    max_mean_valueB(u)=-Inf;
                    NeuronRespProfile = zeros(nT,9); %4 columns plus: ofsett, dir, size, speed, lum.

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
                        NeuronRespProfile(k,9) = C(i*mergeTrials,6);
                        k = k+1;

                    end

                    %figure;imagesc(uM);xline(max_position_Trial(i,2));xline(max_position_Trial(i,2)+window_size(2))
                    NeuronVals(u,:,:) = NeuronRespProfile;
                end

                S.(fieldName).C = C;
                S.(fieldName).Coff = Coff;
                S.(fieldName).NeuronVals = NeuronVals;
                S.(fieldName).stimDur = stimDur;

            end


            colNames = {'Resp','MaxWinTrial','MaxWinBin','RespSubsBaseline','Directions','Offsets','Sizes','Speeds','Luminosities'};

            S.params = params;
            S.colNames = {colNames};
            S.goodU = goodU;
            S.stimInter = stimInter;

            %save results in the right file
            fprintf('Saving results to file.\n');
            save(obj.getAnalysisFileName,'-struct', 'S');
            results = S;
        end

        function results = ShufflingAnalysis(obj, params)

            arguments (Input) %ResponseWindow.mat
                obj
                params.overwrite logical = false
                params.analysisTime = datetime('now')
                params.inputParams = false
                params.trialThreshold = 0.6
                params.N_bootstrap = 2000;
                params.speed = 0; %min =1, middle = 2, max=3;
                params.AllSpeeds = true;
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

            NeuronResp = obj.ResponseWindow;

            if isfile(ResponseWindowPath)
                NeuronResp = load(ResponseWindowPath);
            else
                fprintf('Run ResponseWindow method first!\n');
                return
            end

            if params.AllSpeeds
                speeds = sort(unique(obj.VST.speeds));
            else
                speeds = sort(unique(obj.VST.speeds));
                speeds = speeds(params.speed);
            end

            x = length(speeds);

            % Generate field names for each unique speed
            topFields = arrayfun(@(n) sprintf('Speed%d', n), 1:x, 'UniformOutput', false);

            % Initialize each as an empty struct with subfields
            subFields = {'pvalsResponse','ZScoreU','boot_means'};  % subfield names
            emptySubStruct = cell2struct(cell(1, numel(subFields)), subFields, 2);
            S = cell2struct(repmat({emptySubStruct}, 1, x), topFields, 2);


            for s=1:x %iterate among unique speeds

                fieldName = sprintf('Speed%d', s);

                directimesSorted = NeuronResp.(fieldName).C(:,1)';
                goodU = NeuronResp.goodU;
                p = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);

                trialDivision = numel(directimesSorted)/numel(unique(NeuronResp.(fieldName).C(:,2)))/numel(unique(NeuronResp.(fieldName).C(:,3)))/...
                    numel(unique(NeuronResp.(fieldName).C(:,4)))...
                    /numel(unique(NeuronResp.(fieldName).C(:,5)))/numel(unique(NeuronResp.(fieldName).C(:,6)));

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

                [respVal,respVali]= max(NeuronResp.(fieldName).NeuronVals(:,:,1),[],2);

                Mr = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted)/bin),round((NeuronResp.(fieldName).stimDur+duration)/bin)); %response matrix

                %%% Calculate p-value & Filter out neurons in which max response window is empty for more than
                %%% 60% of trials

                pvalsResponse = zeros(1,nN);
                ZScoreU = zeros(1,nN);

                for u = 1:nN

                    posTr = NeuronResp.(fieldName).NeuronVals(u,respVali(u),2);
                    posBin = NeuronResp.(fieldName).NeuronVals(u,respVali(u),3);

                    maxWindow = squeeze(Mr(posTr*trialDivision-trialDivision+1:posTr*trialDivision,u,posBin:posBin+duration-1));

                    emptyRows = sum(all(maxWindow == 0, 2));

                    pvalsResponse(u) = mean(boot_means(:,u)>respVal(u));
                    ZScoreU(u) = (respVal(u)-mean(boot_means(:,u)))/(std(boot_means(:,u))+1/(params.N_bootstrap*trialDivision));


                    if emptyRows/trialDivision >= params.trialThreshold %%%Check which unit is 34 i respU
                        pvalsResponse(u) = 1;
                    end

                end
                S.(fieldName).pvalsResponse = pvalsResponse;
                S.(fieldName).ZScoreU = ZScoreU;
                S.(fieldName).boot_means = boot_means;

            end
            params.speed = speeds;
            S.params = params;

            %save results in the right file
            fprintf('Saving results to file.\n');
            save(obj.getAnalysisFileName,'-struct','S');
            
            results = S;
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