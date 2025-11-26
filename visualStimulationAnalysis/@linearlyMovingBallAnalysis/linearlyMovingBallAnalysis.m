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
                params.GaussianLength =5
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
                obj.getSessionTime("overwrite",true);
                obj.getDiodeTriggers("extractionMethod",'digitalTriggerDiode','overwrite',true);
                DiodeCrossings = obj.getSyncedDiodeTriggers;
            end

            stimOn = DiodeCrossings.stimOnFlipTimes;
            stimOff = DiodeCrossings.stimOffFlipTimes;

            if  ~isfield(obj.VST,'Luminosities')
                obj.VST.Luminosities = zeros(1,numel(obj.VST.directions))+255;
                
            end

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

            preBase = round(stimInter-params.preBase);
            % Load Kilosort and phy results
            p = obj.dataObj.convertPhySorting2tIc(obj.spikeSortingFolder);
            label = string(p.label');
            goodU = p.ic(:,label == 'good');

            for s = 1:x %iterate among unique speeds

                fieldName = sprintf('Speed%d', s);

                %Select maximum speed
                At =  A(A(:,5) ==speeds(s),:);

                Bt =  B(B(:,5) ==speeds(s),:);

                [C indexS] = sortrows(At,[2 6 3 4 5]); %2 6 3 4 5

                [Coff indexSo] = sortrows(Bt,[2 6 3 4 5]);

                stimOnT = At(:,1);
                stimOffT = Bt(:,1);
               
                stimDur = mean(-stimOnT+stimOffT);

                %4. Sort directions:
                directimesSorted = C(:,1)';
                Mr = BuildBurstMatrix(goodU,round(p.t/params.binRaster),round(directimesSorted/params.binRaster),round((stimDur+ params.durationWindow)/params.binRaster)); %response matrix
                [MrC]=ConvBurstMatrix(Mr,fspecial('gaussian',[1 params.GaussianLength],3),'same');

                MrNorm = MrC./(sum(MrC,3));
                MrNorm(isnan(MrNorm)) = 0;

                [nT,nN,nB] = size(Mr);

                window_size = [1, round(params.durationWindow/params.binRaster)];
                trialDivision = numel(directimesSorted)/numel(unique(C(:,2)))/numel(unique(C(:,3)))/numel(unique(C(:,4)))...
                    /numel(unique(C(:,5)))/numel(unique(C(:,6)));

                %5. Initialize the maximum mean value and its position
                max_position = zeros(nN,2);
                max_mean_value = zeros(1,nN);
                max_mean_valueB = zeros(1,nN);
                NeuronVals = zeros(nN,nT/trialDivision,10);

                %Baseline = size window

                [Mbd] = BuildBurstMatrix(goodU,round(p.t/params.binRaster),round((directimesSorted-preBase)/params.binRaster),round(preBase/params.binRaster)); %Baseline matrix plus

                MbdC = ConvBurstMatrix(Mbd,fspecial('gaussian',[1 params.GaussianLength],3),'same');

                MbNorm = MbdC./(sum(MbdC,3));
                MbNorm(isnan(MbNorm)) = 0;
                %Merge trials:

                %Real data:
                for u =1:nN
                    % Slide the window over the matrix
                    %unit matrix
                    max_mean_value(u) = -Inf; %General max? not needed
                    max_mean_valueB(u)=-Inf;
                    NeuronRespProfile = zeros(nT/trialDivision,10); %5 columns plus: ofsett, dir, size, speed, lum.

                    k =1;
                    max_position_Trial = zeros(nT/trialDivision,2); %Create 2 matrices, for mean inside max window, and for position of window, for each trial category
                    max_position_TrialB = zeros(nT/trialDivision,2);
                    max_mean_value_Trial = zeros(1,nT/trialDivision);
                    max_mean_value_TrialB = zeros(1,nT/trialDivision);
                    for i = 1:trialDivision:nT %Iterate across trials
                        uM = squeeze(MrNorm(i:i+trialDivision-1,u,:));%Create matrix per unique trial conditions
                        uMB = squeeze(MbNorm(i:i+trialDivision-1,u,:));%Create matrix per unique trial conditions
                        %uMb = Mbd2(i,u);
                        max_mean_value_Trial(k) = -inf;
                        max_mean_value_TrialB(k) = -inf;
                        for j = 1:2:size(uM, 2) - window_size(2) + 1 %Iterate across bins
                            % Extract the sub-matrix
                            sub_matrix = uM(:,j:min(j+window_size(2)-1,end)); %Create matrix of size window per bin
                            sub_matrixB = uMB(:,min(j,size(uMB,2)-window_size(2)+1):min(j+window_size(2)-1,end));
                            % Compute the mean value
                            mean_value = mean(sub_matrix(:)); %Compute mean of each window
                            mean_valueB = mean(sub_matrixB(:)); %Compute mean of each window
                           
                            % %Check if within window there are more than 60%
                            % %of trials with spikes
                            % 
                            % emptyRowsPercentage = sum(all(sub_matrix == 0, 2))/trialDivision;

                            % Update the maximum mean value and its position (a
                            % window is selected across each trial division

                            if mean_value >  max_mean_value_Trial(k) % && emptyRowsPercentage<= 1-params.TrialThreshold
                                max_mean_value_Trial(k) = mean_value;
                                max_position_Trial(k,:) = [i j];
                            end

                            if mean_valueB >  max_mean_value_TrialB(k)
                                max_mean_value_TrialB(k) = mean_valueB;
                                max_position_TrialB(k,:) = [i j];
                            end

                            %%%Check if window is valid for pval

                        end
                        %Save across each trial (in a row) the max bin window
                        %1%. Response
                        NeuronRespProfile(k,1) = max_mean_value_Trial(k);

                        %2%. WindowTrial
                        NeuronRespProfile(k,2) = max_position_Trial(k,1);

                        %3%. WindowBin
                        NeuronRespProfile(k,3) = max_position_Trial(k,2);

                        %4% real spike rate of response
                        NeuronRespProfile(k,4) = mean(Mr(max_position_Trial(k,1):max_position_Trial(k,1)+trialDivision-1,u,...
                            max_position_Trial(k,2):max_position_Trial(k,2)+window_size(2)-1),'all');%max_position_Trial(k,2);

                        BaseResp = mean(Mbd(max_position_TrialB(k,1):max_position_TrialB(k,1)+trialDivision-1,u,...
                            max_position_TrialB(k,2):max_position_TrialB(k,2)+window_size(2)-1),'all');
                       
                        %5%. Resp - Baseline
                        NeuronRespProfile(k,5) =  NeuronRespProfile(k,4)-BaseResp;
                        
                        %Assign visual stats to last columns of NeuronRespProfile. Select
                        %according  to trial (d) the appropiate parameter > directions'offsets' sizes' speeds' freq'
                        NeuronRespProfile(k,6) = C(k*trialDivision,2);
                        NeuronRespProfile(k,7) = C(k*trialDivision,3);
                        NeuronRespProfile(k,8) = C(k*trialDivision,4);
                        NeuronRespProfile(k,9) = C(k*trialDivision,5);
                        NeuronRespProfile(k,10) = C(k*trialDivision,6);
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
                params.N_bootstrap = 5000
                params.speed = 0 %min =1, middle = 2, max=3;
                params.AllSpeeds = true
                params.normBaseline = false
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
                resampled_indicesTi = single(randi(nB, [preBase, params.N_bootstrap]));
                kernel = ones(trialDivision, duration) / (trialDivision * duration); % Normalize for mean


                tic
                for i = 1:params.N_bootstrap
                    % Resample trials with replacement
                    resampled_trials = baseline(resampled_indicesTr(:, i), :,resampled_indicesTi(:, i));
                    for ui = 1:nN
                        % Extract the slice for the current unit (t x b matrix)
                        slice = resampled_trials(:, ui, :);
                        slice = squeeze(slice); % Result is t x b


                        if params.normBaseline
                            %Normalize slice trials
                            slice = slice./sum(slice,2);
                            slice(isnan(Norm)) = 0;

                        end

                        % Compute the mean using 2D convolution
                        means = conv2(slice, kernel, 'valid'); % 'valid' ensures the window fits within bounds

                        % Find the maximum mean in this slice
                        boot_means(i, ui) = max(means(:));
                    end
                end
                toc

                %[bootstats] = get_bootstrapped_equalsamples(data,nruns,num_trials,param)

                if params.normBaseline
                    [respVal,respVali]= max(NeuronResp.(fieldName).NeuronVals(:,:,1),[],2);
                else
                    [respVal,respVali]= max(NeuronResp.(fieldName).NeuronVals(:,:,4),[],2);
                end

                Mr = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted)/bin),round((NeuronResp.(fieldName).stimDur+duration)/bin)); %response matrix
                %%% Calculate p-value 

                pvalsResponse = zeros(1,nN);
                ZScoreU = zeros(1,nN);

                for u = 1:nN
                    pvalsResponse(u) = mean(boot_means(:,u)>=respVal(u));
                    ZScoreU(u) = (respVal(u)-mean(boot_means(:,u)))/(std(boot_means(:,u))+1/(params.N_bootstrap*trialDivision));

                    posTr = NeuronResp.(fieldName).NeuronVals(u,respVali(u),2);
                    posBin = NeuronResp.(fieldName).NeuronVals(u,respVali(u),3);

                    maxWindow = squeeze(Mr(posTr:posTr+trialDivision-1,u,posBin:posBin+duration-1));

                    emptyRows = sum(all(maxWindow == 0, 2));

                    if emptyRows/trialDivision >= params.trialThreshold 
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
