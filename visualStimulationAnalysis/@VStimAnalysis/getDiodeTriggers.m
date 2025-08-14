%Extract the frame flips from the diode signal
function results=getDiodeTriggers(obj,params)
arguments (Input)
    obj
    %extractionMethod (1,1) string {mustBeMember(extractionMethod,{'diodeThreshold','digitalTriggerDiode'})} = 'diodeThreshold';
    params.extractionMethod string = 'diodeThreshold' %the method used to extract frame flipes-{'diodeThreshold','digitalTriggerDiode'}
    params.analogDataCh = 1
    params.overwrite logical = false %if true overwrites results
    params.analysisTime = datetime('now') %extract the time at which analysis was performed
    params.inputParams = false %if true - prints out the iput parameters so that it is clear what can be manipulated in the method
end
if params.inputParams,disp(params),return,end

%load previous results if analysis was previuosly performed and there is no need to overwrite otherwise continue
results = obj.isOutputAnalysis(obj.getAnalysisFileName,params.overwrite,nargout==1);
if ~isempty(results), return, end

fprintf('Extracting diode signal from analog channel #%d\n',params.analogDataCh);
switch params.extractionMethod
    case "diodeThreshold"

        if ~any(isempty([obj.sessionStartTime,obj.sessionEndTime]))
            [A,t_ms]=obj.dataObj.getAnalogData(params.analogDataCh,obj.sessionStartTime,obj.sessionEndTime-obj.sessionStartTime); %extract diode data for entire recording

            Th=mean(A(1:100:end));
            diodeUpCross=t_ms(A(1:end-1)<Th & A(2:end)>=Th)+obj.sessionStartTime;
            diodeDownCross=t_ms(A(1:end-1)>Th & A(2:end)<=Th)+obj.sessionStartTime;
        else
            disp('Missing start and end times!!! Please run getSessionTime before extracting triggers');
            return;
        end
    case "digitalTriggerDiode"

        switch obj.trialType

            case 'videoTrials'

                if ~any(isempty([obj.sessionStartTime,obj.sessionEndTime]))

                    if all(obj.trialType == 'videoTrials')
                        expectedFlipsperTrial = unique(obj.VST.nFrames);
                        speeds = obj.VST.speeds;
                        framesNspeed = zeros(2,length(speeds));
                        framesNspeed(1,:) =  speeds;
                        %Works for two speeds
                        framesNspeed(2,speeds == min(speeds)) = max(expectedFlipsperTrial);
                        framesNspeed(2,speeds ==  max(speeds)) = min(expectedFlipsperTrial);
                    else
                        framesNspeed = ones(2,1);
                    end

                    t = obj.dataObj.getTrigger;
                    trialOn = t{3}(t{3} > obj.sessionStartTime & t{3} < obj.sessionEndTime);
                    trialOff = t{4}(t{4} > obj.sessionStartTime & t{4} < obj.sessionEndTime);
                    interDelayMs = obj.VST.interTrialDelay*1000;

                    %[A,t_ms]=obj.dataObj.getAnalogData(params.analogDataCh,trialOn(1)-interDelayMs/2,trialOff(end)-trialOn(1)+interDelayMs); %extract diode data for entire recording
                    [A,t_ms]=obj.dataObj.getAnalogData(params.analogDataCh,trialOn(1),trialOff(end)-trialOn(1)+interDelayMs); %extract diode data for entire recording

                    DiodeCrosses = cell(2,numel(trialOn));
                    moreCross =0;
                    trialMostcross=inf;
                    intTrials =[];
                    iMC =0;
                    intrialsNum = 0;
                    trialFail =0;
                    failedTrials =[];
                    for i =1:length(trialOff)

                        startSnip  = round((trialOn(i)-trialOn(1))*(obj.dataObj.samplingFrequencyNI/1000))+1;
                        endSnip  = round((trialOff(i)-trialOn(1)+interDelayMs)*(obj.dataObj.samplingFrequencyNI/1000));

                        if endSnip>length(A)
                            signal = squeeze(A(startSnip:end));
                            t_msS = t_ms(startSnip:end);
                        else
                            signal =squeeze(A(startSnip:endSnip));
                            t_msS = t_ms(startSnip:endSnip);
                        end
                        fDat=medfilt1(signal,15);
                        Th=mean(fDat(1:100:end));
                        stdS = std(fDat(1:100:end));
                        sdK = 0;
                        upTimes=t_msS(fDat(1:end-1)<Th-sdK*stdS & fDat(2:end)>=Th+sdK*stdS)+trialOn(1);%+interDelayMs/2; %get real recording times
                        downTimes=t_msS(fDat(1:end-1)>Th+sdK*stdS  & fDat(2:end)<=Th-sdK*stdS )+trialOn(1);%+interDelayMs/2;

                        % Filter crossings: Remove those too close together (e.g., < 50 ms)
                        minISI = 2*floor(1000/obj.VST.fps);  % ms
                        filterISI = @(x) x([true, diff(x) > minISI]);

                        try
                            DiodeCrosses{1,i} = filterISI(upTimes);
                            DiodeCrosses{2,i} = filterISI(downTimes);
                        catch
                            DiodeCrosses{1,i} = upTimes;
                            DiodeCrosses{2,i} = downTimes;
                        end

                        if (length(DiodeCrosses{1,i}) + length(DiodeCrosses{2,i}))*1.1 < framesNspeed(2,i)
                            %if the number of calculated frames is less than 10%
                            %then perform an interpolation with the
                            %first and last cross

                            if (length(DiodeCrosses{1,i}) + length(DiodeCrosses{2,i})) < framesNspeed(2,i)*0.5
                                %%Diode failure. Use digital triggers
                                %%and interpolate.
                                diodeAll = zeros(1,2);
                                diodeAll(1) = trialOn(i);
                                diodeAll(2) = trialOff(i);
                                ind = 2;
                                trialFail = trialFail +1;
                                failedTrials = [failedTrials i];
                            else
                                [~, ind]=min([DiodeCrosses{1,i}(1)  DiodeCrosses{2,i}(1)]); %check if trial starts with up or down cross
                                diodeAll = sort([DiodeCrosses{1,i} DiodeCrosses{2,i}]);
                            end

                            %DiodeInterp = linspace(diodeAll(1),diodeAll(end),framesNspeed(2,i));
                            DiodeInterp = diodeAll(1):1000/obj.VST.fps: diodeAll(1) + (framesNspeed(2,i)-1)*(1000/obj.VST.fps);
                            if ind == 2 %Trial starts with down cross
                                DiodeCrosses{2,i} = DiodeInterp(1:2:end);
                                DiodeCrosses{1,i} = DiodeInterp(2:2:end);
                            else
                                DiodeCrosses{2,i} = DiodeInterp(2:2:end);
                                DiodeCrosses{1,i} = DiodeInterp(1:2:end);
                            end

                            intTrials = [intTrials i];
                            intrialsNum = intrialsNum+1;

                        end

                        if (length(DiodeCrosses{1,i}) + length(DiodeCrosses{2,i}))>framesNspeed(2,i)
                            %if there are more crosses than there
                            %should be
                            moreCross = moreCross+1;
                            if trialMostcross>(length(DiodeCrosses{1,i}) + length(DiodeCrosses{2,i})) - framesNspeed(2,i)
                                trialMostcross = (length(DiodeCrosses{1,i}) + length(DiodeCrosses{2,i})) - framesNspeed(2,i);
                                iMC = i;
                            end
                        end
                    end

                    diodeUpCross=cell2mat(DiodeCrosses(1,:));
                    diodeDownCross=cell2mat(DiodeCrosses(2,:));

                    fprintf('%d trials out of %d have little or no diode signal, assuming diode failure but correct fliping in trials:',trialFail,length(trialOff))
                    fprintf('%.0f ', failedTrials);
                    fprintf('\n');
                    fprintf('%d trials have excess crossings out of %d; trial %d has the most excess crossings: %d',moreCross,length(trialOff),iMC,trialMostcross)
                    fprintf('\n');
                    fprintf('%d Interpolated trials (out of %d) with more than 10%% of crosses missing: ',intrialsNum,length(trialOff));
                    fprintf('%.0f ', intTrials);
                    fprintf('\n');
                    %Test
                    figure;plot(squeeze(fDat));
                    hold on;xline((DiodeCrosses{1,i} - trialOn(i))*(obj.dataObj.samplingFrequencyNI/1000))
                    xline((DiodeCrosses{2,i}  - trialOn(i))*(obj.dataObj.samplingFrequencyNI/1000),'r')
                    % %xline((trialOff(i)-trialOn(1))*(obj.dataObj.samplingFrequencyNI/1000),'b')
                    % figure;plot(1:length(DiodeCrosses{1,i}) + length(DiodeCrosses{2,i})-1,diff(sort([DiodeCrosses{1,i} DiodeCrosses{2,i}])))
                else
                    disp('Missing start and end times!!! Please run getSessionTime before extracting triggers');
                end

            case 'imageTrials'

                t = obj.dataObj.getTrigger;
                trialOn = t{3}(t{3} > obj.sessionStartTime & t{3} < obj.sessionEndTime);
                trialOff = t{4}(t{4} > obj.sessionStartTime & t{4} < obj.sessionEndTime);
                interDelayMs = obj.VST.interTrialDelay*1000;

                %[A,t_ms]=obj.dataObj.getAnalogData(params.analogDataCh,trialOn(1)-interDelayMs/2,trialOff(end)-trialOn(1)+interDelayMs); %extract diode data for entire recording
                [A,t_ms]=obj.dataObj.getAnalogData(params.analogDataCh,trialOn(1),trialOff(end)-trialOn(1)+interDelayMs); %extract diode data for entire recording

                DiodeCrosses = cell(2,numel(trialOn));
                moreCross =0;
                trialMostcross=inf;
                intTrials =[];
                iMC =0;
                intrialsNum = 0;
                trialFail =0;
                failedTrials =[];
                for i =1:length(trialOff)

                    startSnip  = round((trialOn(i)-trialOn(1)-100)*(obj.dataObj.samplingFrequencyNI/1000))+1;
                    endSnip  = round((trialOff(i)-trialOn(1)+interDelayMs/2)*(obj.dataObj.samplingFrequencyNI/1000));

                    if endSnip>length(A)
                        signal = squeeze(A(startSnip:end));
                        t_msS = t_ms(startSnip:end);
                    elseif startSnip<1
                        signal =squeeze(A(1:endSnip));
                        t_msS = t_ms(1:endSnip);
                    else
                        signal =squeeze(A(startSnip:endSnip));
                        t_msS = t_ms(startSnip:endSnip);
                    end



                    fDat=medfilt1(signal,(obj.VST.stimDuration/4)*1000);
                    Th=mean(fDat(1:100:end));
                    stdS = std(fDat(1:100:end));
                    sdK = 0;
                    upTimes=t_msS(fDat(1:end-1)<Th-sdK*stdS & fDat(2:end)>=Th+sdK*stdS)+trialOn(1)+100;%+interDelayMs/2; %get real recording times
                    downTimes=t_msS(fDat(1:end-1)>Th+sdK*stdS  & fDat(2:end)<=Th-sdK*stdS )+trialOn(1)+100;%+interDelayMs/2;

                    if length(upTimes) >1 || length(downTimes)>1
                        upTimes=upTimes(1);
                        downTimes = downTimes(1);
                    end

                    if length(upTimes) <1 || length(downTimes)<1
                        2+2
                    end

                    DiodeCrosses{1,i} = upTimes;
                    DiodeCrosses{2,i} = downTimes;


                end

                figure;plot(squeeze(fDat));
                hold on;xline((DiodeCrosses{1,i} - trialOn(i))*(obj.dataObj.samplingFrequencyNI/1000))
                xline((DiodeCrosses{2,i}  - trialOn(i))*(obj.dataObj.samplingFrequencyNI/1000),'r')

                diodeUpCross=cell2mat(DiodeCrosses(1,:));
                diodeDownCross=cell2mat(DiodeCrosses(2,:));

        end
end
results.diodeUpCross = diodeUpCross;
results.diodeDownCross = diodeDownCross;

%save results in the right file
fprintf('Saving results to file.\n');
save(obj.getAnalysisFileName,'params','diodeUpCross','diodeDownCross','Th');
end