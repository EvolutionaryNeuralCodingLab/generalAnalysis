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

        function results = setUpAnalysis(obj, params)

            arguments (Input)
                obj
                params.speedChange string = 'OneSpeed'
                params.overwrite logical = false
                params.analysisTime = datetime('now')
                params.inputParams = false
            end
            if params.inputParams,disp(params),return,end

            obj.getSessionTime;

            try
                DiodeCrossings = obj.getSyncedDiodeTriggers;
            catch
                obj.getDiodeTriggers("extractionMethod",'digitalTriggerDiode','overwrite',true); 
                DiodeCrossings = obj.getSyncedDiodeTriggers;
            end
            allDiodeCross = sort([DiodeCrossings.diodeUpCross DiodeCrossings.diodeDownCross]);
            t = obj.dataObj.getTrigger;
            trialOn = t{3}(t{3} > obj.sessionStartTime & t{3} < obj.sessionEndTime);
            trialOff = t{4}(t{4} > obj.sessionStartTime & t{4} < obj.sessionEndTime);

            stimOn = zeros(1, length(trialOn));
            for i=1:length(trialOn)
                [mn,ind] = min(abs(allDiodeCross - trialOn(i)));
                stimOn(i) = allDiodeCross(ind); 
            end

            % check start and end diode closest to digital trigger on and
            % off
        end

    end

    %static methods that are not directly related to the class are defined here
    methods (Static)
    end
end

%%%1. Get Diode
%%%2. Create A Matrix.
%%%3. Load Kilos0rt and phy results.
%%%4. Create response matrix.
%%%5. Create shuffling analysis