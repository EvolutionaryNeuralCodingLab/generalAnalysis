classdef linearlyMovingBallAnalysis < VStimAnalysis

    properties

    end

    properties (Constant)
        trialType = 'imageTrials'
    end

    methods (Hidden)
        %class constructor - name of class should be identical to the visual stimulation with the addition of Analysis
        function [obj] = linearlyMovingBallAnalysis(dataObj)
            obj = obj.initialize(dataObj);
        end

        %Get diode 
        
    end

    methods

    end
end