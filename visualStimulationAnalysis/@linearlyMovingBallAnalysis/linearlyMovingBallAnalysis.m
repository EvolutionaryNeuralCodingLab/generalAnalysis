classdef linearlyMovingBallAnalysis < VStimAnalysis

    properties

    end

    methods (Hidden)
        %class constructor - name of class should be identical to the visual stimulation with the addition of Analysis
        function [obj] = linearlyMovingBallAnalysis(dataObj)
            obj = obj.initialize(dataObj);
        end
    end

    methods

    end
end