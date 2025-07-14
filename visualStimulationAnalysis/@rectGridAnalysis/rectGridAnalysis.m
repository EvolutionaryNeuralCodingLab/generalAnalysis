classdef rectGridAnalysis < VStimAnalysis

    properties

    end

    properties (Constant)
        trialType = 'imageTrials'
    end

    methods (Hidden)
        %class constructor - name of class should be identical to the visual stimulation with the addition of Analysis
        function [obj] = rectGridAnalysis(dataObj)

            obj = obj.initialize(dataObj);
        end
    end

    methods

    end
end