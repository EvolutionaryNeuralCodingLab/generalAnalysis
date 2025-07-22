classdef rectGridAnalysis < VStimAnalysis

    properties

    end

    properties (Constant)
        trialType = 'imageTrials'
    end

    methods (Hidden)
        %class constructor - name of class should be identical to the visual stimulation with the addition of Analysis
        function [obj] = rectGridAnalysis(dataObj)
            if nargin==0
                dataObj=[];
            end
            % Call superclass constructor
            obj@VStimAnalysis(dataObj);
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

    end
end