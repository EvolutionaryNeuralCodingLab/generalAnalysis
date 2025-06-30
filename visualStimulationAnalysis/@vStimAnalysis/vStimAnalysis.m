classdef vStimAnalysis

    properties
        vStimClass(1,1) string {mustBeMember(vStimClass,...
            {"rectGrid","linearlyMovingBall","linearlyMovingBarWidth","linearlyMovingDotsWrapped","rectNoiseGrid"})} = "";
    end

    methods (Hidden)

        %class constructor
        function obj=vStimAnalysis(vStimClassValue) 
            obj.vStimClass=vStimClassValue;
        end
    end

    methods

        function [upCross,downCross]=getDiodeTriggers(dataObj,analogDataCh,extractionMethod)
            arguments
                dataObj 
                extractionMethod (1,:) string {mustBeMember(extractionMethod,...
                    {"diodeThreshold","digitalTriggerDiode"})} = "diodeThreshold";
            end

            switch extractionMethod
                case "diodeThreshold"
                    [A,t_ms]=dataObj.getAnalogData(analogDataCh,0,analogDataObj.recordingDuration_ms);

                    Th=mean(A(1:100:end));
                    upCross=t_ms(A(1:end-1)<Th & A(2:end)>=Th);
                    downCross=t_ms(A(1:end-1)>Th & A(2:end)<=Th);

                case "digitalTriggerDiode"


            end



        end

    end
end