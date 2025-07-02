classdef vStimAnalysis

    properties
        vStimClass(1,1) string {mustBeMember(vStimClass,{'','rectGrid','linearlyMovingBall','linearlyMovingBarWidth','linearlyMovingDotsWrapped','rectNoiseGrid'})}='';
        diodeUpCross
        diodeDownCross
    end

    methods (Hidden)

        %class constructor
        function obj=vStimAnalysis(vStimClassValue)
            if isempty(vStimClassValue)
                disp('No visual class was entered, please define the vStimClass property.');
            else
                obj.vStimClass=vStimClassValue;
            end
        end
    end

    methods

        function [upCross,downCross]=getDiodeTriggers(obj,dataObj,analogDataCh,extractionMethod)
            arguments (Input)
                obj
                dataObj
                analogDataCh
                extractionMethod (1,:) string {mustBeMember(extractionMethod,...
                    {'diodeThreshold','digitalTriggerDiode'})} = 'diodeThreshold';
            end

            switch extractionMethod
                case "diodeThreshold"
                    [A,t_ms]=dataObj.getAnalogData(analogDataCh,0,dataObj.recordingDuration_ms); %extract diode data for entire recording 

                    Th=mean(A(1:100:end));
                    obj.diodeUpCross=t_ms(A(1:end-1)<Th & A(2:end)>=Th);
                    obj.diodeDownCross=t_ms(A(1:end-1)>Th & A(2:end)<=Th);

                case "digitalTriggerDiode"


            end



        end

    end
end