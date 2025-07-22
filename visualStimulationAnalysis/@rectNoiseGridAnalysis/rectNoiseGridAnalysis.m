classdef rectNoiseGridAnalysis < VStimAnalysis

    properties

    end

    properties (Constant)
        trialType = 'imageTrials'
    end

    methods (Hidden)
        %class constructor - name of class should be identical to the visual stimulation with the addition of Analysis
        function [obj] = rectNoiseGridAnalysis(dataObj)
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
        %list of methods in seperate m files
        plotSpatialTuningLFP(obj,params)

        function results=getReceptiveFields(obj,params)
            arguments
                obj
                params.win = [100,500] %[1x2]-[ms] the time window [start,end] post image flip for analysis
                params.bin = 10 %[ms] - bin size temporal resolution for analysis
                params.overwrite logical = false %if true overwrites results
                params.analysisTime = datetime('now') %extract the time at which analysis was performed   
                params.inputParams = false %if true - prints out the iput parameters so that it is clear what can be manipulated in the method
            end

            T=obj.getSessionTime;
            s=obj.getSpikeTIcData; %load spike sorting
            SD=obj.getSyncedDiodeTriggers; %get synched triggers

            V = ones(size(obj.VST.stimSequence,1),ceil((T.sessionEndTime-T.sessionStartTime+params.win(2))/params.bin))*obj.VST.visualFieldBackgroundLuminance;
            for i=1:(numel(SD.stimOnFlipTimes)-1)
                p=round((SD.stimOnFlipTimes(i)-T.sessionStartTime)/params.bin):round((SD.stimOnFlipTimes(i+1)-T.sessionStartTime)/params.bin);
                V(:,p)=obj.VST.stimSequence(:,i)*ones(1,numel(p));
            end
            binsInFrame=round(obj.VST.stimDuration*1000/params.bin);
            V(:,p(end)+1:p(end)+binsInFrame)=obj.VST.stimSequence(:,end)*ones(1,binsInFrame); %add last frames based on the duration of each noise stim
            V=V-obj.VST.visualFieldBackgroundLuminance;

            delay=0;
            M=squeeze(BuildBurstMatrix(s.ic,round((s.t-delay)/params.bin),round(T.sessionStartTime/params.bin),ceil((T.sessionEndTime-T.sessionStartTime+params.win(2))/params.bin)));
            
            %plot activity with Stims
            plot((1:size(M,2))*params.bin,mean(M,1),'k');hold on;line([SD.stimOnFlipTimes-T.sessionStartTime;SD.stimOnFlipTimes-T.sessionStartTime],ylim,'Color','r');
            
            t=((1:size(V,2))-size(V,2)/2)*params.bin;
            for i=1:size(M,1)
                if sum(M(i,:))>10
                    C=convn(V',M(i,:)',"same")';
                    plot(t,mean(C,1));
                    %imagesc(C);
                    pause;
                end
            end

            pos2X
            pos2Y
            obj.VST.rectData



        end

        %Methods in this m file
        function results=getSyncedDiodeTriggers(obj,varargin)
            %Call superclass function but analysis for this visual stimulation should only be done on the on Flips
            %when using the function for receptive field mapping - in other scenarios there may be off flips in which case
            %a condition should be added here
            results=getSyncedDiodeTriggers@VStimAnalysis(obj,'analyzeOnlyOnFlips',true,'ignoreNLastFlips',1,varargin{:}); 
        end
          
    end
end