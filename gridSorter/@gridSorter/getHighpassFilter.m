function obj=getHighpassFilter(obj,samplingFrequency)
%calculate the highpass filter for sorting
if ~isempty(obj.dataRecordingObj)
    samplingFrequency=obj.dataRecordingObj.samplingFrequency;
elseif nargin==2
    disp(['Building filter with sampling freq. :' num2str(samplingFrequency)]);
    %answer = inputdlg(prompt,dlg_title)
else
    disp('No data recorded object exists, either add object or enter sampling freq.as a second argument');
    return;
end

obj.filterObj=filterData(samplingFrequency);
if obj.simpleButter
    obj.filterObj.highPassCutoff=obj.highPassCutoff;
    obj.filterObj.lowPassCutoff=obj.lowPassCutoff;
    obj.filterObj.filterOrder=obj.filterOrder;
else
    obj.filterObj.highPassPassCutoff=obj.filterHighPassPassCutoff;
    obj.filterObj.highPassStopCutoff=obj.filterHighPassStopCutoff;
    obj.filterObj.lowPassPassCutoff=obj.filterLowPassPassCutoff;
    obj.filterObj.lowPassStopCutoff=obj.filterLowPassStopCutoff;
    obj.filterObj.attenuationInHighpass=obj.filterAttenuationInHighpass;
    obj.filterObj.attenuationInLowpass=obj.filterAttenuationInLowpass;
    obj.filterObj.rippleInPassband=obj.filterRippleInPassband;
end
obj.filterObj.filterDesign=obj.filterDesign;
obj.filterObj.padding=true;
%obj.filterObj.highPassCutoff=obj.filterHighPassPassCutoff;
%obj.filterObj.lowPassCutoff=obj.filterLowPassPassCutoff;
%obj.filterObj.filterOrder=8;

obj.filterObj=obj.filterObj.designBandPass;

end
