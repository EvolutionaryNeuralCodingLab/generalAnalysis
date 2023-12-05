function [obj,metaDataExistFlag] = loadMetaData(obj,dataRecordingFolder) %gridSorter
if nargin==1
    dataRecordingFolder=obj.sortingDir;
end
if ~exist([dataRecordingFolder filesep 'metaData.mat'],'file')
    metaDataExistFlag=0;
    disp('Meta data file does not exist!!!');
else
    disp('Setting gridSorted properties according to metaData file');
    metaData=load([dataRecordingFolder filesep 'metaData.mat']);
    for i=1:size(metaData.props,1)
        if isprop(obj,metaData.props{i,1}) & ~strcmp(metaData.props{i,1},'dataRecordingObj')
            obj.(metaData.props{i,1})=metaData.props{i,2};
        else
            disp(['Property ' metaData.props{i,1} ' not loaded!']);
        end
    end
    metaDataExistFlag=1;
end

end