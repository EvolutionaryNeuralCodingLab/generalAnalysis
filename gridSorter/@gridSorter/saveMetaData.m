function saveMetaData(obj)
[props]=getProperties(obj);
save([obj.sortingDir filesep 'metaData.mat'],'props');