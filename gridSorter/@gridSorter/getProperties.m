function [publicProps]=getProperties(obj)
%get the properties of the gridSorter Class
metaClassData=metaclass(obj);
allPropName={metaClassData.PropertyList.Name}';
allPropSetAccess={metaClassData.PropertyList.SetAccess}';
pPublic=find(strcmp(allPropSetAccess, 'public'));
publicProps(:,1)=allPropName(pPublic);

%collect all prop values
for i=1:size(publicProps,1)
    publicProps{i,2}=obj.(publicProps{i,1});
end