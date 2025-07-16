function folderLocation=findFolderInExperiment(rootFolder,folderNamePart,params)
arguments
    rootFolder char
    folderNamePart char
    params.inputParams = false
    params.nParentFolders2Check = 2
end

folderFound=false;

%find visual stimulation folder
tmpDir=dir([rootFolder filesep folderNamePart '*']);
if isempty(tmpDir) %check if not in the current data folder
    %go one folder back and look for visualStimulation folder
    fileSepTransitions=regexp(rootFolder,filesep); %look for file separation transitions
    if fileSepTransitions(end)==numel(rootFolder) %if last transition appears in the end of the folder remove this transition
        fileSepTransitions(end)=[];
    end
    for i=1:params.nParentFolders2Check %repeat folder search params.nParentFolders2Check folders up
        tmpCurrentFolder=rootFolder(1:fileSepTransitions(end));
        %check parent folder for visual stimulation folder
        tmpDir=dir([tmpCurrentFolder filesep folderNamePart '*']);
        if ~isempty(tmpDir)
            folderLocation=[tmpCurrentFolder filesep tmpDir.name];
            folderFound=true;
        end
        fileSepTransitions(end)=[];
    end
else
    folderLocation=[rootFolder filesep tmpDir.name];
    folderFound=true;
end
if ~folderFound
    fprintf('%s folder was not found!!! Please make a folder in the same path as the data folder. \nNotice the the name of the folder should start with %s\n',folderNamePart,folderNamePart);
    folderLocation=[];
else
    fprintf('%s folder found: %s\n',folderNamePart,folderLocation);
end


end