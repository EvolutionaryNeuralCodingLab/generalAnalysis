function folderLocation=findFolderInExperiment(rootFolder,folderNamePart,params)
%find a specific folder in the experiment based on its name and a root folder (searches also parent folders)
arguments
    rootFolder char % the root folder for search
    folderNamePart char %a string that must exist in the begining of the folder name
    params.nParentFolders2Check = 2 % the number of parent folders to search for the folderNamePart folder relative to rootFolder
    params.inputParams = false% if true - prints out the iput parameters so that it is clear what can be manipulated in the method
end

folderFound=false;

%find visual stimulation folder
tmpDir=dir([rootFolder filesep folderNamePart '*']);
tmpDir = tmpDir([tmpDir.isdir]);  % Keep only folders
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
    fprintf('%s folder was not found!!! Making a folder in the same path as the data folder. \n',folderNamePart);
    mkdir([rootFolder filesep folderNamePart]);
    folderLocation=[rootFolder filesep folderNamePart];
else
    fprintf('%s folder found: %s\n',folderNamePart,folderLocation);
end


end