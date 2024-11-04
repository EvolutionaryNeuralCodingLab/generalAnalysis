function []=imageConvert12bitTo8(d,varargin)
%default options
includeSubDirs=true;
originalImageType='tif';
ouputImageType='jpg';

%print out default arguments and values if no inputs are given
if nargin==0
    defaultArguments=who;
    for i=1:numel(defaultArguments)
        eval(['defaultArgumentValue=' defaultArguments{i} ';']);
        disp([defaultArguments{i} ' = ' num2str(defaultArgumentValue)]);
    end
    return;
end

%Collects all options
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end

[strPath,fileName]=fileparts(d);
allFolders=strsplit(genpath([strPath filesep fileName]),';');

for i=1:numel(allFolders)
    imageFiles=dir([allFolders{i} filesep '*.' originalImageType]);
    for j=1:numel(imageFiles)
        I=imread([allFolders{i} filesep imageFiles(j).name]);
        [~,fileName,ext]=fileparts(imageFiles(j).name);
        newImageName=[fileName '_8Bit'];
        if size(I,3)==4
            I=I(:,:,1:3);
        end
        J = imadjust(I,stretchlim(I),[0 1]);
        imwrite(uint8(double(J)/2^8),[allFolders{i} filesep newImageName '.' ouputImageType]);
    end
end
