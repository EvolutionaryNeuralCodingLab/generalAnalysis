function [data,hand]=getCrossCorr(obj,varargin)
% [data,hand]=getCrossCorr(obj,varargin)
% Function purpose : Calculate spike train cross-correlation in [ms]
%
% Last updated : 20/10/16

obj.checkFileRecording;
parseObj = inputParser;

addParameter(parseObj,'maxLag',500,@isnumeric); %ms
addParameter(parseObj,'bin',1,@isnumeric); %ms
addParameter(parseObj,'hand',[]);
addParameter(parseObj,'fieldPar',[]);
addParameter(parseObj,'saveFileName',[]);
addParameter(parseObj,'overwrite',false,@isnumeric);
addParameter(parseObj,'inputParams',false,@isnumeric);

parseObj.parse(varargin{:});
if parseObj.Results.inputParams
    disp(parseObj.Results);
    return;
end
%evaluate all input parameters in workspace
for i=1:numel(parseObj.Parameters)
    eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
end
%make parameter structure
par=parseObj.Results;
load(obj.files.getSpikeInducedFields,'fieldPar');

[funName] = dbstack;funName=funName.name;funName=strsplit(funName,'.');funName=funName{end};
if isempty(saveFileName)
    saveFileName=obj.files.(funName);
end

obj.checkFileRecording(obj.files.getSpikeInducedFields,'Spike induced fields file missing, please run getSpikeInducedFields');

%load data if existing or contrinue with analysis
if exist(saveFileName,'file') & ~overwrite
    if nargout==1
        data=load(saveFileName);
    else
        disp('Cross-correlation analysis already exists. use overwrite to run again');
    end
    return;
end

%populate grid sorter object
obj=populateGridSorterObj(obj);
S=obj.gridSorterObj.getSortedData({'t','ic'});
nNeu=size(S.ic,2);
%% Main code

%calculate cross correlation
[C]=crossCorrRaster(round(S.t/bin),S.ic,maxLag);

%shuffled cross correlation
nShuffles=5;
CShuf=zeros(nNeu,nNeu,maxLag,nShuffles);
for i=1:nShuffles
    tShuf=S.t(randperm(numel(S.t)));
    for j=1:nNeu
        tShuf(S.ic(3,j):S.ic(4,j))=sort(tShuf(S.ic(3,j):S.ic(4,j)));
    end
    CShuf(:,:,:,i)=crossCorrRaster(round(tShuf/bin),S.ic,maxLag);
end

%save crossCorr CShuf C;

%shuffled cross correlation within sub-populations
pI=find(fieldPar.classIE==2);
pE=find(fieldPar.classIE==3);

%shuffled cross correlation -  only within sub-population of excitatory and inhibitory neurons
CIEShuf=zeros(nNeu,nNeu,maxLag,nShuffles);
for i=1:nShuffles
    tShuf=zeros(1,numel(S.t));
    idx=[];
    for j=pE
        idx=[idx S.ic(3,j):S.ic(4,j)];
    end
    tShuf(idx)=S.t(idx(randperm(numel(idx))));
    
    idx=[];
    for j=pI
        idx=[idx S.ic(3,j):S.ic(4,j)];
    end
    tShuf(idx)=S.t(idx(randperm(numel(idx))));
    for j=1:nNeu
        tShuf(S.ic(3,j):S.ic(4,j))=sort(tShuf(S.ic(3,j):S.ic(4,j)));
    end
    CIEShuf(:,:,:,i)=crossCorrRaster(round(tShuf/bin),S.ic,maxLag);
end


%save all data
fprintf('\nSaving results...');
save(saveFileName,'pI','pE','CIEShuf','C','CShuf','par');

fprintf('Done!\n');

%save(saveFileName,'pMaxField','-append');

