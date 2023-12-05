%function modified by MSI
function export_gridSorter(varargin)
% export_csv_(hObject, event)
% if nargin<2, 
fZeroIndex = 0; %zero-based index export (disable to export as matlab index starting from 1)

% S0 = get(0, 'UserData');
if nargin==2
    [S0, P, S_clu] = get0_();
elseif nargin==1
    P = varargin{1};
    vcFile_prm = P.vcFile_prm;
    S0 = load_cached_(P, 0);
    if isempty(S0), fprintf(2, 'Cannot find _jrc.mat.\n'); return; end %exit if file doesn't exist
    P = S0.P;
end

% vcFile_clu = subsFileExt(P.vcFile_prm, '_clu.mat');
% Sclu = load(vcFile_clu); %load Sclu    
% if isfield(Sclu, 'Sclu'), Sclu = Sclu.Sclu; end

if isfield(S0, 'S_clu')
    viClu = double(S0.S_clu.viClu);
else
    fprintf(2, 'Cannot find S_clu.\n');
end

nClu=numel(S0.S_clu.cviSpk_clu);
t=cell(nClu,1);
ic=[];
neuCounter=1;
spkCounter=0;
prevSite=-1;
for iClu=1:nClu
    if ~isempty(S0.S_clu.cviSpk_clu{iClu})
        viSpk_clu1 = S0.S_clu.cviSpk_clu{iClu};% 
        iSite_clu1 = S0.S_clu.viSite_clu(iClu);
        if iSite_clu1~=prevSite
            neuCounter=1;
        end
        
        vlSpk_clu1 = iSite_clu1 == S0.viSite_spk(viSpk_clu1);
        %viSite_clu1 = S0.S_clu.P.miSites(:,iSite_clu1);
        
        ic=[ic [iSite_clu1;neuCounter;spkCounter+1;spkCounter+numel(viSpk_clu1)]];
        t{iClu} = S0.viTime_spk(viSpk_clu1(vlSpk_clu1));
        spkCounter=spkCounter+numel(t{iClu});
        neuCounter=neuCounter+1;
        prevSite=iSite_clu1;
    end
end
t=double(cell2mat(t)')./(P.sRateHz/1000);

fileName=subsFileExt_(P.vcFile_prm, '_gridSorter.mat');
classification=S0.S_clu.csNote_clu;
avgRawWF=permute(S0.S_clu.trWav_raw_clu,[3 2 1]);
avgSpkWF=permute(S0.S_clu.tmrWav_clu,[3 2 1]);

save(fileName,'t','ic','avgRawWF','avgSpkWF','classification');
fprintf('wrote grid sorter spike data to %s\n', fileName);
end %func
