function ops=makeConfigKiloSort2(sortingDir,nCh,varargin)

ops.datatype = 'bin';  % binary ('dat', 'bin') or 'openEphys'		
ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
ops.parfor              = 1;
ops.verbose             = 1;
ops.showfigures         = 1;

ops.FiltersPerCh        = 3;
ops.nNeighPC            = 12; % visualization only (Phy): number of channnels to mask the PCs, leave empty to skip (12)		
ops.nNeigh              = 16; % visualization only (Phy): number of neighboring templates to retain projections of (16)		

ops.whitening       = 'full'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)		
ops.spkTh           = -4;      % spike threshold in standard deviations (-6)
ops.reorder         = 1;       % whether to reorder batches for drift correction. 
ops.nskip           = 32;  % how many batches to skip for determining spike PCs
ops.whiteningRange      = 32; % number of channels to use for whitening each channel
ops.nSkipCov            = 25; %not sure about this. % compute whitening matrix from every N-th batch (1)

ops.fs = 20000;  
% frequency for high pass filtering (150)
ops.fshigh = 200;   
% minimum firing rate on a "good" channel (0 to skip)
ops.minfr_goodchannels = 0.1; 
% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.Th = [6 12 12];
% ops.Th = [12 6];  
% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot) 
ops.lam = [5 5 5];  

ops.nannealpasses    = 4;        

% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.AUCsplit = 0.9; 
% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.minFR = 1/50; 
% number of samples to average over (annealed from first to second value) 
ops.momentum = [20 400]; 

ops.shuffle_clusters = 1;            % allow merges and splits during optimization (1)		
ops.mergeT           = .1;           % upper threshold for merging (.1)		
ops.splitT           = .1;  

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 30; 
% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8; 
%% danger, changing these settings can lead to fatal errors
% options for determining PCs

ops.criterionNoiseChannels = 0.2; % fraction of "noise" templates allowed to span all channel groups (see createChannelMapFile for more info). 		

% other options for controlling the model and optimization		
ops.Nrank               = 3;    % matrix rank of spike template model (3)		
ops.nfullpasses         = 6;    % number of complete passes through data during optimization (6)		
ops.maxFR               = 20000;  % maximum number of spikes to extract per batch (20000)		
ops.fslow               = 3000;   % frequency for low pass filtering (optional)

ops.Nfilt               = 1024; %64; % max number of clusters
% ops.Nfilt               = round(nCh*ops.FiltersPerCh/32)*32;
ops.nfilt_factor        = 3; % max number of clusters per good channel (even temporary ones)
% ops.Nfilt               = round(nCh*ops.FiltersPerCh/32)*32;  % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)    
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.NT                  = 128*1024+ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). 


ops.scaleproc           = 200;   % int16 scaling of whitened data
ops.nPCs                = 3; % how many PCs to project the spikes into
ops.useRAM              = 0; % not yet available

ops.initialize      = 'fromData'; 
ops.loc_range       = [3  1];  % ranges to detect peaks; plus/minus in time and channel ([3 1])		
ops.long_range      = [30  6]; % ranges to detect isolated peaks ([30 6])		
ops.maskMaxChannels = 5;       % how many channels to mask up/down ([5])		
ops.crit            = .65;     % upper criterion for discarding spike repeates (0.65)		
ops.nFiltMax        = 10000;   % maximum "unique" spikes to consider (10000)		

ops.fracse  = 0.1; % binning step along discriminant axis for posthoc merges (in units of sd)		
ops.epu     = Inf;	
		
ops.ForceMaxRAMforDat   = 20e9; % maximum RAM the algorithm will try to use; on Windows it will autodetect.
dd                  = load('PCspikes2.mat'); % you might want to recompute this from your own data		
ops.wPCA            = dd.Wi(:,1:7);   % PCs 	

ops.trange = [0 Inf];
% ops.NchanTOT    = 252;
ops.NchanTOT    = nCh;
%%
if nargin~=0
    defaultArguments=fields(ops);
    for i=1:numel(defaultArguments)
        eval(['defaultArgumentValue=ops.' defaultArguments{i} ';']);
%         if numel(defaultArgumentValue)==1
%             disp([defaultArguments{i} ' = ' num2str(defaultArgumentValue)]);
%         else
%             fprintf([defaultArguments{i} ' = ']);
%             disp(defaultArgumentValue);
%         end
    end
%     return;
end

%options depending on the input data object
binaryFileName='data4KiloSort.bin';
channelMapFile='chanMap.mat';

ops.fbinary             = fullfile(sortingDir, binaryFileName); % will be created for 'openEphys'		
ops.fproc               = fullfile(sortingDir, 'temp_wh.dat'); % residual from RAM of preprocessed data		
ops.root                = sortingDir; % 'openEphys' only: where raw files are		
% define the channel map as a filename (string) or simply an array		
ops.chanMap             = fullfile(sortingDir, channelMapFile); % make this file using createChannelMapFile.m		
% ops.chanMap = 1:ops.Nchan; % treated as linear probe if unavailable chanMap file		

ops.nt0 = 25;
ops.fig = 0;
%% Overwrite with input values
for i=1:2:length(varargin)
    eval(['ops.' varargin{i} '=' 'varargin{i+1};'])
end

NSKToolBoxMainDir=fileparts(which('identifierOfMainDir4NSKToolBox'));
configText=fscanf(fopen([NSKToolBoxMainDir filesep 'PCspecificFiles' filesep 'kiloSortOpsPath.txt']),'%c');
ops_temp = load(configText);
% ops_gui = ops_temp.saveDat.ops;
ops_gui = ops_temp.ops;
% ops_gui = rmfield(ops_gui, 'fbinary');
ops_gui = rmfield(ops_gui,'fproc');
ops_gui = rmfield(ops_gui,'saveDir');
ops_gui = rmfield(ops_gui,'chanMap');
ops_gui = rmfield(ops_gui,'fshigh');
for fn = fieldnames(ops_gui)'
   ops.(fn{1}) = ops_gui.(fn{1});
end
end
