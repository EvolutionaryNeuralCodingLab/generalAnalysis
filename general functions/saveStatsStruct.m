function outPath = saveStatsStruct(vsObj, statsName, blocks, producerParams, doSave)
% saveStatsStruct  Persist swarm-panel statistics next to the exported figures.
%
%   outPath = saveStatsStruct(vsObj, statsName, blocks, producerParams, doSave)
%
%   Writes <statsName>_statsStruct.mat into the SAME directory printFig writes
%   its images to, so an analysis's numbers travel with its figures. File them
%   into the 'Figure N _ ...' subfolder together and
%   makeSupplementaryStatsTable(N) will find them.
%
%   Mirrors what AllExpAnalysis already does for its own summaryStats, so there
%   is one convention across the pipeline rather than one per analysis.
%
%   INPUTS
%     vsObj          - any VStimAnalysis-derived object, used only to resolve
%                      the output directory the same way printFig does
%     statsName      - base filename, normally matching the figure name passed
%                      to printFig so the pair is obvious on disk
%     blocks         - Nx3 cell array, one row per reported comparison:
%                        {label, stats, unitName}
%                      label    : char, what the panel shows
%                                 (e.g. 'Spatial tuning index (MB vs SB)')
%                      stats    : struct with .groups and .diffs, as returned by
%                                 plotSwarmBootstrapWithComparisons (4th output)
%                                 or assembled via computeSwarmStats
%                      unitName : 'neurons' | 'insertions' | 'animals', what one
%                                 observation is. buildStatsTable reports this
%                                 in the nUnitIs column because n means different
%                                 things in different panels.
%     producerParams - the calling analysis's params struct, stored so the
%                      exclusion-criteria column can be derived from what was
%                      actually run rather than hand-typed
%     doSave         - gate (pass params.PaperFig) so stats are written on
%                      exactly the runs that export figures
%
%   OUTPUT
%     outPath - full path written, or '' when doSave is false
%
%   See also: computeSwarmStats, buildStatsTable, makeSupplementaryStatsTable

outPath = '';

% Only persist on the runs that also export figures
if ~doSave
    return
end

% Shape check: catching this here beats discovering it when the table is built
assert(iscell(blocks) && size(blocks,2) == 3, ...
    'saveStatsStruct: blocks must be an Nx3 cell {label, stats, unitName}.');

% Drop blocks whose stats are empty (an analysis branch that did not run),
% rather than writing rows that would appear as blanks in the table
isEmptyBlock = false(size(blocks,1),1);
for i = 1:size(blocks,1)
    st = blocks{i,2};
    isEmptyBlock(i) = isempty(st) || ~isstruct(st) || ...
        (~isfield(st,'groups') && ~isfield(st,'diffs'));
end
if any(isEmptyBlock)
    fprintf('  saveStatsStruct: skipping %d empty block(s): %s\n', ...
        sum(isEmptyBlock), strjoin(blocks(isEmptyBlock,1)', ', '));
end
blocks = blocks(~isEmptyBlock, :);

if isempty(blocks)
    warning('saveStatsStruct:noBlocks', ...
        'No non-empty stats blocks for "%s"; nothing written.', statsName);
    return
end

% Resolve printFig's output directory the same way printFig does
% (VStimAnalysis.printFig: split visualStimPlotsFolder, keep the first two
% parts, append 'Paper_figs'). Derived rather than duplicated as a literal so
% the two cannot drift apart.
parts  = split(vsObj.visualStimPlotsFolder, filesep);
outDir = [fullfile(parts{1:2}), filesep, 'Paper_figs'];
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Assemble and write
S.blocks         = blocks;                             % the reported comparisons
S.producerParams = producerParams;                     % what was run, for exclusion text
S.producer       = statsName;                          % provenance
S.savedOn        = datetime('now');                    % when

outPath = fullfile(outDir, [statsName '_statsStruct.mat']);
save(outPath, '-struct', 'S');

fprintf('Saved stats (%d block(s)) to: %s\n', size(blocks,1), outPath);

end
