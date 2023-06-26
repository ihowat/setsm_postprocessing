function batch_insertSubset(tileDir, resolution, annualRootDir, boxShp, pShp, varargin)

tileFiles = dir(sprintf('%s/*_%s.mat', tileDir, resolution));
tileNames = {tileFiles.name};
tileFiles = fullfile({tileFiles.folder}, tileNames);

ntiles = length(tileFiles);
for i = 1:ntiles
    tileFile = tileFiles{i};
    fprintf("Working on tile %d of %d: %s\n", i, ntiles, tileFile);
    call_insertSubset(tileFile, annualRootDir, boxShp, pShp, varargin{:});
end
