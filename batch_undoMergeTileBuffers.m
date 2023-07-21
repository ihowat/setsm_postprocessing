function batch_undoMergeTileBuffers(tileDir, resolution)
% undoMergeTileBuffers reverts the tile edge buffers to the original values
%
% undoMergeTileBuffers(fileNames) replaces the merged tile buffers, created
% by mergeTileBuffer, with the original values, setting buffer flags to
% false.

regFiles = dir(sprintf('%s/*_%s_reg.mat', tileDir, resolution));
unregFiles = dir(sprintf('%s/*_%s.mat', tileDir, resolution));
regFiles = fullfile({regFiles.folder}, {regFiles.name});
unregFiles = fullfile({unregFiles.folder}, {unregFiles.name});
%unregFiles = setdiff(unregFiles, strrep(regFiles, '_reg.mat', '.mat'));
tileFiles = [unregFiles regFiles];

ntiles = length(tileFiles);
for i = 1:ntiles
    tileFile = tileFiles{i};
    fprintf("Working on tile %d of %d: %s\n", i, ntiles, tileFile);
    undoMergeTileBuffersSingle(tileFile, 'all', 'ignoreMergedVars',true)
end
