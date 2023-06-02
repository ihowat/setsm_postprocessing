function batch_resetMergeTileBuffers(fileNames)
% undoMergeTileBuffers reverts the tile edge buffers to the original values
%
% undoMergeTileBuffers(fileNames) replaces the merged tile buffers, created
% by mergeTileBuffer, with the original values, setting buffer flags to
% false.

% make cell if not
if ~iscell(fileNames)
    fileNames={fileNames};
end

%%
i=1;
for i =1:length(fileNames)
    check_file = fileNames{i};
    fprintf('\nChecking tile file: %s\n', check_file)
    if ~isfile(check_file)
        error('Source file does not exist')
    end
    resetMergeTileBuffers(check_file)
end
