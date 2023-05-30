function countMergedEdges(tileDir, resolution)

regFiles = dir(sprintf('%s/*_%s_reg.mat', tileDir, resolution));
unregFiles = dir(sprintf('%s/*_%s.mat', tileDir, resolution));
regFiles = fullfile({regFiles.folder}, {regFiles.name});
unregFiles = fullfile({unregFiles.folder}, {unregFiles.name});
unregFiles = setdiff(unregFiles, strrep(regFiles, '_reg.mat', '.mat'));
tileFiles = [unregFiles regFiles];

for i = 1:length(tileFiles)
    tileFile = tileFiles{i};
    [~,tileFname,~] = fileparts(tileFile);

    m = matfile(tileFile);
    m_struct = m;
    m_varlist = who(m_struct);

    mergedTop = ismember('mergedTop', m_varlist) && m_struct.mergedTop;
    mergedBottom = ismember('mergedBottom', m_varlist) && m_struct.mergedBottom;
    mergedLeft = ismember('mergedLeft', m_varlist) && m_struct.mergedLeft;
    mergedRight = ismember('mergedRight', m_varlist) && m_struct.mergedRight;

    count = nnz([mergedTop, mergedBottom, mergedLeft, mergedRight]);

    fprintf('Tile:%s,Count:%d\n', tileFname, count);
end
