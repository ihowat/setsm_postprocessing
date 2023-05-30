function verifyReg(tileDir, resolution)

tileFiles = dir(sprintf('%s/*_%s_reg.mat', tileDir, resolution));
tileNames = {tileFiles.name};
tileFiles = fullfile({tileFiles.folder}, tileNames);

for i = 1:length(tileFiles)
    tileFile = tileFiles{i};
    [~,tileFname,~] = fileparts(tileFile);

    m = matfile(tileFile);
    m_struct = m;
    m_varlist = who(m_struct);

    regToCOP30 = ismember('regToCOP30', m_varlist) && m_struct.regToCOP30;

    if regToCOP30
        fprintf('Tile:%s,IS_REGISTERED\n', tileFname);
    else
        fprintf('Tile:%s,NOT_REGISTERED\n', tileFname);
    end
end
