
%tileFile = '/data4/REMA/region_06_luitpold_coast/mosaic_reg_qc_feather2/40m/36_18_40m_dem.mat';
dbasefile  = 'rema_strips_8m_wqc_cs2bias.mat'; % database file
tileDir ='/data4/REMA/region_27_qml_south/mosaic_reg_qc_feather2/40m'; %directory containing tiles
 
tileFiles = dir([tileDir,'/*dem.mat']);
fileDates = [tileFiles.datenum];
tileFiles = {tileFiles.name};
tileFiles = cellfun( @(x) [tileDir,'/',x], tileFiles,'uniformoutput',0);

% use this if to only look at recently created files
% n = fileDates > datenum('31-Jan-2018 12:00:00'); 
% tileFiles = tileFiles(n);

meta=load(dbasefile);

i=1;
for i=1:length(tileFiles)

    fprintf('tile %s,  %d of %d\n',tileFiles{i},i,length(tileFiles));
    
    qcTile(tileFiles{i},meta);

end