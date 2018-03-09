function batch_qcTiles(tileNames)

%tileFile = '/data4/REMA/region_06_luitpold_coast/mosaic_reg_qc_feather2/40m/36_18_40m_dem.mat';
dbasefile  = 'V:\pgc\data\scratch\claire\repos\setsm_postprocessing_pgc\arcticDEMdatabase_2m.mat'; % database file
tileDir = 'V:\pgc\data\elev\dem\setsm\ArcticDEM\mosaic\2m_v3.1_tileqc'; %directory containing tiles
%tileName = '31_39';

tileNames = strsplit(tileNames,',');

for i=1:length(tileNames);
    tileFiles{i} = [tileDir,'\',char(tileNames(i)),'\',char(tileNames(i)),'_8m_dem.mat'];
end
%fileDates = [tileFiles.datenum];
%tileFiles = {tileFiles.name};
%tileFiles = cellfun( @(x) [tileDir,'\',tileName,'\',x], tileFiles,'uniformoutput',0);

% use this if to only look at recently created files
% n = fileDates > datenum('31-Jan-2018 12:00:00'); 
% tileFiles = tileFiles(n);

fprintf('Loading db\n');
meta=load(dbasefile);

i=1;
for i=1:length(tileFiles)

    if exist(tileFiles{i},'file')
        fprintf('tile %s, %d of %d\n',tileFiles{i},i,length(tileFiles));
        qcTile(tileFiles{i},meta);
    end

end