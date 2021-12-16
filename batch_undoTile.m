function batch_undoTile(tileDir,strt,inc)

fileNames=dir([tileDir,'/*_10m_reg.mat']);
fileNames=cellfun( @(x) [tileDir,'/',x], {fileNames.name},'uniformoutput',0);

i=strt;
for i=strt:inc:length(fileNames)
    fprintf('%d: %s\n',i, fileNames{i})
     undoTile(fileNames{i})
end
fprintf('complete/n')