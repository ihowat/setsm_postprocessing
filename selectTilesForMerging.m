function fileNames= selectTilesForMerging(tileDir)

%tileDir='/Users/ihowat/project/howat.4/rema_mosaic/rema_mosaic_v13e';

unregFiles=dir([tileDir,'/*_10m.mat']);
%unregFiles=dir([tileDir,'/*/*_10m.mat']);
unregFiles=fullfile({unregFiles.folder}, {unregFiles.name});

if ~isempty(unregFiles)
    
    regFiles = strrep(unregFiles,'10m.mat','10m_reg.mat');
    
    n = cellfun( @exist, regFiles);
    n = n==2;
    
    fileNames = [regFiles(n) unregFiles(~n)];
    
else
    fileNames=dir([tileDir,'/*_10m_reg.mat']);
%    fileNames=dir([tileDir,'/*/*_10m_reg.mat']);
    fileNames=fullfile({fileNames.folder}, {fileNames.name});
end

