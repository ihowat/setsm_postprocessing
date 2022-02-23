function fileNames= selectTilesForMerging(tileDir,varargin)

%tileDir='/Users/ihowat/project/howat.4/rema_mosaic/rema_mosaic_v13e';

n = find(strcmpi('resolution',varargin));
if ~isempty(n)
    resolution = varargin{n+1};
else
    resolution = '10m';
end

unregFiles=dir([tileDir,'/*_',resolution,'.mat']);
%unregFiles=dir([tileDir,'/*/*_',resolution,'.mat']);
unregFiles=fullfile({unregFiles.folder}, {unregFiles.name});

if ~isempty(unregFiles)
    
    regFiles = strrep(unregFiles,'10m.mat','10m_reg.mat');
    
    n = cellfun( @exist, regFiles);
    n = n==2;
    
    fileNames = [regFiles(n) unregFiles(~n)];
    
else
    fileNames=dir([tileDir,'/*_',resolution,'_reg.mat']);
%    fileNames=dir([tileDir,'/*/*_',resolution,'_reg.mat']);
    fileNames=fullfile({fileNames.folder}, {fileNames.name});
end

