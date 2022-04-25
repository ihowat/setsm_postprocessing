function fileNames = selectTilesForMerging(tileDir,varargin)

%tileDir='/Users/ihowat/project/howat.4/rema_mosaic/rema_mosaic_v13e';

if ~isfolder(tileDir)
    error("input 'tileDir' folder does not exist: %s", tileDir);
end

org = 'osu';
resolution = '10m';

if ~isempty(varargin)
    org_choices = {'osu', 'pgc'};
    n = find(strcmpi('org', varargin));
    if ~isempty(n)
        org = varargin{n+1};
    end
    if ~ismember(org, org_choices)
        error("'org' must be one of the following, but was '%s': {'%s'}", org, strjoin(org_choices, "', '"));
    end

    resolution_choices = {'10m', '2m'};
    n = find(strcmpi('resolution', varargin));
    if ~isempty(n)
        resolution = varargin{n+1};
    end
    if ~ismember(resolution, resolution_choices)
        error("'resolution' must be one of the following, but was '%s': {'%s'}", resolution, strjoin(resolution_choices, "', '"));
    end
end

if strcmp(org, 'osu')
    unregFiles=dir([tileDir,'/*_',resolution,'.mat']);
else
    unregFiles=dir([tileDir,'/*/*_',resolution,'.mat']);
end
unregFiles=fullfile({unregFiles.folder}, {unregFiles.name});

if ~isempty(unregFiles)
    
    regFiles = strrep(unregFiles, [resolution,'.mat'], [resolution,'_reg.mat']);
    
    n = cellfun(@exist, regFiles);
    n = n==2;
    
    fileNames = [regFiles(n) unregFiles(~n)];
    
else
    if strcmp(org, 'osu')
        fileNames=dir([tileDir,'/*_',resolution,'_reg.mat']);
    else
        fileNames=dir([tileDir,'/*/*_',resolution,'_reg.mat']);
    end
    fileNames=fullfile({fileNames.folder}, {fileNames.name});
end
