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
    regFiles=dir([tileDir,'/*_',resolution,'_reg.mat']);
else
    unregFiles=dir([tileDir,'/*/*_',resolution,'.mat']);
    regFiles=dir([tileDir,'/*/*_',resolution,'_reg.mat']);
end

regFiles = fullfile({regFiles.folder}, {regFiles.name});
unregFiles = fullfile({unregFiles.folder}, {unregFiles.name});
unregFiles = setdiff(unregFiles, strrep(regFiles, '_reg.mat', '.mat'));
fileNames = [unregFiles regFiles];
