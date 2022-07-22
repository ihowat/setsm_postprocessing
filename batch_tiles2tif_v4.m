function batch_tiles2tif_v4(tileDir,projstr,varargin)

%% Parse optional args

% Arg default values
resolution = '10m';
metaOnly = false;
noMeta = false;

if ~isempty(varargin)
    n = find(strcmpi('resolution', varargin));
    if ~isempty(n)
        resolution = varargin{n+1};
    end

    n = find(strcmpi('metaOnly', varargin));
    if ~isempty(n)
        metaOnly = true;
    end

    n = find(strcmpi('noMeta', varargin));
    if ~isempty(n)
        noMeta = true;
    end
end

%% Verify args

if ~isfolder(tileDir)
    error("input 'tileDir' folder does not exist: %s", tileDir);
end

resolution_choices = {'10m', '2m'};
if ~ismember(resolution, resolution_choices)
    error("'resolution' must be one of the following, but was '%s': {'%s'}", resolution, strjoin(resolution_choices, "', '"));
end


%% Gather input tiles to process

unregFiles=dir([tileDir,'/*_',resolution,'.mat']);
unregFiles=fullfile({unregFiles.folder}, {unregFiles.name});

if ~isempty(unregFiles)
    regFiles = strrep(unregFiles, [resolution,'.mat'], [resolution,'_reg.mat']);

    n = cellfun(@exist, regFiles);
    n = n==2;

    fileNames = [regFiles(n) unregFiles(~n)];
else
    fileNames=dir([tileDir,'/*_',resolution,'_reg.mat']);
    fileNames=fullfile({fileNames.folder}, {fileNames.name});
end


%% Main processing loop

num_tiles = length(fileNames);
for tile_idx = 1:num_tiles
    tile_file = fileNames{tile_idx};
    fprintf("\nProcessing tile (%d/%d): %s\n\n", tile_idx, num_tiles, tile_file);
    if ~metaOnly
        writeTileToTifv4(tile_file, projstr, varargin{:});
    end
    if ~noMeta
        tileMetav4(tile_file, varargin{:});
    end
end
