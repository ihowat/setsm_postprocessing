function batchRegisterTiles(tileDir,is2DirOrFile,varargin)
% set paths

%     %   tileDir='/Users/ihowat/project/howat.4/rema_mosaic/rema_mosaic_v13e';
%     is2DirOrFile='/Users/ihowat/data4/REMA/altimetryByTile';

%     %    tileDir='/fs/project/howat.4/howat.4/rema_mosaic_v13e/rema_19_victoria_land';
%     is2DirOrFile='/fs/byo/howat-data4/REMA/altimetryByTile';

strt=1;
inc=1;
resolution='10m';
overwrite=false;
overwrite_flag='';
skipDzfit=false;
dzfitMinPoints=[];

is2Dir = [];
arg_is2TileFile = [];
if exist(is2DirOrFile, 'dir') == 7
    is2Dir = is2DirOrFile;
elseif exist(is2DirOrFile, 'file') == 2
    arg_is2TileFile = is2DirOrFile;
elseif endsWith(is2DirOrFile, '.mat')
    error('is2 tile file does not exist: %s', is2DirOrFile);
else
    error('is2 tile dir does not exist: %s', is2DirOrFile);
end

if length(varargin) == 2 && ~isnan(str2double(varargin{1})) && ~isnan(str2double(varargin{2}))
    strt=varargin{1};
    inc= varargin{2};

elseif ~isempty(varargin)
    n=find(strcmpi(varargin,'start'));
    if ~isempty(n)
        strt=varargin{n+1};
    end

    n=find(strcmpi(varargin,'inc'));
    if ~isempty(n)
        inc=varargin{n+1};
    end

    resolution_choices = {'10m', '2m'};
    n = find(strcmpi('resolution', varargin));
    if ~isempty(n)
        resolution = varargin{n+1};
    end
    if ~ismember(resolution, resolution_choices)
        error("'resolution' must be one of the following, but was '%s': {'%s'}", resolution, strjoin(resolution_choices, "', '"));
    end

    if any(strcmpi(varargin,'overwrite'))
        overwrite=true;
        overwrite_flag='overwrite';
    end

    if any(strcmpi(varargin,'skipDzfit'))
        skipDzfit=true;
    end

    n = find(strcmpi('dzfitMinPoints', varargin));
    if ~isempty(n)
        dzfitMinPoints = varargin{n+1};
    end
end

fprintf('Processing every %d tile, starting at %d\n',inc,strt)

tileFiles = dir(sprintf('%s/*_%s.mat', tileDir, resolution));
%tileFiles = dir([tileDir,'/*.mat']);
tileNames = {tileFiles.name};
tileFiles = fullfile({tileFiles.folder}, tileNames);
tileNames10m = cellfun(@(x) get10mTileName(x), tileNames, 'UniformOutput',false);

tileFiles=tileFiles(:);
tileNames10m=tileNames10m(:);

i=strt;
for i=strt:inc:length(tileFiles)

    % get this tile name and set is2 file name from it
    tileFile=tileFiles{i};
    if endsWith(tileFile, '_reg.mat')
        continue
    end
    tileName10m = tileNames10m{i};
    if ~isempty(arg_is2TileFile)
        is2TileFile = arg_is2TileFile;
    else
        is2TileFile = [is2Dir,'/',tileName10m,'_is2.mat'];
    end
    
    fprintf('Working on tile %d of %d: %s\n',i,length(tileFiles),tileFile);

    regTileFile=strrep(tileFile,'.mat','_reg.mat');
    unregTileFile=strrep(tileFile,'.mat','_unreg.mat.bak');

    if exist(unregTileFile,'file') ~= 2
        fprintf('Creating backup copy of unaltered tile matfile: %s\n',unregTileFile);
        eval(['!cp --preserve=timestamps ',tileFile,' ',unregTileFile]);
    end

    missing_aux_file = false;
    if exist(is2TileFile,'file') ~= 2
        fprintf('icesat registration file does not exist: %s\n',is2TileFile);
        missing_aux_file = true;
    end
    if missing_aux_file
        fprintf('skipping registration and/or dzfit: %s\n',unregTileFile);
        continue;
    end

    if exist(regTileFile,'file') ~= 2 || overwrite
        registerTileToIS2(tileFile,is2TileFile)
        applyRegistration(tileFile,[],overwrite_flag)
    end

    if exist(regTileFile,'file') ~= 2
        fprintf('registered tile was not created, skipping dzfit\n')
    else
        if skipDzfit
            fprintf('skipping dzfit due to provided skipDzfit flag\n')
        else
            fprintf('calculating dzfit\n')
            fit2is2(regTileFile,is2TileFile,'dzfitMinPoints',dzfitMinPoints)
        end
    end
end
