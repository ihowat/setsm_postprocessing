function batchRegisterTilesToCOP30(tileDir,refDemDirOrFile,waterMaskDirOrFile,varargin)
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

refDemDir = [];
arg_refDemTif = [];
if exist(refDemDirOrFile, 'dir') == 7
    refDemDir = refDemDirOrFile;
elseif exist(refDemDirOrFile, 'file') == 2
    arg_refDemTif = refDemDirOrFile;
elseif endsWith(refDemDirOrFile, '.mat')
    error('refDemTif tile file does not exist: %s', refDemDirOrFile);
else
    error('refDemTif tile dir does not exist: %s', refDemDirOrFile);
end

waterMaskDir = [];
arg_waterMaskTif = [];
if exist(waterMaskDirOrFile, 'dir') == 7
    waterMaskDir = waterMaskDirOrFile;
elseif exist(waterMaskDirOrFile, 'file') == 2
    arg_waterMaskTif = waterMaskDirOrFile;
elseif endsWith(waterMaskDirOrFile, '.mat')
    error('waterMask tile file does not exist: %s', waterMaskDirOrFile);
else
    error('waterMask tile dir does not exist: %s', waterMaskDirOrFile);
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

    % get this tile name and set aux tile file names from it
    tileFile=tileFiles{i};
    if endsWith(tileFile, '_reg.mat')
        continue
    end
    tileName10m = tileNames10m{i};
    if ~isempty(arg_refDemTif)
        refDemTif = arg_refDemTif;
    else
        refDemTif = [refDemDir,'/',tileName10m,'_10m_cop30_wgs84.tif'];
    end
    if ~isempty(arg_waterMaskTif)
        waterMaskTif = arg_waterMaskTif;
    else
        waterMaskTif = [waterMaskDir,'/',tileName10m,'_10m_cover.tif'];
    end
    
    fprintf('Working on tile %d of %d: %s\n',i,length(tileFiles),tileFile);

    unregTileFile=strrep(tileFile,'.mat','_unreg.mat.bak');
    regTempFile=strrep(tileFile,'.mat','_reg.mat.temp');
    regTileFile=strrep(tileFile,'.mat','_reg.mat');

    if exist(unregTileFile,'file') ~= 2
        fprintf('Creating backup copy of unaltered tile matfile: %s\n',unregTileFile);
        eval(['!cp --preserve=timestamps ',tileFile,' ',unregTileFile]);
    end
    
    if exist(regTileFile,'file') == 2
        if overwrite
            fprintf('reg matfile already exists, will overwrite: %s\n',regTileFile);
        else
            fprintf('reg matfile already exists, skipping: %s\n',regTileFile);
            continue;
        end
    end

    missing_aux_file = false;
    if exist(refDemTif,'file') ~= 2
        fprintf('reference DEM file does not exist: %s\n',refDemTif);
        missing_aux_file = true;
    end
    if exist(waterMaskTif,'file') ~=2
        fprintf('watermask file does not exist: %s\n',waterMaskTif);
        missing_aux_file = true;
    end
    if missing_aux_file
        fprintf('skipping registration: %s\n',unregTileFile);
        continue;
    end

    if exist(regTileFile,'file') ~= 2 || overwrite
        try
            if exist(regTileFile, 'file') == 2
                fprintf('Removing existing reg.mat file: %s\n', regTileFile);
                delete(regTileFile);
            end
            if exist(regTempFile, 'file') == 2
                fprintf('Removing existing reg.mat.tmp file: %s\n', regTempFile);
                delete(regTempFile);
            end

            fprintf('Calculating vertical registration\n');
            [z_reg,~,~,~] = registerTileToCOP30(unregTileFile,[],[],refDemTif,waterMaskTif,'registerBlobs');

            fprintf('Copying unreg.mat file to new reg.mat.tmp file: %s\n',regTempFile);
            eval(['!cp ',unregTileFile,' ',regTempFile]);

            fprintf('Writing registered z to reg.mat.tmp file\n');
            m1=matfile(regTempFile);
            m1.Properties.Writable = true;
            m1.z = z_reg;
            m1.regToCOP30 = true;
            m1.aux_regCOP30_refDemTif = refDemTif;

            fprintf('Renaming reg.mat.tmp file to reg.mat: %s\n',regTileFile);
            eval(['!mv ',regTempFile,' ',regTileFile]);
        catch ME
            if exist(regTileFile, 'file') == 2
                fprintf('Detected error, removing reg.mat file: %s\n',regTileFile);
                delete(regTileFile);
            end
            if exist(regTempFile, 'file') == 2
                fprintf('Detected error, removing reg.mat.tmp file: %s\n',regTempFile);
                delete(regTempFile);
            end
            rethrow(ME);
        end
    end

end
