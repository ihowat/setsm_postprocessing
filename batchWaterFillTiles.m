function batchWaterFillTiles(tileDir,refDemDirOrFile,waterMaskDirOrFile,varargin)
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

regFiles = dir(sprintf('%s/*_%s_reg.mat', tileDir, resolution));
unregFiles = dir(sprintf('%s/*_%s.mat', tileDir, resolution));
regFiles = fullfile({regFiles.folder}, {regFiles.name});
unregFiles = fullfile({unregFiles.folder}, {unregFiles.name});
unregFiles = setdiff(unregFiles, strrep(regFiles, '_reg.mat', '.mat'));
tileFiles = [unregFiles regFiles];

tileNames10m = cellfun(@(x) get10mTileName(x), tileFiles, 'UniformOutput',false);

tileFiles=tileFiles(:);
tileNames10m=tileNames10m(:);

i=strt;
for i=strt:inc:length(tileFiles)

    % get this tile name and set aux tile file names from it
    tileFile=tileFiles{i};
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

    unfillTileFile=strrep(tileFile,'.mat','_unfill.mat.bak');
    fillTempFile=strrep(tileFile,'.mat','_fill.mat.temp');
    fillTileFile=strrep(tileFile,'.mat','_fill.mat');

    if exist(unfillTileFile,'file') ~= 2
        fprintf('Creating backup copy of unaltered tile matfile: %s\n',unfillTileFile);
        eval(['!cp --preserve=timestamps ',tileFile,' ',unfillTileFile]);
    end
    
    if exist(fillTileFile,'file') == 2
        if overwrite
            fprintf('fill matfile already exists, will overwrite: %s\n',fillTileFile);
        else
            fprintf('fill matfile already exists, skipping: %s\n',fillTileFile);
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
        fprintf('skipping fill: %s\n',unfillTileFile);
        continue;
    end

    if exist(fillTileFile,'file') ~= 2 || overwrite
        try
            if exist(fillTileFile, 'file') == 2
                fprintf('Removing existing fill.mat file: %s\n', fillTileFile);
                delete(fillTileFile);
            end
            if exist(fillTempFile, 'file') == 2
                fprintf('Removing existing fill.mat.tmp file: %s\n', fillTempFile);
                delete(fillTempFile);
            end

            fprintf('Calculating waterfill\n');
            [z_fill,waterfill_mask,~,~,~] = fillWater(tileFile,waterMaskTif,refDemTif);

            fprintf('Copying unfill.mat file to new fill.mat.tmp file: %s\n',fillTempFile);
            eval(['!cp ',unfillTileFile,' ',fillTempFile]);

            fprintf('Writing filled z to fill.mat.tmp file\n');
            m1=matfile(fillTempFile);
            m1.Properties.Writable = true;
            m1.z = z_fill;
            m1.waterFillMask = waterfill_mask;
            m1.waterFilled = true;

            fprintf('Renaming fill.mat.tmp file to fill.mat: %s\n',fillTileFile);
            eval(['!mv ',fillTempFile,' ',fillTileFile]);
        catch ME
            if exist(fillTileFile, 'file') == 2
                fprintf('Detected error, removing fill.mat file: %s\n',fillTileFile);
                delete(fillTileFile);
            end
            if exist(fillTempFile, 'file') == 2
                fprintf('Detected error, removing fill.mat.tmp file: %s\n',fillTempFile);
                delete(fillTempFile);
            end
            rethrow(ME);
        end
    end

end
