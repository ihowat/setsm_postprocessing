function [meta,landTile,buffer,minStripOverlap,filterFlag] = initializeMosaic(project,tileName,outDir,varargin)

% set to true to remove the subdirectory of subtile .mat files after
% completion:
removeSubtileFlag=false;

subTileDir='';
res=10; % ouput mosaic resolution in meters
projection='';
tileDefFile='';
stripDatabaseFile='';
stripsDirectory='';
waterTileDir='';
refDemFile='';
tileqcDir='';
qcFile='';
tileParamListFile='';
returnMetaOnly=false;
dateFiltStart='';
dateFiltEnd='';

% Set Paths and filenames
if strcmpi(project,'rema')
    
    projection='polar stereo south';
    
    if ismac
        tileDefFile = '/Users/ihowat/gdrive/projects/rema2/rema_tile_definitions.mat';
        stripDatabaseFile = '/Users/ihowat/rema_strips.shp';
        stripsDirectory='/Users/ihowat/project/REMA';
        refDemFile='/Users/ihowat/project/REMA/TanDEM-X_Antarctica_90m/TanDEM_Antarctica_Mosaic.tif';
       %subTileDir = ['/Users/ihowat/project/howat.4/rema_mosaic/',tileName];
        %qcFile = '/Users/ihowat/rema_strip_v4_qc.mat';
        qcFile = '';
    else
        tileDefFile = '/fs/byo/howat-data4/REMA/rema_tile_definitions.mat';
        stripDatabaseFile = '/fs/byo/howat-data4/REMA/metadata_rema_active/rema_strips.shp';
        stripsDirectory='/fs/byo/howat-data4/REMA';
        refDemFile='/fs/byo/howat-data4/REMA/TanDEM-X_Antarctica_90m/TanDEM_Antarctica_Mosaic.tif';
        qcFile = ['/fs/byo/howat-data4/REMA/rema_strip_automosaic_qc/rema_strip_automosaic_qc_',tileName,'.mat'];
        tileParamListFile='/fs/byo/howat-data4/REMA/tileParamList.txt';
    end
    
elseif strcmpi(project,'arcticdem')
    
    projection='polar stereo north';
    
    if ismac
        tileDefFile = '/Users/ihowat/unity-home/earthdem/PGC_Imagery_Mosaic_Tiles_Arctic.mat';
        stripDatabaseFile = '/Users/ihowat/gdrive/projects/earthdem/earthdem_database_unf.mat';
        waterTileDir='~/project/howat.4/earthdem/masks';
        %  refDemFile='~/TDM90_Greenland/TDM1_DEM__30_mosaic.tif';
    else
        tileDefFile = '~/earthdem/PGC_Imagery_Mosaic_Tiles_Arctic.mat'; %PGC/NGA Tile definition file
        stripDatabaseFile = '/fs/byo/howat-data5/ArcticDEM/metadata_arcticdem/arcticdem_strips.shp';
        stripsDirectory='/fs/byo/howat-data5/ArcticDEM';
        waterTileDir='/fs/project/howat.4/howat.4/earthdem/masks';
        %refDemFile='/fs/project/howat.4/EarthDEM/tandemx_alaska_mosaic_3413_tap90m.tif';
        refDemFile='/fs/project/howat.4/EarthDEM/TDM1_DEM__30_mosaic.tif';
      %  subTileDir = ['/fs/project/howat.4/howat.4/arcticdem_mosaic/',tileName];
    end
    
end


fprintf('project = %s\n',project)
fprintf('tileName = %s\n',tileName)
fprintf('outDir = %s\n',outDir)

if ~isempty(project)

    subTileDir = [outDir,'/',tileName];

    outName=[subTileDir,'_10m.mat'];
    if exist(outName,'file') == 2
        fprintf('%s exists, skipping\n',outName)
        return
    end

else
    fprintf('subTileDir = %s\n',subTileDir)
end


% Parse varargins
n = find(strcmpi(varargin,'res'));
if ~isempty(n)
    res = varargin{n+1};
end
fprintf('res = %d\n',res)

n = find(strcmpi(varargin,'projection'));
if ~isempty(n)
    projection = varargin{n+1};
end
fprintf('projection = %s\n',projection)

n = find(strcmpi(varargin,'subTileDir'));
if ~isempty(n)
    subTileDir = varargin{n+1};
end
fprintf('subTileDir = %s\n',subTileDir)

n = find(strcmpi(varargin,'tileDefFile'));
if ~isempty(n)
    tileDefFile = varargin{n+1};
end
fprintf('tileDefFile = %s\n',tileDefFile)

n = find(strcmpi(varargin,'stripDatabaseFile'));
if ~isempty(n)
    stripDatabaseFile = varargin{n+1};
end
fprintf('stripDatabaseFile = %s\n',stripDatabaseFile)

n = find(strcmpi(varargin,'stripsDirectory'));
if ~isempty(n)
    stripsDirectory = varargin{n+1};
end
fprintf('stripsDirectory = %s\n',stripsDirectory)

n = find(strcmpi(varargin,'waterTileDir'));
if ~isempty(n)
    waterTileDir = varargin{n+1};
end
fprintf('waterTileDir = %s\n',waterTileDir)

n = find(strcmpi(varargin,'refDemFile'));
if ~isempty(n)
    refDemFile = varargin{n+1};
end
fprintf('refDemFile = %s\n',refDemFile)

n = find(strcmpi(varargin,'tileqcDir'));
if ~isempty(n)
    tileqcDir = varargin{n+1};
end
fprintf('tileqcDir = %s\n',tileqcDir)

n = find(strcmpi(varargin,'tileParamListFile'));
if ~isempty(n)
    tileParamListFile = varargin{n+1};
end
fprintf('tileParamListFile = %s\n',tileParamListFile)

n = find(strcmpi(varargin,'returnMetaOnly'));
if ~isempty(n)
    returnMetaOnly = varargin{n+1};
end
fprintf('returnMetaOnly = %d\n',returnMetaOnly)

n = find(strcmpi(varargin,'dateFiltStart'));
if ~isempty(n)
    dateFiltStart = varargin{n+1};
end
fprintf('dateFiltStart = %s\n',dateFiltStart)

n = find(strcmpi(varargin,'dateFiltEnd'));
if ~isempty(n)
    dateFiltEnd = varargin{n+1};
end
fprintf('dateFiltEnd = %s\n',dateFiltEnd)


% attempt to set uninitialized qc file location
if isempty(qcFile) && ~isempty(tileqcDir)
    qcFile = [tileqcDir,'/',tileName,'.mat'];
    fprintf('qcFile = %s\n',qcFile)
end


% load tile definition file into structure
tileDefs=load(tileDefFile);

if startsWith(tileName,'utm')
    sl = split(tileName,'_');
    tilePrefix = [sl{1},'_'];
    tileName_in_tileDef = strjoin(sl(2:3),'_');
else
    tilePrefix = '';
    tileName_in_tileDef = tileName;
end

% if projection argument not given, check tile definition
tileDefs_projection = '';
if isfield(tileDefs, 'projection')
    tileDefs_projection = tileDefs.projection;
    tileDefs = rmfield(tileDefs, 'projection');
elseif isfield(tileDefs, 'projstr')
    tileDefs_projection = tileDefs.projstr;
    tileDefs = rmfield(tileDefs, 'projstr');
end
if isempty(projection) && ~isempty(tileDefs_projection)
    projection = tileDefs_projection;
end

% find index of tile in tile def database
tileInd = strcmp(tileDefs.I,tileName_in_tileDef);

% crop tileDefs structure to just this tile
tileDefs = structfun( @(x) x(tileInd), tileDefs,'uniformoutput',0);

% get tile boundaries with buffer for making land mask and strip search
%res=10; % ouput mosaic resolution in meters
buffer=100; % size of tile/subtile boundary buffer in meters
minStripOverlap=0.1;
filterFlag= true;

if ~isempty(tileParamListFile)
    if exist(tileParamListFile,'file') == 2
        [~,tileParams]=system(['grep ',tileName,' ',tileParamListFile]);
        paramReadFailFlag=true;
        if ~isempty(isempty(tileParams))
            tileParams=strtrim(tileParams);
            tileParams=strrep(tileParams,'_',' ');
            tileParams=str2num(tileParams);

            if length(tileParams) == 4
                buffer=tileParams(3);
                minStripOverlap=tileParams(4);
                paramReadFailFlag=false;
            end
        end

        if  paramReadFailFlag == true
            fprintf('Could not find %s in %s, using defaults\n',...
                tileName,tileParamListFile)
        end
    else
        fprintf('tileParamListFile %s does not exist, using defaults\n', tileParamListFile)
    end
end

x0=tileDefs.x0-buffer;
y0=tileDefs.y0-buffer;
x1=tileDefs.x1+buffer;
y1=tileDefs.y1+buffer;

% make land/water(true/false) mask
if isfield(tileDefs,'coastlinePolyshape')
    [landTile.z,landTile.x,landTile.y] = polyshape2maskmap(...
        tileDefs.coastlinePolyshape,[x0 x1 y0 y1],res);
else
    landTile = waterTileDir;
end

% convert argument/tileDefs projection to abbreviated tile projection "tileProjName"
tileProjName = '';
if startsWith(tileName,'utm')
    tileName_parts = strsplit(tileName,'_');
    tileProjName = tileName_parts{1};
end
if ~isempty(projection)

    if strcmpi(projection, 'polar stereo north')
        tileProjName = 'psn';
    elseif strcmpi(projection, 'polar stereo south')
        tileProjName = 'pss';
    elseif startsWith(projection,'utm', 'IgnoreCase',true)
        if ~startsWith(tileProjName, projection, 'IgnoreCase',true)
            error("'projection' argument (%s) and tileName (%s) prefix (%s) do not match", ...
                projection, tileName, tileProjName);
        end
    else
        error("No matching tileProjName for 'projection' argument: '%s'", tileProjName);
    end
end

%%  load the strip database file and search for overlapping strips
if endsWith(stripDatabaseFile,'.shp', 'IgnoreCase',true)

    fprintf('reading strip database shapefile: %s\n',stripDatabaseFile)
    S=shaperead(stripDatabaseFile);

    % (the following isn't necessary because a strip database in shapefile format
    %  will only contain strip records that are in the same projection)
%    if isfield(S,'projection')
%        % trim strip database to only strips with a projection matching the tile
%        if isempty(projection)
%            error("'projection' argument must be provided when strip dbase has 'projection' field");
%        end
%        n = strcmpi({S.projection},tileProjName);
%        if ~any(n)
%            error("no matches between tileProjName='%s' and strip dbase 'projection' field", tileProjName);
%        end
%        S = S(n);
%    end

    % search for overlapping strips
    n = stripSearch({S.X},{S.Y},x0,x1,y0,y1);
    if isempty(n)
        error("no overlapping strips")
    end

    % crop structure to this tile
    S = S(n);

    %% Create meta input file from shp structure

    % make full path filanames for each segment without extension
    meta.fileName = cellfun(@(w,x,y,z) [stripsDirectory,'/',w,'/strips_v4.1/2m/',regexprep(x,'^SETSM_s2s\d{3}_',''),'_',y,'/',x,'_seg',num2str(z),'_dem_10m.tif'],{S.region},{S.strip},{S.version},...
        {S.seg_id},'uniformoutput',0);

    meta.strip = {S.strip};

    % ensure strip filename scheme is as expected
    if isempty(regexp(meta.fileName{1}, 'SETSM_s2s041'))
        error("Expected meta.fileName items to contain 'SETSM_s2s041' substring, but it is not present: %s", meta.fileName{1});
    end

    % convert shapefile vertices to cells in meta struct, removing nans
    meta.x = cellfun( @(x) x(~isnan(x)), {S.X}, 'uniformoutput', 0);
    meta.y = cellfun( @(x) x(~isnan(x)), {S.Y}, 'uniformoutput', 0);

    % make mean scene2strip aligment rmse field
    meta.scene_alignment_meanrmse = [S.rmse]';

    meta.scene_alignment_meanrmse(meta.scene_alignment_meanrmse <= 0) = NaN;


elseif endsWith(stripDatabaseFile,'.mat', 'IgnoreCase',true)

    fprintf('reading strip database .mat file: %s\n',stripDatabaseFile)
    meta=load(stripDatabaseFile);

    % ensure strip filename scheme is as expected
    if isempty(regexp(meta.fileName{1}, 'SETSM_s2s041'))
        error("Expected meta.fileName items to contain 'SETSM_s2s041' substring, but it is not present: %s", meta.fileName{1});
    end

    % trim strip database to only strips with a projection matching the tile
    if isfield(meta,'strip_projection_name')
        if isempty(projection)
            error("'projection' argument must be provided when strip dbase has 'strip_projection_name' field");
        end
        n = strcmpi(meta.strip_projection_name,tileProjName);
        if ~any(n)
            error("no matches between tileProjName='%s' and strip dbase 'strip_projection_name' field", tileProjName);
        end
        meta = structfun(@(x) x(n), meta, 'uniformoutput',0);
    end

    % search for overlapping strips
    n = stripSearch(meta.x,meta.y,x0,x1,y0,y1);
    if isempty(n)
        error("No strips overlap tile: %s", tileName)
    end
    meta = structfun(@(x) x(n), meta, 'uniformoutput',0);

    % make mean scene2strip aligment rmse field
    if ~isfield(meta,'scene_alignment_meanrmse')
        if isfield(meta,'avg_rmse')
            % this is from an old version of the meta files that used 0 in mean
            error('avg_rmse field in meta structure, needs to be updated')
        elseif isfield(meta,'scene_alignment')
            % need to rm zeros (first scene) and nans (unused redundant scenes),
            % strips w/ 1 scene will be NaN
            meta.scene_alignment_meanrmse = cellfun(@(x)...
                mean(x.rmse(x.rmse~=0 & ~isnan(x.rmse))), meta.scene_alignment);
        else
            error('missing scene alignment field in meta structure')
        end
    end

end

if ~isfield(meta,'A')
    % get strip areas and alignment stats for quality selection
    meta.A = cellfun(@(x,y) polyarea(x,y), meta.x,meta.y);
end


if ~isempty(dateFiltStart) || ~isempty(dateFiltEnd)
    fprintf('filtering strip overlap to records that match date filter\n')
    fprintf('dateFiltStart = %s\n',dateFiltStart)
    fprintf('dateFiltEnd = %s\n',dateFiltEnd)

    if ~isfield(meta,'stripDate')
        meta.stripDate=cellfun(@(x) datenum(parsePairnameDatestring(x),'yyyymmdd'), meta.strip, 'uniformoutput',0);
        meta.stripDate=cell2mat(meta.stripDate);
    end

    n = [];
    if ~isempty(dateFiltStart)
        n_dateFiltStart = meta.stripDate >= datenum(dateFiltStart,'yyyymmdd');
    end
    if ~isempty(dateFiltEnd)
        n_dateFiltEnd = meta.stripDate <= datenum(dateFiltEnd,'yyyymmdd');
    end
    if ~isempty(dateFiltStart) && ~isempty(dateFiltEnd)
        n = n_dateFiltStart & n_dateFiltEnd;
    elseif ~isempty(dateFiltStart)
        n = n_dateFiltStart;
    elseif ~isempty(dateFiltEnd)
        n = n_dateFiltEnd;
    end

    fprintf('%d of %d overlapping strips match date filter, applying filter\n', nnz(n), length(meta.stripDate))
    meta = structfun(@(x) x(n), meta, 'uniformoutput',0);
end


%% load qc file if it exists and add to meta struct
if exist(qcFile,'file') == 2
    
    fprintf('reading qc file: %s\n',qcFile)
    qc = load(qcFile);

    % ensure qc stripID scheme is as expected
    if isempty(regexp(qc.stripID{1}, 'SETSM_s2s041'))
        error("Expected qc.stripID items to contain 'SETSM_s2s041' substring, but it is not present: %s", qc.stripID{1});
    end
    
    A = cellfun( @(x,y) [x,'_',num2str(y)], qc.stripID, num2cell(qc.seg), 'uniformoutput',0);
    if exist('S','var')
        B= cellfun( @(x,y) [x,'_',num2str(y)], {S.strip}, {S.seg_id}, 'uniformoutput',0);
    else
        B= cellfun( @(x) stripFileNameToSegID(x), meta.fileName, 'uniformoutput',0);
    end
    
    [~,IA,IB] =  intersect(A,B);
    if isempty(IA)
        error("Could not match qc.stripID entries to any meta.fileName records");
    end
    
    meta.qc.flag = zeros(size(meta.fileName),'uint8');
    meta.qc.x = cell(size(meta.fileName));
    meta.qc.y = cell(size(meta.fileName));
    
    meta.qc.flag(IB) = qc.flag(IA);
    meta.qc.x(IB) = qc.x(IA);
    meta.qc.y(IB) = qc.y(IA);
    
else
    fprintf('no qc file found\n')
end


if returnMetaOnly; return; end


%% Build and mosaic subtiles
buildSubTiles(tileName,subTileDir,tileDefs,meta,'landTile',landTile,...
    'refDemFile',refDemFile,'buffer',buffer,'minStripOverlap',...
    minStripOverlap,'filter',filterFlag,'projection',projection)

call_mosaicSubTiles(subTileDir,10,projection,tileDefFile);

if removeSubtileFlag
    rmdir(subTileDir, 's')
end


function strip_seg_id = stripFileNameToSegID(fileName)
pname_seg_re = '^.+/(?<pairname>SETSM_s2s041_[A-Z0-9]{4}_[0-9]{8}_[A-F0-9]{16}_[A-F0-9]{16})_[^/]+/[^/]+_seg(?<seg_id>\d+)_[^/]+$';
[tokens, match_idx] = regexp(fileName, pname_seg_re, 'names');
if isempty(match_idx)
    error("Cannot parse strip fileName parts with regex '%s' from input fileName string: %s", tilename_re, tilename);
end
strip_seg_id = [tokens.pairname,'_',tokens.seg_id];




