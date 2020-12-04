function initializeMosaic(project,tileName,outDir)

% set to true to remove the subdirectory of subtile .mat files after
% completion:
removeSubtileFlag=true;

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

subTileDir = [outDir,'/',tileName];

outName=[subTileDir,'_10m.mat'];
if exist(outName,'file')
    fprintf('%s exists, skipping\n',outName)
    return
end

% load tile definition file into structure
tileDefs=load(tileDefFile);

% find index of tile in tile def database
tileInd = strcmp(tileDefs.I,tileName);

% crop tileDefs structure to just this tile
tileDefs = structfun( @(x) x(tileInd), tileDefs,'uniformoutput',0);

% get tile boundaries with buffer for making land mask and strip search
% get tile boundaries with buffer for making land mask and strip search
res=10; % ouput mosaic resolution in meters
buffer=100; % size of tile/subtile boundary buffer in meters
minStripOverlap=0.1;
filterFlag= true;

if exist(tileParamListFile,'file')
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

%%  load the strip database file and search for overlapping strips

fprintf('reading database file: %s\n',stripDatabaseFile)
S=shaperead(stripDatabaseFile);

% search for overlapping strips
n = stripSearch({S.X},{S.Y},x0,x1,y0,y1);

% crop structure to this tile
S = S(n);

%% Create meta input file from shp structure

% make full path filanames for each segment without extension
meta.fileName = cellfun(@(w,x,y,z) [stripsDirectory,'/',w,'/strips_v4/2m/',x,'_',y,'/',x,'_seg',num2str(z),'_dem_10m.tif'],{S.region},{S.strip},{S.version},...
    {S.seg_id},'uniformoutput',0);

% convert shapefile vertices to cells in meta struct, removing nans
meta.x = cellfun( @(x) x(~isnan(x)), {S.X}, 'uniformoutput', 0);
meta.y = cellfun( @(x) x(~isnan(x)), {S.Y}, 'uniformoutput', 0);

% make mean scene2strip aligment rmse field
meta.scene_alignment_meanrmse = [S.rmse]';

meta.scene_alignment_meanrmse(meta.scene_alignment_meanrmse <= 0) = NaN;

%% load qc file if it exists and add to meta struct
if exist(qcFile,'file')
    
    fprintf('reading qc file: %s\n',qcFile)
    qc = load(qcFile);
    
    A = cellfun( @(x,y) [x,'_',num2str(y)], qc.stripID, num2cell(qc.seg), 'uniformoutput',0);
    B= cellfun( @(x,y) [x,'_',num2str(y)], strrep({S.strip},'_2m_lsf',''), {S.seg_id}, 'uniformoutput',0);
    
    [~,IA,IB] =  intersect(A,B);
    
    meta.qc.flag = zeros(size(meta.fileName),'uint8');
    meta.qc.x = cell(size(meta.fileName));
    meta.qc.y = cell(size(meta.fileName));
    
    meta.qc.flag(IB) = qc.flag(IA);
    meta.qc.x(IB) = qc.x(IA);
    meta.qc.y(IB) = qc.y(IA);
    
else
    fprintf('no qc file found\n')
end

%% Build and mosaic subtiles
buildSubTiles(tileName,subTileDir,tileDefs,meta,'landTile',landTile,...
    'refDemFile',refDemFile,'buffer',buffer,'minStripOverlap',...
    minStripOverlap,'filter',filterFlag)

call_mosaicSubTiles(subTileDir,10,projection);

if removeSubtileFlag
    rmdir(subTileDir, 's')
end






