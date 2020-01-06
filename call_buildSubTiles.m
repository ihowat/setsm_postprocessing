function call_buildSubTiles(tileName)
%paths/files
if ismac
    tileDefFile = '/Users/ihowat/unity-home/earthdem/PGC_Imagery_Mosaic_Tiles_Arctic.mat';
    databaseFile = '/Users/ihowat/gdrive/projects/earthdem/earthdem_database_unf.mat';
    outDir = ['/Users/ihowat/project/howat.4/earthdem/earthdem_mosaic_testing_1km/',tileName];
    addpath('/Users/ihowat/unity-home/demtools');
    waterTileDir='~/data/pgc_projects/ak_water_rasters_v2';
    refDemFile='~/tandemx_alaska_mosaic_3413_tap90m.tif';
    % coastlinePolyFile='/Users/ihowat/gdrive/projects/earthdem/gshhg_237_alaska_coastline_3413.mat';
    % lakePolyFile='/Users/ihowat/gdrive/projects/earthdem/gshhg_237_alaska_lakes_3413.mat';
else
    tileDefFile = '~/earthdem/PGC_Imagery_Mosaic_Tiles_Arctic.mat'; %PGC/NGA Tile definition file
    databaseFile = '~/earthdem/earthdem_database_unf.mat';
    outDir = ['/home/howat.4/project/earthdem/earthdem_mosaic_testing_1km/',tileName];
    addpath('/home/howat.4/demtools');
    waterTileDir='/fs/byo/howat-data/pgc_projects/ak_water_rasters_v2';
    refDemFile='/fs/project/howat.4/EarthDEM/tandemx_alaska_mosaic_3413_tap90m.tif';
    %     coastlinePolyFile='gshhg_237_alaska_coastline_3413.mat';
    %     lakePolyFile='gshhg_237_alaska_lakes_3413.mat';
end


buildSubTiles(tileName,outDir,tileDefFile,databaseFile,waterTileDir,refDemFile)