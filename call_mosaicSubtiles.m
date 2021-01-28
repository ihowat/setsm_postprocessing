function call_mosaicSubtiles(tileName,dx,projection,tileDefFile)
%tileName='48_13';
%dx=2;

if ismac
    subTileDir=['/Users/ihowat/project/howat.4/earthdem/earthdem_mosaic_testing_1km/',tileName];
    addpath('/Users/ihowat/unity-home/demtools')
%    tileDefFile = '/Users/ihowat/unity-home/earthdem/PGC_Imagery_Mosaic_Tiles_Arctic.mat';
else
    
    subTileDir=['/fs/project/howat.4/howat.4/earthdem/earthdem_mosaic_testing_1km/',tileName];
    addpath('/home/howat.4/demtools')
%    tileDefFile = '~/earthdem/PGC_Imagery_Mosaic_Tiles_Arctic.mat'; %PGC/NGA Tile definition file
end

if dx == 10
    outName=[subTileDir,'_10m.mat'];
    
    [x0,x1,y0,y1]=getTileExtents(tileName,tileDefFile);
    
    mosaicSubTiles(subTileDir,10,outName,'extent',[x0,x1,y0,y1]);
    
else
    
    outName=[subTileDir,'_1_1_2m.mat'];
    [x0,x1,y0,y1]=getTileExtents(tileName,tileDefFile,'quadrant','1_1');
    mosaicSubTiles(subTileDir,2,outName,'projection',projection,'quadrant','1_1','extent',[x0,x1,y0,y1]);
    
    outName=[subTileDir,'_1_2_2m.mat'];
    mosaicSubTiles(subTileDir,2,outName,'projection',projection,'quadrant','1_2','extent',[x0,x1,y0,y1]);
    
    outName=[subTileDir,'_2_1_2m.mat'];
    mosaicSubTiles(subTileDir,2,outName,'projection',projection,'quadrant','2_1','extent',[x0,x1,y0,y1]);
    
    outName=[subTileDir,'_2_2_2m.mat'];
    mosaicSubTiles(subTileDir,2,outName,'projection',projection,'quadrant','2_2','extent',[x0,x1,y0,y1]);
    
end

