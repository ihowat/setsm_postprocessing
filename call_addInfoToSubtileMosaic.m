function call_addInfoToSubtileMosaic(tileName,dx)
%tileName='48_13';
%dx=2;

if ismac
    subTileDir=['/Users/ihowat/project/howat.4/earthdem/earthdem_mosaic_testing_1km/',tileName];
    addpath('/Users/ihowat/unity-home/demtools')
else
    
    subTileDir=['/fs/project/howat.4/howat.4/earthdem/earthdem_mosaic_testing_1km/',tileName];
    addpath('/home/howat.4/demtools')
end

if dx == 10
    outName=[subTileDir,'_10m.mat'];
    
    addInfoToSubtileMosaic(subTileDir,10,outName);
    
else
    
    outName=[subTileDir,'_1_1_2m.mat'];
    addInfoToSubtileMosaic(subTileDir,2,outName,'1_1');
    
    outName=[subTileDir,'_1_2_2m.mat'];
    addInfoToSubtileMosaic(subTileDir,2,outName,'1_2');
    
    outName=[subTileDir,'_2_1_2m.mat'];
    addInfoToSubtileMosaic(subTileDir,2,outName,'2_1');
    
    outName=[subTileDir,'_2_2_2m.mat'];
    addInfoToSubtileMosaic(subTileDir,2,outName,'2_2');
    
end

