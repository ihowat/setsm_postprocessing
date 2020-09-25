function call_mosaicSubTiles(subTileDir,res,projection)
%tileName='48_13';
%res=2;

if res == 10
    outName=[subTileDir,'_10m.mat'];

	if exist(outName,'file')
		fprintf('%s exists, skipping\n',outName)
		return
	end
    
    mosaicSubTiles(subTileDir,10,outName,'projection',projection);
    
else
    
    outName=[subTileDir,'_1_1_2m.mat'];
    mosaicSubTiles(subTileDir,2,outName,'1_1','projection',projection);
    
    outName=[subTileDir,'_1_2_2m.mat'];
    mosaicSubTiles(subTileDir,2,outName,'1_2','projection',projection);
    
    outName=[subTileDir,'_2_1_2m.mat'];
    mosaicSubTiles(subTileDir,2,outName,'2_1','projection',projection);
    
    outName=[subTileDir,'_2_2_2m.mat'];
    mosaicSubTiles(subTileDir,2,outName,'2_2','projection',projection);
    
end

