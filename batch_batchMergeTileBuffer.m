function batch_batchMergeTileBuffer(outdir,tileList,res)

%% Blend tile edges
% get list of tiles recursively where the parent dir matches an item in the
% cell attay tileList
% Usage: batch_batchMergeTileBuffer('/path/to/tile/dirs',{'14_51_1_1','14_51_1_2'},'2m')
% Usage: batch_batchMergeTileBuffer('/path/to/tile/dirs',{'14_51','14_52'},'10m')

tilef={};

for i=1:length(tileList)
    tile=tileList{i};

    tile_parts=split(tile,'_');
    if startsWith(tile,'utm')
        supertile = strjoin(tile_parts(1:3),'_');
    else
        supertile = strjoin(tile_parts(1:2),'_');
    end

    p=[outdir,'/',supertile,'/',tile,'_',res,'_reg.mat'];
    p2=[outdir,'/',supertile,'/',tile,'_',res,'.mat'];
    if exist(p)
        fprintf('Found reg tile %s\n',p)   ;
        tilef = [tilef, p];
    elseif exist(p2)
        fprintf('Found tile %s\n',p2);
        tilef = [tilef, p2];
    else
        fprintf('Tile not found %s\n', tile);
    end
end

batchMergeTileBuffer_old(tilef);

end
