function land = getTileWaterMask(waterTileDir,tileName,x0,x1,y0,y1)
% getTileWaterMask retrieves raster water mask for given tile and range
%
% land = getTileWaterMask(waterTileDir,tileName,x0,x1,y0,y1) returns land.x
% land.y land.z mask of ones for land, zeros for water using the water
% mask tiles in  waterTileDir for the given tileName ('rr_cc') and the
% coordinate range x0,x1,y0,y1 - reads neigboring tiles and pulls needed
% data for buffers.

tileNum = strsplit(tileName,'_');
tileCol = str2num(tileNum{1});
tileRow = str2num(tileNum{2});
tileCol = [tileCol;tileCol  ;tileCol  ;tileCol+1; tileCol-1; tileCol+1; tileCol+1; tileCol-1; tileCol-1];
tileRow = [tileRow;tileRow+1;tileRow-1;tileRow  ; tileRow  ; tileRow+1; tileRow-1; tileRow-1; tileRow+1];

land.x = x0:10:x1;
land.y = y1:-10:y0;
land.z = false(length(land.y),length(land.x));

i=1;
for i=1:length(tileRow)
    
    % get this tile
    waterTileName=[waterTileDir,'/',num2str(tileCol(i)),'_',num2str(tileRow(i)),'_water.tif'];
    land0 = readGeotiff(waterTileName,'map_subset',[x0 x1 y0 y1]);
    
    % get corner indexes for this grid in the master
    col0  = find(land0.x(1) == land.x);
    col1  = find(land0.x(end) == land.x);
    row0  = find(land0.y(1) == land.y);
    row1  = find(land0.y(end) == land.y);
    
    land.z(row0:row1,col0:col1) = land0.z;
    
end

land.z = land.z==0;