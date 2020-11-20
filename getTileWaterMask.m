function land = getTileWaterMask(waterTileDir,tileName,x0,x1,y0,y1,dx,varargin)
% getTileWaterMask retrieves raster water mask for given tile and range
%
% land = getTileWaterMask(waterTileDir,tileName,x0,x1,y0,y1,dx) returns land.x
% land.y land.z mask of ones for land, zeros for water using the water
% mask tiles in  waterTileDir for the given tileName ('rr_cc') and the
% coordinate range x0,x1,y0,y1 - reads neigboring tiles and pulls needed
% data for buffers.

if startsWith(tileName,'utm')
    sl = split(tileName,'_');
    tilePrefix = [sl{1},'_'];
    tileName_in_tileDef = strjoin(sl(2:3),'_');
else
    tilePrefix = '';
    tileName_in_tileDef = tileName;
end

tileNum = strsplit(tileName_in_tileDef,'_');
tileCol = str2num(tileNum{1});
tileRow = str2num(tileNum{2});
tileCol = [tileCol;tileCol  ;tileCol  ;tileCol+1; tileCol-1; tileCol+1; tileCol+1; tileCol-1; tileCol-1];
tileRow = [tileRow;tileRow+1;tileRow-1;tileRow  ; tileRow  ; tileRow+1; tileRow-1; tileRow-1; tileRow+1];

land.x = x0:dx:x1;
land.y = y1:-dx:y0;
land.z = false(length(land.y),length(land.x));

if ~isfolder(waterTileDir)
    fprintf('Error: waterTileDir does not exist: %s\n', waterTileDir)
    fprintf('Assuming all area is land, please fix!\n')
    land.z = true(length(land.y),length(land.x));
    return
end

includeIceFlag = false;
if any(strcmpi(varargin,'includeIce'))
    includeIceFlag = true;
end
%%
i=1;
for i=1:length(tileRow)
    
    waterFlag=false;
    
    % get this tile
    waterTileName=[waterTileDir,'/',tilePrefix,sprintf('%02d',tileCol(i)),'_',...
        sprintf('%02d',tileRow(i)),'_land.tif'];
    
 
        iceTileName=[waterTileDir,'/',tilePrefix,sprintf('%02d',tileCol(i)),'_',...
            sprintf('%02d',tileRow(i)),'_ice.tif'];
    
    if ~exist(waterTileName,'file')
        
        waterTileName=[waterTileDir,'/',tilePrefix,sprintf('%02d',tileCol(i)),'_',...
            sprintf('%02d',tileRow(i)),'_water.tif'];
       
        if ~exist(waterTileName,'file')
            continue
        end
        
         waterFlag=true;
    end

    land0 = readGeotiff(waterTileName,'map_subset',[x0 x1 y0 y1]);
    
    if includeIceFlag && exist(iceTileName,'file')
           ice0 = readGeotiff(iceTileName,'map_subset',[x0 x1 y0 y1]);
           land0.z = land0.z == 1 | ice0.z == 1;
    end
           
    
    % get corner indexes for this grid in the master
    col0  = find(land0.x(1) == land.x);
    col1  = find(land0.x(end) == land.x);
    row0  = find(land0.y(1) == land.y);
    row1  = find(land0.y(end) == land.y);
    
    if isempty(col0) || isempty(col1) || isempty(row0) || isempty(row1) 
        continue
    end

    if ~waterFlag
        land.z(row0:row1,col0:col1) = land0.z; 
    else
        land.z(row0:row1,col0:col1) = ~land0.z;
    end
        
end

%land.z = land.z==0;
