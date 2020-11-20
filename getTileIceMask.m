function ice = getTileIceMask(iceTileDir,tileName,x0,x1,y0,y1,dx)
% getTileWaterMask retrieves raster water mask for given tile and range
%
% ice = getTileWaterMask(iceTileDir,tileName,x0,x1,y0,y1,dx) returns ice.x
% ice.y ice.z mask of ones for ice, zeros for water using the water
% mask tiles in  iceTileDir for the given tileName ('rr_cc') and the
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

ice.x = x0:dx:x1;
ice.y = y1:-dx:y0;
ice.z = false(length(ice.y),length(ice.x));

if ~isfolder(iceTileDir)
    fprintf('Error: iceTileDir does not exist: %s\n', iceTileDir)
    fprintf('Assuming all area is not ice, please fix!\n')
    return
end

%%
i=1;
for i=1:length(tileRow)
    
 
    iceTileName=[iceTileDir,'/',tilePrefix,sprintf('%02d',tileCol(i)),'_',...
        sprintf('%02d',tileRow(i)),'_ice.tif'];

       
        if ~exist(iceTileName,'file')
            continue
        end
        


    ice0 = readGeotiff(iceTileName,'map_subset',[x0 x1 y0 y1]);
   
    
    % get corner indexes for this grid in the master
    col0  = find(ice0.x(1) == ice.x);
    col1  = find(ice0.x(end) == ice.x);
    row0  = find(ice0.y(1) == ice.y);
    row1  = find(ice0.y(end) == ice.y);
    
    if isempty(col0) || isempty(col1) || isempty(row0) || isempty(row1) 
        continue
    end


        ice.z(row0:row1,col0:col1) = ice0.z; 

        
end

%ice.z = ice.z==0;
