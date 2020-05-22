function dzdt = getTileDzdt(dzdtTileDir,tileName,x0,x1,y0,y1,dx)
% getTileWaterMask retrieves raster water mask for given tile and range
%
% dzdt = getTileWaterMask(dzdtTileDir,tileName,x0,x1,y0,y1,dx) returns dzdt.x
% dzdt.y dzdt.z mask of ones for dzdt, zeros for water using the water
% mask tiles in  dzdtTileDir for the given tileName ('rr_cc') and the
% coordinate range x0,x1,y0,y1 - reads neigboring tiles and pulls needed
% data for buffers.

tileNum = strsplit(tileName,'_');
tileCol = str2num(tileNum{1});
tileRow = str2num(tileNum{2});
tileCol = [tileCol;tileCol  ;tileCol  ;tileCol+1; tileCol-1; tileCol+1; tileCol+1; tileCol-1; tileCol-1];
tileRow = [tileRow;tileRow+1;tileRow-1;tileRow  ; tileRow  ; tileRow+1; tileRow-1; tileRow-1; tileRow+1];

dzdt.x = x0:dx:x1;
dzdt.y = y1:-dx:y0;
dzdt.z = zeros(length(dzdt.y),length(dzdt.x),'single');

%%
i=1;
for i=1:length(tileRow)
    
 
    dzdtTileName=[dzdtTileDir,'/',sprintf('%02d',tileCol(i)),'_',...
        sprintf('%02d',tileRow(i)),'_dzdt.tif'];

       
        if ~exist(dzdtTileName,'file')
            continue
        end
        


    dzdt0 = readGeotiff(dzdtTileName,'map_subset',[x0 x1 y0 y1]);
   
    
    % get corner indexes for this grid in the master
    col0  = find(dzdt0.x(1) == dzdt.x);
    col1  = find(dzdt0.x(end) == dzdt.x);
    row0  = find(dzdt0.y(1) == dzdt.y);
    row1  = find(dzdt0.y(end) == dzdt.y);
    
    if isempty(col0) || isempty(col1) || isempty(row0) || isempty(row1) 
        continue
    end


        dzdt.z(row0:row1,col0:col1) = dzdt0.z; 

        
end

%dzdt.z = dzdt.z==0;
