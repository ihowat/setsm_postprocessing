function data = getDataFromTileAndNeighbors(tileFile,x0,x1,y0,y1,dx,nodataVal)
% getDataFromTileAndNeighbors retrieves raster data for given tile and range
%
% data = getDataFromTileAndNeighbors(tileDir,tileName,tileSuffix,x0,x1,y0,y1,dx,nodataVal) returns
% data.x data.y data.z using tiles in tileDir for the given tileName ('rr_cc') and the
% coordinate range x0,x1,y0,y1 - reads neigboring tiles and pulls needed
% data for buffers.

if ~isfile(tileFile)
    error('tileFile does not exist: %s\n', tileFile)
end

[tileDir,tileBasename,tileExt] = fileparts(tileFile);
tileFname = [tileBasename,tileExt];

sl = split(tileFname,'_');
if startsWith(tileFname,'utm')
    tilePrefix = [sl{1},'_'];
    tileName = strjoin(sl(2:3),'_');
    tileSuffix = ['_',strjoin(sl(4:end),'_')];
else
    tilePrefix = '';
    tileName = strjoin(sl(1:2),'_');
    tileSuffix = ['_',strjoin(sl(3:end),'_')];
end

tileNum = strsplit(tileName,'_');
tileCol = str2num(tileNum{1});
tileRow = str2num(tileNum{2});
tileCol = [tileCol; tileCol  ; tileCol  ; tileCol+1; tileCol-1; tileCol+1; tileCol+1; tileCol-1; tileCol-1];
tileRow = [tileRow; tileRow+1; tileRow-1; tileRow  ; tileRow  ; tileRow+1; tileRow-1; tileRow-1; tileRow+1];

data.x = x0:dx:x1;
data.y = y1:-dx:y0;
if isnan(nodataVal)
    data.z = NaN(length(data.y), length(data.x));
else
    data.z = ones(length(data.y), length(data.x)) * nodataVal;
end

if ~isfolder(tileDir)
    error('tileDir does not exist: %s\n', tileDir)
end

%%
i=1;
for i=1:length(tileRow)
    
    tileFile=[tileDir,'/',tilePrefix,sprintf('%02d',tileCol(i)),'_',...
        sprintf('%02d',tileRow(i)),tileSuffix];
        
    if ~exist(tileFile,'file')
        if i == 1
            error('Reconstructed argument tileFile does not exist: %s', tileFile)
        end
        continue
    end
    
    data0 = readGeotiff(tileFile,'map_subset',[x0 x1 y0 y1]);
     
    % get corner indexes for this grid in the master
    col0 = find(data0.x(1) == data.x);
    col1 = find(data0.x(end) == data.x);
    row0 = find(data0.y(1) == data.y);
    row1 = find(data0.y(end) == data.y);
    
    if isempty(col0) || isempty(col1) || isempty(row0) || isempty(row1) 
        continue
    end
    
    data.z(row0:row1,col0:col1) = data0.z;
        
end
