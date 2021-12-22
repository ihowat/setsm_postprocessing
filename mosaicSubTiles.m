function mosaicSubTiles(varargin)
% mosaicSubTiles mosaic subtiles and write mat and geotiff output
%
% mosaicSubTiles(subTileDir,dx,outName) mosaics all of the
% subtile .mat files in the directory subTileDir into a mosaic with grid
% resolution dx (2 or 10) and writes matfile output to outName, with tiff
% output as strrep(outName,'.mat','.tif'). Mosaic extend deteremined
% from subtile extents or from extents provided as input args. Subtile
% names must be in format row_col_*.mat with the row and column number of
% the subtile in the tile.
%
% Subtile alignment from universal LSQ adjustment vertical offsets at
% subtile buffers. Large errors/offsets are filtered. Subtiles that fail
% adjustment are added afterward with the median offset from the mosaic
% removed.
%
% %%%%% ONLY COMPATIBLE WITH 1km x 1km subtiles!!!!!

subTileDir = varargin{1};
dx = varargin{2};
outName = varargin{3};

fprintf('subTileDir = %s\n',subTileDir)
fprintf('outName = %s\n',outName)

if ~exist(subTileDir,'dir')
    error('Subtiles folder does not exist: %s', subTileDir)
end

projection = '';
n = find(strcmpi('projection',varargin));
if ~isempty(n)
    projection = varargin{n+1};
else
    error("'projection' argument must be provided");
end

n = find(strcmpi('version',varargin));
if ~isempty(n)
    version = varargin{n+1};
    fprintf('Version: %s\n', version)
end

n = find(strcmpi('exportTif',varargin));
if ~isempty(n)
    exportTif = varargin{n+1};
else
    exportTif = true;
end

n = find(strcmpi('extent',varargin));
if ~isempty(n)
    x0 = varargin{n+1}(1);
    x1 = varargin{n+1}(2);
    y0 = varargin{n+1}(3);
    y1 = varargin{n+1}(4);
end

%if output directory doesnt already exist, make it
[outDir,~,~] = fileparts(outName);
if ~exist(outDir,'dir')
    mkdir(outDir)
end

fprintf('Indexing subtiles\n')

% make a cellstr of resolved subtile filenames
subTileFiles=dir([subTileDir,'/*_',num2str(dx),'m.mat']);
if isempty(subTileFiles)
    fprintf('No files found matching %s, skipping\n',...
        [subTileDir,'/*_',num2str(dx),'m.mat']);
    return
end

subTileFiles=cellfun( @(x) [subTileDir,'/',x],{subTileFiles.name},...
    'uniformoutput',0);

% Get rows and cols of subtile from file names - assumes the subtile
% names is *_{row}_{col}_{res}.mat
[~,subTileName] = cellfun(@fileparts,subTileFiles,'uniformoutput',0);
subTileName=cellfun(@(x) strsplit(x,'_'),subTileName,'uniformoutput',0);
subTileCol = cellfun(@(x) str2num(x{end-1}),subTileName);
subTileRow = cellfun(@(x) str2num(x{end-2}),subTileName);

Ncols = max(subTileCol);
Nrows = max(subTileRow);

% Get tile projection information, esp. from UTM tile name
%[tileProjName,projection] = getProjName(subTileName{1}{1},projection);

% column-wise subtile number
subTileNum = sub2ind([Nrows,Ncols],subTileRow,subTileCol);

% % sort subtilefiles by ascending subtile number order
[subTileNum,n] = sort(subTileNum);
subTileCol = subTileCol(n);
subTileRow= subTileRow(n);
subTileFiles = subTileFiles(n);

NsubTileFiles = length(subTileFiles);

% find buffer size from neighboring tiles
n = diff(subTileNum);
n(mod(subTileNum(1:end-1),Nrows) == 0) = 0; %DEPENDS ON SUBTILE ROW/COLS
n = find(n == 1,1,'first');
if ~isempty(n)
    buffcheck1=load(subTileFiles{n},'y');
    buffcheck2=load(subTileFiles{n+1},'y');
    buff = (length(buffcheck2.y)-find(buffcheck1.y(1) == buffcheck2.y))/2;
    buff = round(buff);
else
    fprintf('no neighboring tiles, returning\n')
    return
end

fprintf('performing coregistration & adjustment between adjoining subtiles\n')
dZ = getOffsets(subTileFiles,subTileNum,buff);

% if extent of mosaic not specidied, get from tiles at edges
if ~exist('x0','var')

    % left
    [~,n] = min(subTileCol);
    m = matfile(subTileFiles{n(1)});
    x0 = min(m.x);
    
    % right
    [~,n] = max(subTileCol);
    m = matfile(subTileFiles{n(1)});
    x1 = max(m.x);
    
     % bottom
    [~,n] = min(subTileRow);
    m = matfile(subTileFiles{n(1)});
    y0 = min(m.y);
    
    % top
    [~,n] = max(subTileRow);
    m = matfile(subTileFiles{n(1)});
    y1 = max(m.y);
    
end

% make a polyshape out of boundary for checking subtile overlap
tilePoly = polyshape([x0 x0 x1 x1]',[y0 y1 y1 y0]');

% build tile coordinate vectors
x = x0:dx:x1;
y = y1:-dx:y0;

% build tile output arrays
z = nan(length(y),length(x),'single');
N = zeros(length(y),length(x),'uint8');
Nmt = zeros(length(y),length(x),'uint8');
z_mad = z;
tmax = zeros(length(y),length(x),'uint16');
tmin = zeros(length(y),length(x),'uint16');


% initialize subtile count for use in n-weighted alignment
subtile_n=1;

% initialize count of pixels with data in mosaic for error checking
Nn=0;
for filen=1:NsubTileFiles
    
    % first only add subtiles with adjustements
    if isnan(dZ(filen))
        continue
    end
    
    fprintf('adding subtile %d of %d: %s\n',filen,NsubTileFiles,subTileFiles{filen})
    
    % get list of variables within this subtile mat file
    mvars  = who('-file',subTileFiles{filen});
    
    if ~any(ismember(mvars,'za_med'))
        fprintf('no za_med found, skipping\n')
        continue
    end
    
    % Check if subtile overlaps mosaic boundary:
    % open matfile to load coordinate vectors withoutloading arrays
    m=matfile(subTileFiles{filen});
    
    % make a polyshape out of this subtile boundary
    filenPoly = polyshape([min(m.x) min(m.x) max(m.x) max(m.x)]',...
        [min(m.y) max(m.y) max(m.y) min(m.y)]');
    
    % test overlap
    if ~overlaps(tilePoly,filenPoly)
        fprintf('subtile out of bounds, skipping\n')
        continue
    end
    
    % load subtile into structure
    zsub=load(subTileFiles{filen},'x','y','za_med','land','N','Nmt','za_mad','tmax','tmin');
    
    if ~any(~isnan(zsub.za_med(:)))
        fprintf('all nans, skipping\n')
        continue
    end
    
    % add adjustment offset
    zsub.za_med = zsub.za_med - dZ(filen);
    
    % get corner indexes for this grid in the master
    col0  = find(zsub.x(1) == x);
    col1  = find(zsub.x(end) == x);
    row0  = find(zsub.y(1) == y);
    row1  = find(zsub.y(end) == y);
    
    if isempty(col0) || isempty(col1) || isempty(row0) || isempty(row1)
        fprintf('subtile and tile coordinates dont match, subtile may span boundary,skipping\n')
        continue
    end
    
    % subset the current mosaic by range of subtile
    z1 = z(row0:row1,col0:col1);
    
    z_mad1 = z_mad(row0:row1,col0:col1);
    
    % find overlapping non-nan and non-water pixels
    n_overlap = ~isnan(z1(:)) & ~isnan(zsub.za_med(:)) & zsub.land(:);
    
    % check if overlapping pixels exist to determine if blending is needed
    if any(n_overlap)
        
        % fill missing data in this subtile with current data in mosaic
        % subset
        zsub.za_med(isnan(zsub.za_med) & ~isnan(z1)) =...
            z1(isnan(zsub.za_med) & ~isnan(z1));
        
        
        % also blend za_mad
        zsub.za_mad(isnan(zsub.za_mad) & ~isnan(z_mad1)) =...
            z1(isnan(zsub.za_mad) & ~isnan(z_mad1));
        
        % Create the blending array by setting zeros at the far edge of
        % subtile/tile overlap and ones at the other edge, linearly
        % interpolating between the two.
        % find where data is missing in both subtile and tile
        buffA = single(~(~isnan(zsub.za_med) & ~isnan(z1)));
        
        % set pixels with data in both subtile and tiles to NaN as an
        % interpolation flag for inpaint_nans
        buffA(~buffA) = NaN;
        
        % set boundaries of blend array to zero
        buffA(1,isnan(buffA(1,:))) = 0;
        buffA(end,isnan(buffA(end,:))) = 0;
        
        buffA(isnan(buffA(:,1)),1) = 0;
        buffA(isnan(buffA(:,end)),end) = 0;
        
        % interpolate linearly across NaNs
        buffA=inpaint_nans(double(buffA),2);
        
        % find where there is data in the tile subset
        notMissing = ~isnan(z1);
        
        % blend the subtile and mosaic subset, applying the edge-distance
        % weighting
        zsub.za_med(notMissing) = zsub.za_med(notMissing).*buffA(notMissing) +...
            z1(notMissing).*(1- buffA(notMissing));
        
        
        zsub.za_mad(notMissing) = zsub.za_mad(notMissing).*buffA(notMissing) +...
            z_mad1(notMissing).*(1- buffA(notMissing));
        
    else
        
        % if no pixels overlap, just add the subset data into the subtile,
        % replacing the NaNs
        zsub.za_med(~isnan(z1(:))) = z1(~isnan(z1(:)));
        
        zsub.za_mad(~isnan(z1(:))) = z_mad1(~isnan(z1(:)));
        
    end
    
    % place the belended substile into the tile
    z(row0:row1,col0:col1) = zsub.za_med;
    
    % place N grid into the tile just within the tile borders
    N(row0+buff:row1-buff,col0+buff:col1-buff) =...
        zsub.N(buff+1:end-buff,buff+1:end-buff);
    
    Nmt(row0+buff:row1-buff,col0+buff:col1-buff) =...
        zsub.Nmt(buff+1:end-buff,buff+1:end-buff);
    z_mad(row0:row1,col0:col1) = zsub.za_mad;
    tmax(row0+buff:row1-buff,col0+buff:col1-buff) =...
        zsub.tmax(buff+1:end-buff,buff+1:end-buff);
    tmin(row0+buff:row1-buff,col0+buff:col1-buff) =...
        zsub.tmin(buff+1:end-buff,buff+1:end-buff);
    
    % count the number of pixels with data after this merge
    Nn1 = sum(~isnan(z(:)));
    
    % if the new number of pixels with data is now less than before,
    % then data was overwritten with NaNs, which is an error.
    if Nn1 < Nn
        error('more nans in mosaic with this iteration\n');
    end
    
    % reset count of pixels with data for next iteration
    Nn = Nn1;
    
    % update count of subtiles added to mosaic
    subtile_n= subtile_n+1;
    
end

%% Add dems with nan dZ's

% index vector of offsets with nans. We will remove these as they are added
% to the mosaic (or are found to not contain useable data).
nf = find(isnan(dZ));

length_nf = length(nf);
% count of which dem in the nf vector to attempt to add to the mosaic. If
% there is no overlapping data currently in the mosaic, the count will
% increase to the next dem and so forth, rptating back to 1 once it cycles
% through. If the count goes through a full cycle of 1:length(nf), the
% remaining nf dems will be added without registration.
count=1;

% flag to indicate if that addition of all dems have failed due to lack of
% overlap and to add without registration.
noRegFlag = false;

while ~isempty(nf)
    
    if count > length(nf) % cycle complete, reset counter
        count = 1;
        
        % check if no more subtiles added this cycle
        if length(nf) == length_nf
            noRegFlag = true; % turn off registration if no overlap
        end
        
        % reset length nf of this cycle to compare with next
        length_nf = length(nf);
        
    end
    
    filen = nf(count);
    
    fprintf('%d remainging subtiles, attempting to add: %s\n',length(nf),subTileFiles{nf(count)})
    
    % get list of variables within this subtile mat file
    mvars  = who('-file',subTileFiles{filen});
    
    if ~any(ismember(mvars,'za_med'))
        nf(count) = [];
        fprintf('no za_med found, skipping\n')
        continue
    end
    
    % Check if subtile overlaps mosaic boundary:
    % open matfile to load coordinate vectors withoutloading arrays
    m=matfile(subTileFiles{filen});
    
    % make a polyshape out of this subtile boundary
    filenPoly = polyshape([min(m.x) min(m.x) max(m.x) max(m.x)]',...
        [min(m.y) max(m.y) max(m.y) min(m.y)]');
    
    % test overlap
    if ~overlaps(tilePoly,filenPoly)
        nf(count) = [];
        fprintf('subtile out of bounds, skipping\n')
        continue
    end
    
    % load subtile into structure
    zsub=load(subTileFiles{filen},'x','y','za_med','land','N','Nmt','za_mad','tmax','tmin');
    
    if ~any(~isnan(zsub.za_med(:)))
        nf(count) = [];
        fprintf('all nans, skipping\n')
        continue
    end
    
    % get corner indexes for this grid in the master
    col0  = find(zsub.x(1) == x);
    col1  = find(zsub.x(end) == x);
    row0  = find(zsub.y(1) == y);
    row1  = find(zsub.y(end) == y);
    
    if isempty(col0) || isempty(col1) || isempty(row0) || isempty(row1)
        nf(count) = [];
        fprintf('subtile and tile coordinates dont match, subtile may span boundary,skipping\n')
        continue
    end
    
    % subset the current mosaic by range of subtile
    z1 = z(row0:row1,col0:col1);
    z_mad1 = z_mad(row0:row1,col0:col1);
    
    % find overlapping non-nan and non-water pixels
    n_overlap = ~isnan(z1(:)) & ~isnan(zsub.za_med(:)) & zsub.land(:);
    
    % check if overlapping pixels exist to determine if blending is needed
    if any(n_overlap)
        
        %get median difference between this subtile and mosaic subset
        dz_med = median(zsub.za_med(n_overlap) - z1(n_overlap));
        
        zsub.za_med = zsub.za_med - dz_med; % shift the subtile
        
        % fill missing data in this subtile with current data in mosaic
        % subset
        zsub.za_med(isnan(zsub.za_med) & ~isnan(z1)) =...
            z1(isnan(zsub.za_med) & ~isnan(z1));
        
        % also blend za_mad
        zsub.za_mad(isnan(zsub.za_mad) & ~isnan(z_mad1)) =...
            z1(isnan(zsub.za_mad) & ~isnan(z_mad1));
        
        % Create the blending array by setting zeros at the far edge of
        % subtile/tile overlap and ones at the other edge, linearly
        % interpolating between the two.
        % find where data is missing in both subtile and tile
        buffA = single(~(~isnan(zsub.za_med) & ~isnan(z1)));
        
        % set pixels with data in both subtile and tiles to NaN as an
        % interpolation flag for inpaint_nans
        buffA(~buffA) = NaN;
        
        % set boundaries of blend array to zero
        buffA(1,isnan(buffA(1,:))) = 0;
        buffA(end,isnan(buffA(end,:))) = 0;
        
        buffA(isnan(buffA(:,1)),1) = 0;
        buffA(isnan(buffA(:,end)),end) = 0;
        
        % interpolate linearly across NaNs
        buffA=inpaint_nans(double(buffA),2);
        
        % find where there is data in the tile subset
        notMissing = ~isnan(z1);
        
        % blend the subtile and mosaic subset, applying the edge-distance
        % weighting
        zsub.za_med(notMissing) = zsub.za_med(notMissing).*buffA(notMissing) +...
            z1(notMissing).*(1- buffA(notMissing));
        
        zsub.za_mad(notMissing) = zsub.za_mad(notMissing).*buffA(notMissing) +...
            z_mad1(notMissing).*(1- buffA(notMissing));
        
    else
        
        % no overlap with mosaic, add subtile without registration if
        % already cycled through without additions, as indicated by the
        % noRegFlag
        if noRegFlag
            % if no pixels overlap, just add the subset data into the subtile,
            % replacing the NaNs
            zsub.za_med(~isnan(z1(:))) = z1(~isnan(z1(:)));
            
            zsub.za_mad(~isnan(z1(:))) = z_mad1(~isnan(z1(:)));
            
        else
            % skip this subtile for this cycle
            count = count + 1;
            continue
        end
        
    end
    
    % place the belended substile into the tile
    z(row0:row1,col0:col1) = zsub.za_med;
    
    % place N grid into the tile just within the tile borders
    N(row0+buff:row1-buff,col0+buff:col1-buff) =...
        zsub.N(buff+1:end-buff,buff+1:end-buff);
    
    Nmt(row0+buff:row1-buff,col0+buff:col1-buff) =...
        zsub.Nmt(buff+1:end-buff,buff+1:end-buff);
    z_mad(row0:row1,col0:col1) = zsub.za_mad;
    tmax(row0+buff:row1-buff,col0+buff:col1-buff) =...
        zsub.tmax(buff+1:end-buff,buff+1:end-buff);
    tmin(row0+buff:row1-buff,col0+buff:col1-buff) =...
        zsub.tmin(buff+1:end-buff,buff+1:end-buff);
    
    
    % count the number of pixels with data after this merge
    Nn1 = sum(~isnan(z(:)));
    
    % if the new number of pixels with data is now less than before,
    % then data was overwritten with NaNs, which is an error.
    if Nn1 < Nn
        error('more nans in mosaic with this iteration\n');
    end
    
    % reset count of pixels with data for next iteration
    Nn = Nn1;
    
    % remove this subtile from index
    nf(count) = [];
    
end

%% Write Output
% save matfile outputs
if exist('version','var')
    save(outName,'x','y','z','N','Nmt','z_mad','tmax','tmin','version','-v7.3')
else
    save(outName,'x','y','z','N','Nmt','z_mad','tmax','tmin','-v7.3')
end

% use the tile cropping logic from writeTileToTifv4
writeTileToTifv4(outName, projection, 'browseOnly', ~exportTif)

%% write tiff files
%z(isnan(z)) = -9999;
%outNameTif = strrep(outName,'.mat','_dem.tif');
%writeGeotiff(outNameTif,x,y,z,4,-9999,projection)
%
%gdalpath =[]; %set to the path of the gdal binary if not in system path.
%if ismac
%    gdalpath = '/Library/Frameworks/GDAL.framework/Versions/Current/Programs/';
%end
%% make browse hillshade at 10m resolution
%if dx == 2
%    outNameTemp = strrep(outNameTif,'_dem.tif','_temp.tif');
%    system(['gdal_translate -q -tr 10 10 -r bilinear -co bigtiff=if_safer -a_nodata -9999 ',...
%            outNameTif,' ', outNameTemp]);
%else
%    outNameTemp = outNameTif;
%end
%system([gdalpath ,'gdaldem hillshade -z 3 -compute_edges  -co TILED=YES -co BIGTIFF=IF_SAFER -co COMPRESS=LZW ',...
%   outNameTemp,' ',strrep(outNameTif,'_dem.tif','_browse.tif')]);
%if dx == 2
%    delete(outNameTemp);
%end
%
%outNameTif = strrep(outName,'.mat','_count.tif');
%writeGeotiff(outNameTif,x,y,N,1,0,projection)
%
%outNameTif = strrep(outName,'.mat','_countmt.tif');
%writeGeotiff(outNameTif,x,y,Nmt,1,0,projection)
%
%z_mad(isnan(z_mad)) = -9999;
%outNameTif = strrep(outName,'.mat','_mad.tif');
%writeGeotiff(outNameTif,x,y,z_mad,4,-9999,projection)
%
%outNameTif = strrep(outName,'.mat','_maxdate.tif');
%writeGeotiff(outNameTif,x,y,tmax,2,0,projection)
%
%outNameTif = strrep(outName,'.mat','_mindate.tif');
%writeGeotiff(outNameTif,x,y,tmin,2,0,projection)

% Build mosaic centent list and write meta.txt
fprintf('Gathering subtile information for tile meta.txt file\n')
if exist('quadrant','var')
    addInfoToSubtileMosaic(subTileDir,dx,outName,quadrant);
else
    addInfoToSubtileMosaic(subTileDir,dx,outName);
end
%if ~browseOnly
    fprintf('Writing meta.txt\n')
    tileMetav4(outName)
%end
