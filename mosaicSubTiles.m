function mosaicSubTiles(varargin)
% mosaicSubTiles mosaic subtiles and write mat and geotiff output
%
%%%%% ONLY COMPATIBLE WITH 10K, 1km x 1km subtiles!!!!!
%
%
% mosaicSubTiles(subTileDir,dx,outName) mosaics all of the
% subtile .mat files in the directory subTileDir into a mosaic with grid
% resolution dx (2 or 10) and writes matfile output to outName, with tiff
% output as strrep(outName,'.mat','.tif'). Mosaic extend deteremined
% from subtile extents.
%
% mosaicSubTiles(...,'quadrant') only mosaics the subtiles with the quadrant of
% the full tile given by the x_y quad string '1_1','1_2','2_1','2_2'
%
% version 2:
% Subtile alignment from universal LSQ adjustment vertical offsets at
% subtile buffers. Large errors/offsets are filtered. Subtiles that fail
% adjustment are added afterward with the median offset from the mosaic
% removed.

subTileDir = varargin{1};
dx = varargin{2};
outName = varargin{3};

n = find(strcmpi('quadrant',varargin));
if ~isempty(n)
    quadrant = varargin{n+1};
end

n = find(strcmpi('version',varargin));
if ~isempty(n)
    version = varargin{n+1};
    fprintf('Version: %s\n', version)
end

projection = '';
n = find(strcmpi('projection',varargin));
if ~isempty(n)
    projection = varargin{n+1};
else
    error("'projection' argument must be provided");
end

n = find(strcmpi('extent',varargin));
if ~isempty(n)
    x0 = varargin{n+1}(1);
    x1 = varargin{n+1}(2);
    y0 = varargin{n+1}(3);
    y1 = varargin{n+1}(4);
end

fprintf('Indexing subtiles\n')

% make a cellstr of resolved subtile filenames
subTileFiles=dir([subTileDir,'/*_',num2str(dx),'m.mat']);
if isempty(subTileFiles)
    fprintf('No files found matching %s, skipping\n',...
        [subTileDir,'/*_',num2str(dx),'m.mat']);
    return
end
subTileFiles=cellfun( @(x) [subTileDir,'/',x],{subTileFiles.name},'uniformoutput',0);

% Get column-wise number of subtile from file names - assumes the subtile
% names is {tilex}_{tily}_{subtilenum}_....
[~,subTileName] = cellfun(@fileparts,subTileFiles,'uniformoutput',0);
subTileName=cellfun(@(x) strsplit(x,'_'),subTileName,'uniformoutput',0);
if length(subTileName) > 0 && startsWith(subTileName{1}{1},'utm')
    subTileNum = cellfun(@(x) str2num(x{4}),subTileName);
else
    subTileNum = cellfun(@(x) str2num(x{3}),subTileName);
end

%% Get tile projection information, esp. from UTM tile name
%[tileProjName,projection] = getProjName(subTileName{1}{1},projection);

% sort subtilefiles by ascending subtile number order
[subTileNum,n] = sort(subTileNum);
subTileFiles = subTileFiles(n);

if exist('quadrant','var')
    fprintf('selecting subtiles in quadrant %s\n',quadrant)
    % select subtiles by quadrant
    switch quadrant
        case lower('1_1') %lower left quadrant
            n = mod(subTileNum,100) > 0 & mod(subTileNum,100) <= 51 &...
                subTileNum <= 5151;
        case lower('2_1') %upper left quadrant
            n =  (mod(subTileNum,100) == 0 | mod(subTileNum,100) >= 50) &...
                subTileNum <= 5200;
        case lower('1_2') %lower right quadrant
            n = mod(subTileNum,100) > 0 & mod(subTileNum,100) <= 51 &...
                subTileNum >= 5001;
        case lower('2_2') %upper right quadrant
            n =  (mod(subTileNum,100) == 0 | mod(subTileNum,100) >= 50) & subTileNum >= 5050;
        otherwise
            error('quadrant string not recognized')
    end
    
    if ~any(n)
        fprintf('No subtiles within quadrant %s\n',quadrant)
        return
    end
    
    subTileNum = subTileNum(n);
    subTileFiles = subTileFiles(n);
end

NsubTileFiles = length(subTileFiles);

% find buffer size from neighboring tiles
n = diff(subTileNum);
n(mod(subTileNum(1:end-1),100) == 0) = 0;
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
dZ = getOffsets(subTileFiles,subTileNum,buff,outName);


if ~exist('x0','var')
    % get extent of mosaic from tiles at edges
    % lower-left subtile
    m = matfile(subTileFiles{1});
    x0 = m.x(1,1);

    % upper-right subtile
    m = matfile(subTileFiles{end});
    x1 = m.x(1,end);

    % find upper right subtile from each 100th file
    mod100subTileNum = mod(subTileNum,100);
    mod100subTileNum(mod100subTileNum == 0) = 100;

    [~,nMinRow] =  min(mod100subTileNum);
    m = matfile(subTileFiles{nMinRow});
    y0 = m.y(end,1);

    [~,nMaxRow] =  max(mod100subTileNum);
    m = matfile(subTileFiles{nMaxRow});
    y1 = m.y(1,1);

end

% make a polyshape out of boundary for checking subtile overlap
tilePoly = polyshape([x0 x0 x1 x1]',[y0 y1 y1 y0]');

% build tile coordinate vectors
x = x0:dx:x1;
y = y1:-dx:y0;

% setup dataqueue
globalz = nan(length(y),length(x),'single');
globalN = zeros(length(y),length(x),'uint8');
globalNmt = zeros(length(y),length(x),'uint8');
globalz_mad = globalz;
globaltmax = zeros(length(y),length(x),'uint16');
globaltmin = zeros(length(y),length(x),'uint16');
globalNn = 0;
q = parallel.pool.DataQueue;
li = q.afterEach(@(lz) mergemosaic(lz));

% build tile output arrays
spmd
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
rem = mod(NsubTileFiles, numlabs);
lN = fix(NsubTileFiles/numlabs);
if(labindex <= rem)
    start = (labindex-1) * (lN+1)+1;
    stop = start + lN;
else
    start = (labindex-1) * lN + rem+1;
    stop = start + lN-1;
end

for filen=start:stop
    
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
    
    % place the blended substile into the tile
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
    
end % for filen=1:NsubTileFiles
send(q, {z, N, Nmt, z_mad, tmax, tmin});
end %spmd

while q.QueueLength > 0
    pause(0.1);
end
delete(li);
z = globalz;
z_mad = globalz_mad;
N = globalN;
Nmt = globalNmt;
tmax = globaltmax;
tmin = globaltmin;
Nn = globalNn;

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
    
    fprintf('%d remaining subtiles, attempting to add: %s\n',length(nf),subTileFiles{nf(count)})
    
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
    
    % place the blended substile into the tile
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
    
end % while ~isempty(nf)

%% Write Output
% save matfile outputs
if exist('version','var')
    save(outName,'x','y','z','N','Nmt','z_mad','tmax','tmin','version','-v7.3')
else
    save(outName,'x','y','z','N','Nmt','z_mad','tmax','tmin','-v7.3')
end

% write tiff files
z(isnan(z)) = -9999;
outNameTif = strrep(outName,'.mat','_dem.tif');
writeGeotiff(outNameTif,x,y,z,4,-9999,projection)

gdalpath =[]; %set to the path of the gdal binary if not in system path.
if ismac
    gdalpath = '/Library/Frameworks/GDAL.framework/Versions/Current/Programs/';
end
system([gdalpath ,'$BWPY_PREFIX gdaldem hillshade -z 4 -compute_edges  -co TILED=YES -co BIGTIFF=IF_SAFER -co COMPRESS=LZW ',...
   outNameTif,' ',strrep(outNameTif,'_dem.tif','_browse.tif')]);

outNameTif = strrep(outName,'.mat','_count.tif');
writeGeotiff(outNameTif,x,y,N,1,0,projection)

outNameTif = strrep(outName,'.mat','_countmt.tif');
writeGeotiff(outNameTif,x,y,Nmt,1,0,projection)

z_mad(isnan(z_mad)) = -9999;
outNameTif = strrep(outName,'.mat','_mad.tif');
writeGeotiff(outNameTif,x,y,z_mad,4,-9999,projection)

outNameTif = strrep(outName,'.mat','_maxdate.tif');
writeGeotiff(outNameTif,x,y,tmax,2,0,projection)

outNameTif = strrep(outName,'.mat','_mindate.tif');
writeGeotiff(outNameTif,x,y,tmin,2,0,projection)

% Build mosaic centent list and write meta.txt
fprintf('Building meta.txt\n')

if exist('quadrant','var')
    addInfoToSubtileMosaic(subTileDir,dx,outName,quadrant);
else
    addInfoToSubtileMosaic(subTileDir,dx,outName);
end

tileMetav4(outName)

function mergemosaic(x)
lz = x{1};
lN = x{2};
lNmt = x{3};
lz_mad = x{4};
ltmax = x{5};
ltmin = x{6};
n_overlap = ~isnan(globalz(:)) & ~isnan(lz(:)); % & zsub.land(:);
if any(n_overlap)
    lz(isnan(lz) & ~isnan(globalz)) = globalz(isnan(lz) & ~isnan(globalz));
    buffA = single(~(~isnan(globalz) & ~isnan(lz)));
    buffA(~buffA) = NaN;
    buffA(1, isnan(buffA(1,:))) = 0;
    buffA(end,isnan(buffA(end,:))) = 0;
    buffA(isnan(buffA(:,1)),1) = 0;
    buffA(isnan(buffA(:,end)),end) = 0;

    buffA = inpaint_nans(double(buffA), 2);

    notMissing = ~isnan(globalz);

    lz(notMissing) = lz(notMissing).*buffA(notMissing)+...
        globalz(notMissing).*(1-buffA(notMissing));
    lz_mad(notMissing) = lz_mad(notMissing).*buffA(notMissing)+...
        globalz_mad(notMissing).*(1-buffA(notMissing));
    globalz = lz;
    globalz_mad = lz_mad;
else
    globalz(~isnan(lz)) = lz(~isnan(lz));
    globalz_mad(~isnan(lz(:))) = lz_mad(~isnan(lz(:)));
end % if any(n_overlap)
globalN(globalN==0 & lN ~= 0) = lN(globalN==0 & lN~=0);
globalNmt(globalNmt==0 & lNmt ~= 0) = lN(globalNmt==0 & lNmt~=0);
globaltmax(globaltmax==0 & ltmax ~= 0) = ltmax(globaltmax==0 & ltmax~=0);
globaltmin(globaltmin==0 & ltmin ~= 0) = ltmin(globaltmin==0 & ltmin~=0);
Nn1 = sum(~isnan(globalz(:)));
if Nn1 < globalNn
    error('more nans in mosaic with this iteration');
end
globalNn = Nn1;
end % mergemosaic
end % mosaicSubTiles

function dZ = getOffsets(subTileFiles,subTileNum,buff,outName)

NsubTileFiles = length(subTileFiles);

% output file for coregistration offsets so that we don't need to
% recalculate them all if we want to change something
regFile=strrep(outName,'.mat','tileReg.mat');

% check if coregistration file exists and load it if so
if exist(regFile,'file')
    load(regFile)
else
    % no coregistration file, so calculate all suntile offsets
    
    % iteration output variables
    nrt = nan(NsubTileFiles-1,1);
    dzup = nan(NsubTileFiles-1,1);
    dzrt = nan(NsubTileFiles-1,1);
    dzup_mad = nan(NsubTileFiles-1,1);
    dzrt_mad = nan(NsubTileFiles-1,1);
    N=nan(NsubTileFiles,1);
    
    parfor n = 1:NsubTileFiles-1
        
        fprintf('subtile %d of %d ',n,NsubTileFiles-1)
        
        % make sure this subtile has a za_med var
        if  ~ismember('za_med',who('-file',subTileFiles{n}))
            fprintf('variable za_med doesn''t exist, skipping\n')
            continue
        end
        
        m0 = matfile(subTileFiles{n});
        N(n) = max(max(m0.N));
        
        % check if not top row and the up neighbor subtile exists
        if mod(subTileNum(n),100) ~= 0 && ...
                subTileNum(n)+1 == subTileNum(n+1)
            
            % check up neighbor has za_med var
            if ismember('za_med',who('-file',subTileFiles{n+1}))
                
                % load top buffer of the bottom (nth) subtile of pair
                % check to make sure buffer is > 10% land
                l0 = m0.land(2:2*buff,2:end-1);
                if sum(l0(:))/numel(l0) > 0.1
                    z0 = m0.za_med(2:2*buff,2:end-1);
                    z0(~l0) = NaN;
                    y0 = m0.y(2:2*buff,:);
                    x0 = m0.x(:,2:end-1);
                    
                    % load bottom buffer of the top (nth+1) subtile of pair
                    m1 = matfile(subTileFiles{n+1});
                    z1 = m1.za_med(end-(2*buff-1):end-1,2:end-1);
                    y1 = m1.y(end-(2*buff-1):end-1,:);
                    x1 = m1.x(:,2:end-1);
                    
                    % check dimensions of z0 and z0 buffers consistent
                    if ~any(y0 ~= y1) && ~any(x0 ~= x1)
                        dzn = z0(:)-z1(:);
                        if sum(~isnan(dzn))/numel(dzn) > 0.1
                            dzup(n) = mynanmedian(dzn);
                            dzup_mad(n) = mad(dzn,1);
                        end
                    else
                        warning('inconsistent buffers, offset not computed')
                    end
                end
            end
        end
        
        % check right neighbor exists
        nrtn = find(subTileNum(n)+100 ==  subTileNum);
        
        if isempty(nrtn)
            fprintf('\n')
            continue
        end
        
        nrt(n) = nrtn;
        
        % check right neighbor has za_med var
        if ismember('za_med',who('-file',subTileFiles{nrt(n)}))
            
            % check to make sure buffer is > 10% land
            l0 = m0.land(2:end-1,end-(2*buff-1):end-1);
            if sum(l0(:))/numel(l0) > 0.1
                
                % load right buffer of the left (nth) subtile of pair
                z0 = m0.za_med(2:end-1,end-(2*buff-1):end-1);
                z0(~l0) = NaN;
                y0 = m0.y(2:end-1,:);
                x0 = m0.x(:,end-(2*buff-1):end-1);
                
                % load left buffer of the right subtile of pair
                m1 = matfile(subTileFiles{nrt(n)});
                z1 = m1.za_med(2:end-1,2:(2*buff));
                y1 = m1.y(2:end-1,:);
                x1 = m1.x(:,2:(2*buff));
                
                % check dimensions of z0 and z0 buffers consistent
                if ~any(y0 ~= y1) && ~any(x0 ~= x1)
                    dzn = z0(:)-z1(:);
                    if sum(~isnan(dzn))/numel(dzn) > 0.1
                        dzrt(n) = mynanmedian(dzn);
                        dzrt_mad(n) = mad(dzn,1);
                    end
                else
                    warning('inconsistent buffers, offset not computed')
                end
            end
        end
        
        fprintf('\n')
    end
    
    % get N for the last file
    n = NsubTileFiles;
    if ismember('N',who('-file',subTileFiles{n}))
        m0 = matfile(subTileFiles{n});
        N(n) = max(max(m0.N));
    end
    
    save(regFile,'nrt','dzup','dzrt','dzup_mad','dzrt_mad','N')
    
end

% adjustment
%build pair indexes
n1 = [(1:NsubTileFiles-1)';(1:NsubTileFiles-1)'];
n2 = [(2:NsubTileFiles)';nrt];

% concoctenate offsets/errors
dz = [dzup;dzrt];
dze = [dzup_mad;dzrt_mad];

% remove nan pairs or with large offsets and/or mad values
n = any(isnan([n1 n2 dz dze]),2);
n = n | abs(dz) > 50 | dze > 2;

n1(n) =[];
n2(n) =[];
dz(n) =[];
dze(n) =[];

Npairs = length(n1);

% Build design and weight matrices
A = zeros(Npairs,NsubTileFiles); % initialize design matrix

linearInd = sub2ind([Npairs NsubTileFiles], (1:Npairs)', n1);
A(linearInd) = 1;
linearInd = sub2ind([Npairs NsubTileFiles], (1:Npairs)', n2);
A(linearInd) = -1;

% locate missing tiles
n_missing  = ~any(A) | isnan(N)';
%
% remove missing tiles
A(:,n_missing) = [];

% % add delta=0 lines
A = [A;diag(ones(1,size(A,2)))];
dz = [dz;zeros(size(A,2),1)];
% dze = [dze;ones(size(A,2),1).*4];

dze = [ones(size(dze));ones(size(A,2),1).*100];

% dze = ones(size(dze),'single');
% dze = [dze;1./sqrt(N(~n_missing))];

% calculate weights
wz = 1./dze.^2;

% build offset vector
dZ = nan(NsubTileFiles,1);

fprintf('performing LSQ adjustment\n')

dZ(~n_missing) = (wz.*A)\(wz.*dz);

save(strrep(outName,'.mat','tileReg.mat'),'dZ','-append');

% test correction
% dz = z0-z1;
%
% dzn = (z0-dZ0)-(z1-dZ1);
%
% dzn = z0-dZ0-z1+dZ1;
%
% dzn = (z0-z1)-dZ0+dZ1;
%
% dzn = dz-dZ0+dZ1;
%
% dzn = dz - dZ(n1) + dZ(n2);

end