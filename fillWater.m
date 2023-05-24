function [z_filled,water_mask,dz_filled,water_height,M0_waterblobs_forced] = fillWater(tileMatFile,maskFile,refDemFile,varargin)
%input args
% tileMatFile
% maskefile
% refDemFile

% output5
% z_filled

% %externtal functs
% readGeotiff
% tile
% retile
% inpaint_nans

M0_waterblobs_forced = [];

logMessage('Loading data at beginning of fillWater');

z = [];
z_mad = [];
N = [];
x = [];
y = [];
slope_filter_mask = [];
fillWaterInterpMethod = 2;
if length(varargin) >= 1
    z = varargin{1};
end
if length(varargin) >= 3
    z_mad = varargin{2};
    N = varargin{3};
end
if length(varargin) >= 5
    x = varargin{4};
    y = varargin{5};
end
if length(varargin) >= 6
    slope_filter_mask = varargin{6};
end
if length(varargin) >= 7
    fillWaterInterpMethod = varargin{7};
end

if any(cellfun(@isempty,{z, z_mad, N, x, y}))
    % read tile DEM mat file
    if isstr(tileMatFile) % it's a string, might be a filename
        if exist(tileMatFile,'file') % yup its a file
            m = matfile(tileMatFile); % load it
        else
            error('tileMatFile does not exist: %s', tileMatFile);
        end
    elseif isvalid(tileMatFile) % not a string, is it a valid file handle?
        m = tileMatFile;
    else
        error('tileMatFile arg must be a filename or valid matfile handle');
    end
    if isempty(z)
        z = m.z;
    end
    if isempty(z_mad)
        z_mad = m.z_mad;
    end
    if isempty(N)
        N = m.N;
    end
    if isempty(x)
        x = m.x;
    end
    if isempty(y)
        y = m.y;
    end
end

dx = x(2)-x(1);

% load mask
C = readGeotiff(maskFile);

% set no data (which should be ocean) to water
C.z(C.z == 0) = 80;

% make sure mask is same grid as dem (just do nn interp)
C = interp2(C.x,C.y(:),C.z,x,y(:),'*nearest');

cluster_size_px = ceil(500 / dx^2);

% define water mask using classification and MAD and repeat criteria
M0 = ( C == 80 | C == 0 );
M = M0 & ( (z_mad > 0.3 & N > 4 ) | N <= 4 | isnan(z) );
M0 = bwareaopen(M0, cluster_size_px);
M = bwareaopen(M, cluster_size_px);

clear C z_mad N;

% remove small non water clusters
M = ~bwareaopen(~M, cluster_size_px);
M0 = M0 | M;
saved_good_z_over_water = xor(M, M0);

% read reference DEM
R = readGeotiff(refDemFile);
R.z(R.z < -600) = nan;

% make sure reference DEM is same grid as dem
R = interp2(R.x,R.y(:),R.z,x,y(:),'*bilinear');

dz = z - R;

% Define a buffer zone near the edge of water in the water mask.
water_mask_uncertain_px = ceil(100 / dx);
water_mask_uncertain_zone = xor(imdilate(M, ones(water_mask_uncertain_px*2+1)), imerode(M, ones(water_mask_uncertain_px*2+1)));

% Get the height of water from the reference surface
% (copernicus generally has good level water surfaces)
% by sampling the water height a short distance out from land
% and extending that sample back to shore.
get_water_height_dist_from_shore_px = ceil(150 / dx);
avg_water_height_interp_dist_px = ceil(1000/dx);
get_water_height_mask = imerode(M, circleMaskSE(get_water_height_dist_from_shore_px));
water_height = nan(size(M));
water_height(get_water_height_mask) = R(get_water_height_mask);
water_height = movmean(movmean(water_height, avg_water_height_interp_dist_px*2+1, 1, 'omitnan'), avg_water_height_interp_dist_px*2+1, 2, 'omitnan');
water_height(get_water_height_mask) = R(get_water_height_mask);
water_height = min(water_height, R, 'includenan');
clear get_water_height_mask;

% Add clusters of land to the water mask where the land is within
% a short distance from the water's edge and the z height is below
% both water level and the reference surface.
fill_water_minarea_px = ceil(3500 / dx^2);
fill_water_fillhole_px = ceil(500 / dx^2);
fill_water = M | (water_mask_uncertain_zone & z < water_height & z < R);
fill_water = bwareaopen(fill_water, fill_water_minarea_px);
fill_water = ~bwareaopen(~fill_water, fill_water_fillhole_px);
fill_water = fill_water & water_mask_uncertain_zone;
M(fill_water) = 1;
M0(fill_water) = 1;
clear fill_water water_mask_uncertain_zone;

% apply water and no data (ocean) mask to z
z_masked = z;
z_masked(M) = NaN;
z_nonwater_nans = isnan(z) & ~M;
water_mask = M;
clear z;

% Fill the DEM with reference DEM using the Delta method

% difference map between z_masked and reference DEM with NaNs in voids
dz = z_masked - R;
clear z_masked;

% Force water fill to meet reference surface at a set distance from land.
interp_short_px = ceil(50 / dx);
interp_long_px = ceil(1000 / dx);
long_interp_saved_z_dist_from_shore_px = ceil(100 / dx);

force_ref_zone = ~imdilate(~M, circleMaskSE(interp_short_px));
long_interp_area = saved_good_z_over_water & imerode(M0, circleMaskSE(long_interp_saved_z_dist_from_shore_px));
long_interp_area = imdilate(long_interp_area, circleMaskSE(interp_long_px));

shoreline_area = xor(M, imdilate(M, ones(max(3, floor(10/dx)*2+1))));
dz_shoreline = nan(size(M));
dz_shoreline(shoreline_area) = dz(shoreline_area);
dz_interp_short = abs(dz_shoreline) >= 15.0;
dz_interp_short = bwareaopen(dz_interp_short, floor(60 / dx^2));
dz_interp_short_keep = dz_interp_short;

shoreline_area = xor(M, imdilate(M, ones(3)));
dz_shoreline = nan(size(M));
dz_shoreline(shoreline_area) = dz(shoreline_area);
dz_interp_short = abs(dz_shoreline) >= 15.0;
dz_interp_short = bwareaopen(dz_interp_short, 3);
dz_interp_short_keep = dz_interp_short_keep | dz_interp_short;

clear shoreline_area dz_shoreline dz_interp_short;

long_interp_area(imdilate(dz_interp_short_keep, circleMaskSE(interp_short_px*2))) = 0;
force_ref_zone(long_interp_area) = 0;
dz(force_ref_zone) = 0;

clear force_ref_zone long_interp_area dz_interp_short_keep;
clear saved_good_z_over_water;


% interpolate differences across void, breaking dz into tiles to speed up
A = tile(dz,5,5,200);
clear dz;

M_tiles = tile(M,5,5,200);

% initiate filled cell array of tiles
B = cell(size(A));

% loop through cells to interpolate voids
i=1;
for i=1:numel(A)
    fprintf('%d of %d\n',i, numel(A));
    if all(M_tiles{i}(:))  % all water, fill with mean difference from reference
        B{i} = 0 * ones(size(A{i}),'single');
    elseif ~any(M_tiles{i}(:))  % no water, skip
        B{i} = A{i};
    else  % mix of water and no water, interpolate
%        B{i} = single(inpaint_nans(double(A{i}),2));
        water_border = ~imerode(ones(size(M_tiles{i})), ones(7)) & M_tiles{i};
        water_border(imdilate(~M_tiles{i}, circleMaskSE(floor(300/dx)))) = 0;
        Ai = A{i};
        Ai(water_border) = 0;
        B{i} = single(inpaint_nans(double(Ai),fillWaterInterpMethod));
        clear Ai water_border;
    end
end
clear A M_tiles;

% mosiac filled dz tiles
dz_filled = retile(B,200,'linear');
clear B;

% calculate filled DEM
z_filled = R + dz_filled;
%clear dz_filled;

% Flatten to water level divets in the filled DEM caused by interpolation.
underwater = M & z_filled < water_height;
z_filled(underwater) = water_height(underwater);
M(underwater) = 1;
%clear underwater water_height;
clear underwater;

z_nonwater_nans = z_nonwater_nans & ~M;

% Set back to NaN all original NoData that was outside water mask
% at the beginning of this routine (includes data removed by slope filter).
z_filled(z_nonwater_nans) = nan;

% Remove blobs of filled water that are surrounded by original NoData.
waterblobs_in_nodata = zeros(size(M), 'logical');
kept_z_or_waterblobs = ~z_nonwater_nans;
CC_1 = bwconncomp(kept_z_or_waterblobs);
CC_M = bwconncomp(M);
mapM = containers.Map(cellfun(@(x) x(1), CC_M.PixelIdxList), 1:CC_M.NumObjects);
blob_idxlist = [];
for i = 1:CC_1.NumObjects
    blob_idxlist = CC_1.PixelIdxList{i};
    blob_startidx = blob_idxlist(1);
    if isKey(mapM, blob_startidx) && isequal(blob_idxlist, CC_M.PixelIdxList{mapM(blob_startidx)})
        waterblobs_in_nodata(blob_idxlist) = 1;
    end
end
clear kept_z_or_waterblobs z_nonwater_nans;
clear CC_1 CC_M mapM blob_idxlist;

water_mask = M;
clear M M0;

z_filled(waterblobs_in_nodata) = nan;
water_mask(waterblobs_in_nodata) = 0;
