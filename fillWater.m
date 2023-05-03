function [z_filled,water_mask] = fillWater(tileMatFile,maskFile,refDemFile,varargin)
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

z = [];
z_mad = [];
N = [];
x = [];
y = [];
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

if any(cellfun(@isempty,{z, z_mad, N, x, y}))
    % read tile DEM mat file
    m = matfile(tileMatFile);
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

cluster_size_px = double(int32(500 / dx^2));

% define water mask using classification and MAD and repeat criteria
M = ( C == 80 | C == 0 );
M = bwareaopen(M, cluster_size_px);
M = M & ( (z_mad > 0.3 & N > 4 ) | N <= 4 | isnan(z) );

% remove small non water clusters
M = ~bwareaopen(~M, cluster_size_px);

% read reference DEM
R = readGeotiff(refDemFile);
R.z(R.z < -600) = NaN;

% make sure reference DEM is same grid as dem
R = interp2(R.x,R.y(:),R.z,x,y(:),'*bilinear');

dz = z - R;

% Define a buffer zone near the edge of water in the water mask.
water_mask_uncertain_px = int32(100 / dx);
water_mask_uncertain_zone = xor(imdilate(M, ones(water_mask_uncertain_px * 2)), imerode(M, ones(water_mask_uncertain_px * 2)));

% Get the height of water from the reference surface
% (copernicus generally has good level water surfaces)
% by sampling the water height a short distance out from land
% and extending that sample back to shore.
get_water_height_erode_px = int32(500 / dx);
get_water_height_mask = imerode(M, ones(get_water_height_erode_px * 2));
water_height = nan(size(M));
water_height(get_water_height_mask) = R(get_water_height_mask);
water_height = movmean(movmean(water_height, (get_water_height_erode_px+water_mask_uncertain_px+1)*2, 1, 'omitnan'), (get_water_height_erode_px+water_mask_uncertain_px+1)*2, 2, 'omitnan');

% Add clusters of land to the water mask where the land is within
% a short distance from the water's edge and the z height is below
% both water level and the reference surface.
fill_water_minarea_px = double(int32(3500 / dx^2));
fill_water_fillhole_px = double(int32(500 / dx^2));
fill_water = M | (water_mask_uncertain_zone & z < water_height & z < R);
fill_water = bwareaopen(fill_water, fill_water_minarea_px);
fill_water = fill_water | (water_mask_uncertain_zone & ~bwareaopen(~fill_water, fill_water_fillhole_px));
M(water_mask_uncertain_zone & fill_water) = 1;

% apply water and no data (ocean) mask to z
z_masked = z;
z_masked(M) = NaN;
z_nonwater_nans = isnan(z) & ~M;
water_mask = M;
clear z;

% Fill the DEM with reference DEM using the Delta method

% difference map between z_masked and reference DEM with NaNs in voids
dz = z_masked - R;

% Force water fill to meet reference surface at a set distance from land.
interp_buff_px = int32(50 / dx);
force_ref_zone = ~imdilate(~M, ones(interp_buff_px * 2));
dz(force_ref_zone) = 0;

% interpolate differences across void, breaking dz into tiles to speed up
A = tile(dz,5,5,200);
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
        B{i} = single(inpaint_nans(double(A{i}),2));
    end
end
clear A;

% mosiac filled dz tiles
dz_filled = retile(B,200,'linear');

% calculate filled DEM
z_filled = R + dz_filled;

% Flatten to water level divets in the filled DEM caused by interpolation.
flatten_water_zone = xor(M, imerode(M, ones(water_mask_uncertain_px * 2)));
underwater = flatten_water_zone & z_filled < water_height;
z_filled(underwater) = water_height(underwater);
M(underwater) = 1;

% Set back to NaN all original NoData that was outside water mask
% at the beginning of this routine (includes data removed by slope filter).
z_nonwater_nans = z_nonwater_nans & ~M;
z_filled(z_nonwater_nans) = NaN;
