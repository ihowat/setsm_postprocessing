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
        logMessage('Loading m.z');
        z = m.z;
    end
    if isempty(z_mad)
        logMessage('Loading m.z_mad');
        z_mad = m.z_mad;
    end
    if isempty(N)
        logMessage('Loading m.N');
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
logMessage(sprintf('Loading water cover maskfile: %s', maskFile));
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
logMessage(sprintf('Loading reference DEM: %s', refDemFile));
R = readGeotiff(refDemFile);
R.z(R.z < -600) = nan;

% make sure reference DEM is same grid as dem
logMessage('Interpolating reference DEM to z grid');
R = interp2(R.x,R.y(:),R.z,x,y(:),'*bilinear');


% Get the height of water from the reference surface
% (copernicus generally has good level water surfaces)
% by sampling the water height a short distance out from land
% and extending that sample back to shore.
logMessage('Getting "water height" from reference surface');
%get_water_height_dist_from_shore_px = ceil(150 / dx);
%avg_water_height_interp_dist_px = ceil(1000 / dx);
get_water_height_dist_from_shore_px = ceil(600 / dx);
avg_water_height_interp_dist_px = ceil(1500 / dx);
get_water_height_mask = imerode(M0, circleMaskSE(get_water_height_dist_from_shore_px));
water_height = nan(size(M0));
water_height(get_water_height_mask) = R(get_water_height_mask);
water_height = movmean(movmean(water_height, avg_water_height_interp_dist_px*2+1, 1, 'omitnan'), avg_water_height_interp_dist_px*2+1, 2, 'omitnan');
%water_height(get_water_height_mask) = R(get_water_height_mask);
water_height = min(water_height, R, 'includenan');
clear get_water_height_mask;


M_waterblobs_forced = [];
M0_waterblobs_forced = [];

add_to_water_mask = true;

get_M_waterblobs_forced = true;
get_M0_waterblobs_forced = true;
if get_M0_waterblobs_forced
    get_M_waterblobs_forced = true;
end

for k = 1:2
    logMessage(sprintf('Calculating waterfill short and long interpolation zones (it=%d)', k));

    % apply water and no data (ocean) mask to z
    z_masked = z;
    z_masked(M) = nan;
    z_nonwater_nans = isnan(z) & ~M;
    if k == 2
        clear z;
    end

    % Fill the DEM with reference DEM using the Delta method

    % difference map between z_masked and reference DEM with NaNs in voids
    dz = z_masked - R;
    clear z_masked;

    % Force water fill to meet reference surface at a set distance from land.
    interp_short_px = ceil(30 / dx);
    interp_long_px = ceil(1200 / dx);
    long_interp_saved_z_dist_from_shore_px = ceil(100 / dx);
    force_ref_short = ~imdilate(~M, circleMaskSE(interp_short_px));

    long_interp_area = saved_good_z_over_water & imerode(M0, circleMaskSE(long_interp_saved_z_dist_from_shore_px));
    long_interp_area = imdilate(long_interp_area, circleMaskSE(interp_long_px));

    shoreline_area = xor(M, imdilate(M, circleMaskSE(max(3, floor(10/dx)))));
    dz_shoreline = nan(size(M));
    dz_shoreline(shoreline_area) = dz(shoreline_area);
    dz_interp_short = abs(dz_shoreline) > 8.0;
    dz_interp_short = bwareaopen(dz_interp_short, floor(30 / dx^2));
    dz_interp_short_keep = dz_interp_short;

    shoreline_area = xor(M, imdilate(M, ones(3)));
    dz_shoreline = nan(size(M));
    dz_shoreline(shoreline_area) = dz(shoreline_area);
    dz_interp_short = abs(dz_shoreline) >= 8.0;
    dz_interp_short = bwareaopen(dz_interp_short, 3);
    dz_interp_short_keep = dz_interp_short_keep | dz_interp_short;
    long_interp_area(imdilate(dz_interp_short_keep, circleMaskSE(interp_short_px*2))) = 0;

    clear shoreline_area dz_shoreline dz_interp_short dz_interp_short_keep;

    force_ref_zone = force_ref_short;
    force_ref_zone(long_interp_area) = 0;
    dz(force_ref_zone) = 0;
    water_interp_nans = isnan(dz) & M;

    clear force_ref_short long_interp_area force_ref_zone;

%    % Force water fill to meet reference surface at a set distance from land.
%    interp_short_px = ceil(50 / dx);
%    interp_long_px = ceil(1200 / dx);
%    force_ref_short = ~imdilate(~M, circleMaskSE(interp_short_px));
%    force_ref_long = ~imdilate(~M, circleMaskSE(interp_long_px));
%
%    shoreline_area = xor(M, imdilate(M, circleMaskSE(max(3, floor(10/dx)))));
%    dz_shoreline = zeros(size(M));
%    dz_shoreline(shoreline_area) = dz(shoreline_area);
%    dz_interp_short = abs(dz_shoreline) >= 5.0 | isnan(dz_shoreline);
%    dz_interp_short = bwareaopen(dz_interp_short, floor(30 / dx^2));
%    clear shoreline_area dz_shoreline;
%
%    force_ref_zone = force_ref_long;
%    force_ref_zone(force_ref_short & imdilate(dz_interp_short, circleMaskSE(interp_short_px*2))) = 1;
%    clear force_ref_short force_ref_long dz_interp_short;
%
%    dz(force_ref_zone) = 0;
%    water_interp_nans = isnan(dz) & M;
%    clear force_ref_zone;

    if get_M_waterblobs_forced
        logMessage(sprintf('Running bwconncomp to determine fully interpolated waterblobs in M (it=%d)', k));

        waterblobs_fully_interp = zeros(size(M), 'logical');
        CC_1 = bwconncomp(water_interp_nans);
        CC_M = bwconncomp(M);
        mapM = containers.Map(cellfun(@(x) x(1), CC_M.PixelIdxList), 1:CC_M.NumObjects);
        blob_idxlist = [];
        for i = 1:CC_1.NumObjects
            blob_idxlist = CC_1.PixelIdxList{i};
            blob_startidx = blob_idxlist(1);
            if isKey(mapM, blob_startidx) && isequal(blob_idxlist, CC_M.PixelIdxList{mapM(blob_startidx)})
                waterblobs_fully_interp(blob_idxlist) = 1;
            end
        end
        clear CC_1 CC_M mapM blob_idxlist;

        M_waterblobs_forced = M & ~waterblobs_fully_interp;

        clear waterblobs_fully_interp;
    end

    if get_M0_waterblobs_forced
        logMessage('Running bwconncomp to determine fully interpolated waterblobs in M0');

        M0_waterblobs_not_forced = zeros(size(M0), 'logical');
        M0_minus_forced = M0 & ~M_waterblobs_forced;
        CC_M0 = bwconncomp(M0);
        CC_1 = bwconncomp(M0_minus_forced);
        map1 = containers.Map(cellfun(@(x) x(1), CC_1.PixelIdxList), 1:CC_1.NumObjects);
        blob_idxlist = [];
        for i = 1:CC_M0.NumObjects
            blob_idxlist = CC_M0.PixelIdxList{i};
            blob_startidx = blob_idxlist(1);
            if isKey(map1, blob_startidx) && isequal(blob_idxlist, CC_1.PixelIdxList{map1(blob_startidx)})
                M0_waterblobs_not_forced(blob_idxlist) = 1;
            end
        end
        clear CC_M0 CC_1 map1 blob_idxlist;

        M0_waterblobs_forced = M0 & ~M0_waterblobs_not_forced;

        clear M0_waterblobs_not_forced M0_minus_forced;
    end

    M_added = false;

    if k == 1 && add_to_water_mask
        logMessage(sprintf('Checking for pixels below "water height" to add to watermask (it=%d)', k));

        % Define a buffer zone near the edge of water in the water mask.
        M_for_fill = M0_waterblobs_forced;
        water_mask_uncertain_dist_px = ceil(50 / dx);
        water_mask_uncertain_zone = xor(imdilate(M_for_fill, circleMaskSE(water_mask_uncertain_dist_px)), imerode(M_for_fill, circleMaskSE(water_mask_uncertain_dist_px)));

        % Add clusters of land to the water mask where the land is within
        % a short distance from the water's edge and the z height is below
        % both water level and the reference surface.
        fill_water_minarea_px = ceil(3500 / dx^2);
        fill_water_fillhole_px = ceil(500 / dx^2);
        fill_water_keepisland_px = ceil(10000 / dx^2);

        fill_water_check = (M_for_fill | water_mask_uncertain_zone) & z < (water_height - 0.5);
        fill_water_check = bwareaopen(fill_water_check, fill_water_minarea_px);
        fill_water_check = ~bwareaopen(~fill_water_check, fill_water_fillhole_px);

        fill_water = fill_water_check & water_mask_uncertain_zone;
        fill_water_total = fill_water;

        fill_water = fill_water_check & ~water_mask_uncertain_zone;
        fill_water_decide = fill_water & ~M_waterblobs_forced;
        fill_water_decide_dontfill = bwareaopen(fill_water_decide, fill_water_keepisland_px);
        fill_water(fill_water_decide_dontfill) = 0;
        fill_water_total = fill_water_total | fill_water;
        clear fill_water_decide fill_water_decide_dontfill;

        M_to_add = fill_water_total & ~M;
        M_added = any(M_to_add(:));

        if M_added
            M(fill_water_total) = 1;
            M0(fill_water_total) = 1;
        end

        clear M_for_fill M_to_add fill_water_check fill_water fill_water_total water_mask_uncertain_zone;
    end

    if M_added
        logMessage('Added pixels below "water height" to watermask');
    else
        break;
    end
end
clear z saved_good_z_over_water;


logMessage('Splitting arrays into processing tiles for delta fill interplation');

% interpolate differences across void, breaking dz into tiles to speed up
A = tile(dz,5,5,200);
clear dz;

M_tiles = tile(M,5,5,200);

% initiate filled cell array of tiles
B = cell(size(A));

% loop through cells to interpolate voids
i=1;
for i=1:numel(A)
    logMessage(sprintf('Interpolating voids in processing tile %d of %d', i, numel(A)));
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


% Get the height of water from the filled surface
% by sampling the water height a short distance out from land
% and extending that sample back to shore.
logMessage('Getting "water height" from z_filled');
get_water_height_dist_from_shore_px = ceil(100 / dx);
avg_water_height_interp_dist_px = ceil(300 / dx);
get_water_height_mask = imerode(M0, circleMaskSE(get_water_height_dist_from_shore_px));
water_height = nan(size(M0));
water_height(get_water_height_mask) = z_filled(get_water_height_mask);
water_height = movmean(movmean(water_height, avg_water_height_interp_dist_px*2+1, 1, 'omitnan'), avg_water_height_interp_dist_px*2+1, 2, 'omitnan');
% Lower the water height in attempt to save good ArcticDEM data
% in riverbeds and other low-laying areas.
water_height = water_height - 0.15;
water_height = max(water_height, z_filled, 'includenan');
clear get_water_height_mask;

% Flatten to water level divets in the filled DEM caused by interpolation.
underwater_cluster_px = ceil(300 / dx^2);
underwater = M0_waterblobs_forced & z_filled < water_height;
underwater = bwareaopen(underwater, underwater_cluster_px);
if any(underwater(:))
    logMessage('Flattening z_filled below "water height" (it=1)');
    feather_dist_px = ceil(10 / dx);
    underwater = underwater & imerode(M0 & ~isnan(water_height), circleMaskSE(feather_dist_px));
    circle_feather_mask = circleFeatherMask(feather_dist_px, 0);
    flatten_water_weight = conv2(double(underwater), circle_feather_mask, 'same') / sum(circle_feather_mask(:));
    flatten_water_weight(isnan(water_height)) = 0;
    flatten_water_weight = min(1, max(0, flatten_water_weight));
    z_feathered = water_height.*flatten_water_weight + z_filled.*(1 - flatten_water_weight);
    z_feathered(isnan(z_feathered)) = z_filled(isnan(z_feathered));
    z_filled = z_feathered;
    M(flatten_water_weight > 0) = 1;
    clear flatten_water_weight z_feathered;
end
clear underwater;

%underwater = M0_waterblobs_forced & z_filled < (water_height - 0.5) & z_filled < R;
%underwater = bwareaopen(underwater, underwater_cluster_px);
%if any(underwater(:))
%    logMessage('Flattening z_filled below "water height" (it=2)');
%
%    shoreline_area = xor(imdilate(M0, circleMaskSE(max(3, floor(10/dx)))), imerode(M, circleMaskSE(max(3, floor(10/dx)))));
%    dz_shoreline = nan(size(M));
%    dz_shoreline(shoreline_area) = dz_filled(shoreline_area);
%    dz_divots = dz_shoreline < -3;
%    dz_divots = bwareaopen(dz_divots, floor(30 / dx^2));
%    dz_divots_keep = dz_divots;
%
%    shoreline_area = xor(imdilate(M, circleMaskSE(2)), imerode(M, circleMaskSE(2)));
%    dz_shoreline = nan(size(M));
%    dz_shoreline(shoreline_area) = dz_filled(shoreline_area);
%    dz_divots = dz_shoreline < -8;
%    dz_divots_keep = dz_divots_keep | dz_divots;
%
%    clear shoreline_area dz_shoreline dz_divots;
%
%    force_flatten_dist_px = ceil(150 / dx);
%    force_flatten = imdilate(dz_divots_keep, circleMaskSE(force_flatten_dist_px));
%
%    feather_dist_px = ceil(10 / dx);
%    circle_feather_mask = circleFeatherMask(feather_dist_px, 0);
%    flatten_water_weight = conv2(double(force_flatten), circle_feather_mask, 'same') / sum(circle_feather_mask(:));
%    flatten_water_weight(isnan(water_height)) = 0;
%    flatten_water_weight = min(1, max(0, flatten_water_weight));
%
%    flatten_water_weight(~underwater) = 0;
%    z_feathered = water_height.*flatten_water_weight + z_filled.*(1 - flatten_water_weight);
%    z_feathered(isnan(z_feathered)) = z_filled(isnan(z_feathered));
%    z_filled = z_feathered;
%    M(flatten_water_weight > 0) = 1;
%    clear flatten_water_weight z_feathered;
%end
%clear underwater;

%clear water_height;

z_nonwater_nans = z_nonwater_nans & ~M;

% Set back to NaN all original NoData that was outside water mask
% at the beginning of this routine (includes data removed by slope filter).
z_filled(z_nonwater_nans) = nan;

logMessage('Running bwconncomp to remove blobs of waterfill surrounded by NoData');


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
