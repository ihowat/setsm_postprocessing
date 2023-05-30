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

% define water mask using classification and MAD and repeat criteria
M0 = ( C == 80 | C == 0 );

% don't fill within 150m of (at least 300m wide) snow/ice chunks
snowice_mask = (C == 70);
keep_snowice = snowice_mask & imdilate(imerode(snowice_mask, circleMaskSE(ceil(150 / dx))), circleMaskSE(ceil(300 / dx)));
keep_snowice = imdilate(keep_snowice, circleMaskSE(ceil(150 / dx)));
M0 = M0 & (~keep_snowice | isnan(z));
clear snowice_mask keep_snowice;

% only fill where repeat is too low or MAD is high
M = M0 & ( (z_mad > 0.3 & N > 4 ) | N <= 4 | isnan(z) );

% remove small water clusters
cluster_size_px = ceil(750 / dx^2);
M0 = bwareaopen(M0, cluster_size_px);
M = bwareaopen(M, cluster_size_px);

clear C z_mad N;

if ~any(M(:))
    z_filled = z;
    water_mask = zeros(size(z), 'logical');
    dz_filled = [];
    water_height = [];
    M0_waterblobs_forced = [];
    return;
end

% remove small non water clusters
cluster_size_px = ceil(1000 / dx^2);
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


% The following is disabled after we started sourcing water height
% from the post- inpaint_nans 'z_filled'.
if false
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
end


force_ref_zone = [];
M_waterblobs_forced = [];
M0_waterblobs_forced = [];

% At end of weeks of waterfill development for ArcticDEM v4.1 mosaic,
% some large ideas here did not make the final cut and are disabled.
% Use these toggles to reenable further exploration of those ideas.
add_to_water_mask = false;
get_M_waterblobs_forced = false;
get_M0_waterblobs_forced = false;
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

    % Don't force water bodies smaller than 0.5sqkm to meet reference surface.
    small_waterbodies = xor(M, bwareaopen(M, ceil((0.5*1000^2) / dx^2)));
    force_ref_zone(small_waterbodies) = 0;

    dz(force_ref_zone) = 0;
    water_interp_nans = isnan(dz) & M;

    clear force_ref_short long_interp_area small_waterbodies;

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
    if all(M_tiles{i}(:))  % all water, fill with zeros
        B{i} = 0 * ones(size(A{i}),'single');
    elseif ~any(M_tiles{i}(:))  % no water, skip
        B{i} = A{i};
    else  % mix of water and no water, interpolate
%        B{i} = single(inpaint_nans(double(A{i}),2));
        water_border = ~imerode(ones(size(M_tiles{i})), ones(7)) & M_tiles{i};
        water_border(imdilate(~M_tiles{i}, circleMaskSE(floor(300/dx)))) = 0;
        Ai = A{i};
        Ai(water_border) = 0;
        % Had an idea that inpaint_nans could run faster if you tell it which "known" (non-NaN)
        % pixels it should use for NaN interpolation, and have it use just the pixels within 200m
        % of the NaN gaps. It didn't work on the first attempt, and would need more testing.
%        Ai_sampleloc = imdilate(isnan(Ai), circleMaskSE(floor(200/dx)));
%        B{i} = single(inpaint_nans(double(Ai),fillWaterInterpMethod, 'sampleloc',Ai_sampleloc));
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
dz_filled = z_filled;


% Disabled the following and replaced with a modified approach below it
% that essentially does the same thing, but is "topology aware" in that
% it performs this erode-dilate process independently for each blob in the watermask.
if false
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
end

% Get the height of water from the filled surface
% by sampling the water height a short distance out from land
% and extending that sample back to shore.
logMessage('Getting "water height" from z_filled');
get_water_height_dist_px_short = ceil(100 / dx);
avg_water_height_dist_px_short = ceil(300 / dx);
get_water_height_dist_px_long = ceil(200 / dx);
avg_water_height_dist_px_long = ceil(600 / dx);

get_water_height_mask_short = imerode(M0, circleMaskSE(get_water_height_dist_px_short));
get_water_height_mask_long  = imerode(M0, circleMaskSE(get_water_height_dist_px_long));
get_water_height_mask_long = get_water_height_mask_long & imdilate(imerode(get_water_height_mask_long, circleMaskSE(ceil(50 / dx))), circleMaskSE(ceil(100 / dx)));
water_height_short = nan(size(M));
water_height_long  = nan(size(M));
water_height_short(get_water_height_mask_short) = z_filled(get_water_height_mask_short);
water_height_long(get_water_height_mask_long)   = z_filled(get_water_height_mask_long);
%water_height_short(~M) = nan;
%water_height_long(~M)  = nan;

% In glaciated areas, we have seen the ESA WorldCover watermask label as water
% pixels in Copernicus GLO-30 that are heights clearly representing a glacial surface.
% Avoid using those glacial surface heights in sourcing the "water height", because
% that water height data will be smoothed before it is used to perform final water flattening.
water_height_check_block_len = ceil(3000 / dx);
if mod(water_height_check_block_len, 2) == 0
    water_height_check_block_len = water_height_check_block_len + 1;
end
water_height_med = movmedian(movmedian(water_height_short, water_height_check_block_len, 1, 'omitnan'), water_height_check_block_len, 2, 'omitnan');
water_height_mad =       movmad(movmad(water_height_short, water_height_check_block_len, 1, 'omitnan'), water_height_check_block_len, 2, 'omitnan');
ignore_areas = (abs(water_height_short - water_height_med) > 10) | (water_height_mad > 30);
clear water_height_med water_height_mad;

get_water_height_mask_short(ignore_areas) = 0;
get_water_height_mask_long(ignore_areas) = 0;
water_height_short(ignore_areas) = nan;
water_height_long(ignore_areas) = nan;
clear ignore_areas;

water_height = nan(size(M));

if any(get_water_height_mask_short(:))
    CC = bwconncomp(M0);
    blob_idxlist = [];
    for i = 1:CC.NumObjects
        blob_idxlist = CC.PixelIdxList{i};
        if ~any(get_water_height_mask_short(blob_idxlist))
            continue;
        end

    %    water_height_inst = nan(size(M));
    %    water_height_inst(blob_idxlist) = water_height_short(blob_idxlist);

        [idx_row,idx_col] = ind2sub(size(M), blob_idxlist);
        row_min = min(idx_row);
        row_max = max(idx_row);
        col_min = min(idx_col);
        col_max = max(idx_col);
        window_nrows = row_max - row_min + 1;
        window_ncols = col_max - col_min + 1;
        window_size = [window_nrows, window_ncols];
        window_idx = sub2ind(window_size, idx_row-row_min+1, idx_col-col_min+1);

        use_long_height = false;
        % Disabled the "long" water height sourcing because rivers with a width
        % inbetween the short and long sourcing distances, if they are also
        % connected to a large body of water, can end up with gaps in the water height
        % throughout the length of the river.
%        if any(get_water_height_mask_long(blob_idxlist))
%            M_blob = zeros(window_size, 'logical');
%            M_blob(window_idx) = 1;
%            M_blob_long = zeros(window_size, 'logical');
%            M_blob_long(window_idx) = get_water_height_mask_long(blob_idxlist);
%            M_blob_dilated = M_blob & imdilate(M_blob_long, squareMask(avg_water_height_dist_px_long));
%            if (nnz(M_blob_dilated) / nnz(M_blob)) > 0.95
%                use_long_height = true;
%            end
%            clear M_blob M_blob_long M_blob_dilated;
%        end
        if use_long_height
            water_height_source = water_height_long;
            dilate_px = avg_water_height_dist_px_long;
        else
            water_height_source = water_height_short;
            dilate_px = avg_water_height_dist_px_short;
        end

        water_height_inst = nan(window_size);
        water_height_inst(window_idx) = water_height_source(blob_idxlist);

        water_height_inst = movmean(movmean(water_height_inst, dilate_px*2+1, 1, 'omitnan'), dilate_px*2+1, 2, 'omitnan');

        water_height(blob_idxlist) = water_height_inst(window_idx);
    end
    clear CC blob_idxlist;
    clear M_blob water_height_inst;
end

clear get_water_height_mask_short get_water_height_mask_long;
clear water_height_short water_height_long;

% Lower the water height in attempt to save good ArcticDEM data
% in riverbeds and other low-laying areas.
water_height = water_height - 0.15;
water_height = max(water_height, z_filled, 'includenan');


% Flatten to water level divets in the filled DEM caused by interpolation.
underwater_cluster_px = ceil(300 / dx^2);
underwater = M & z_filled < water_height;
underwater = bwareaopen(underwater, underwater_cluster_px);
if any(underwater(:))
    logMessage('Flattening z_filled below "water height" (it=1)');
    feather_dist_px = ceil(10 / dx);
    underwater = underwater & imerode(M & ~isnan(water_height), circleMaskSE(feather_dist_px));
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
