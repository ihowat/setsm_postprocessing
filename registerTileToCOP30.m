function [z,z_at_zr_res,I_ref,C] = registerTileToCOP30(matfile_or_z,x,y,refDemTif,coverTif,varargin)

water_buffer_meters = 200;
ice_buffer_meters = 100;
max_offset_meters = 10;
ref_stdev_thresh = 1.5;

n = find(strcmpi('registerBlobs',varargin));
if ~isempty(n)
    registerBlobs = true;
else
    registerBlobs = false;
end

n = find(strcmpi('registerBlobsSkipregShp',varargin));
if ~isempty(n)
    registerBlobsSkipregShp = varargin{n+1};
else
    registerBlobsSkipregShp = [];
end

n = find(strcmpi('reportOffsetOnly',varargin));
if ~isempty(n)
    reportOffsetOnly = true;
else
    reportOffsetOnly = false;
end

n = find(strcmpi('debugTifProjstr',varargin));
if ~isempty(n)
    debugTifProjstr = varargin{n+1};
else
    debugTifProjstr = 'polar stereo north';
end

if ~isempty(registerBlobsSkipregShp)
    fprintf('Loading registerBlobs "skipreg" shapefile: %s\n', registerBlobsSkipregShp);
    skipreg_mapstruct = shaperead(registerBlobsSkipregShp);
    skipreg_polyshape_arr = arrayfun(@(feat) polyshape(feat.X, feat.Y), skipreg_mapstruct);
else
    skipreg_mapstruct = [];
    skipreg_polyshape_arr = [];
end

fprintf('Loading reference DEM and ESA WorldCover raster\n');
[z,x,y,z_at_zr_res,I_ref,C] = loadSlopefiltWaterfillArrays(matfile_or_z,x,y,refDemTif,coverTif);
zr_dx = I_ref.x(2)-I_ref.x(1);

if ~reportOffsetOnly
    z_at_zr_res_unreg = z_at_zr_res;
else
    z_at_zr_res_unreg = [];
end

z_masked = z_at_zr_res;
zr_masked = I_ref.z;

% Apply water mask
M = (C.z == 0 | C.z == 80);
uncertain_dist_px = ceil(water_buffer_meters / zr_dx);
M = imdilate(M, circleMask(uncertain_dist_px));
z_masked(M) = nan;
zr_masked(M) = nan;
watermask_dilated = M;
clear M;

% Apply snow/ice mask
M = (C.z == 70);
uncertain_dist_px = ceil(ice_buffer_meters / zr_dx);
M = imdilate(M, circleMask(uncertain_dist_px));
z_masked(M) = nan;
zr_masked(M) = nan;
clear M;

% Apply slope filter
fprintf('Calculating and applying slope filter\n');
avoidFilteringWaterFlag = true;
M = slopeDifferenceFilter(I_ref.x,I_ref.y,z_at_zr_res,I_ref.z,C,avoidFilteringWaterFlag);
z_masked(~M) = nan;
clear M;

% Check if any islands from watermask are actually flat filled (water) surface
% in the reference DEM, and set them to NaN so they won't be used in offset calc.
CC = bwconncomp(~isnan(zr_masked));
for i = 1:CC.NumObjects
    blob_idxlist = CC.PixelIdxList{i};
    zr_values = zr_masked(blob_idxlist);
    if std(zr_values(:), 'omitnan') < ref_stdev_thresh
        zr_masked(blob_idxlist) = nan;
    end
end
clear CC;

tile_avg_offset = median(z_masked - zr_masked, 'all', 'omitnan');
fprintf('Median whole-tile offset with water-masked reference DEM: %g\n', tile_avg_offset);
if ~registerBlobs
    if isnan(tile_avg_offset)
        fprintf('Adjusting offset from NaN to zero\n');
        tile_avg_offset = 0;
    elseif abs(tile_avg_offset) > max_offset_meters
        fprintf('Adjusting absolute offset > %d meters to zero\n', max_offset_meters);
        tile_avg_offset = 0;
    elseif std(zr_masked(:), 'omitnan') < ref_stdev_thresh
        fprintf('Adjusting offset to zero because reference DEM stdev < %d\n', ref_stdev_thresh);
        tile_avg_offset = 0;
    end
end

avg_offset_array_zr_res = [];

if registerBlobs && ~reportOffsetOnly
    % "blobs" are contiguous chunks of the 1x1km subtiles that were processed.
    z_data = ~isnan(z_at_zr_res);
    z_data_land = z_data;
    z_data_land(watermask_dilated) = 0;

    z_data_copy = z_data;

    CC = bwconncomp(z_data);
%    if CC.NumObjects == 1
%        registerBlobs = false;
    if true
        avg_offset_array = nan(size(z_at_zr_res));
        for i = 1:CC.NumObjects
            blob_idxlist = CC.PixelIdxList{i};
            skipreg = false;

            if ~isempty(registerBlobsSkipregShp)
                [idx_row,idx_col] = ind2sub(size(z_data), blob_idxlist);
                B = bwtraceboundary(z_data, [idx_row(1) idx_col(1)], 'E');
                B_row = B(:,1);
                B_col = B(:,2);
                blob_poly = polyshape(I_ref.x(B_col), I_ref.y(B_row));
                skipreg_polys_overlapped = overlaps(blob_poly, skipreg_polyshape_arr);
                skipreg_polys_overlapped_ms = skipreg_mapstruct(skipreg_polys_overlapped);
                if ~isempty(skipreg_polys_overlapped_ms)
                    blob_M = [];
                    for skipreg_poly_ms_i = 1:length(skipreg_polys_overlapped_ms)
                        skipreg_poly_ms = skipreg_polys_overlapped_ms(skipreg_poly_ms_i);
                        if skipreg_poly_ms.skipreg == 1
                            if isempty(blob_M)
                                blob_M = zeros(size(z_data), 'logical');
                                blob_M(blob_idxlist) = 1;
                            end
                            poly_x = skipreg_poly_ms.X;
                            poly_y = skipreg_poly_ms.Y;
                            poly_nan = find(isnan(poly_x) | isnan(poly_y));
                            poly_x(poly_nan) = [];
                            poly_y(poly_nan) = [];
                            skipreg_poly_M = roipoly(I_ref.x,I_ref.y,z_data_copy,poly_x,poly_y);
                            blob_skipreg_match = blob_M & skipreg_poly_M;
                            if any(blob_skipreg_match(:))
                                skipreg = true;
                                break;
                            end
                        end
                    end
                end
            end

            if skipreg
                avg_offset = 0;
            else
                z_values = z_masked(blob_idxlist);
                zr_values = zr_masked(blob_idxlist);
                land_data = z_data_land(blob_idxlist);
                if std(zr_values, 'omitnan') < ref_stdev_thresh
                    % Assume we are looking at a flat filled (water) surface in the reference DEM
                    % where there may be good data in z, so make this blob offset zero.
                    avg_offset = 0;
                else
                    dz_values = z_values - zr_values;
                    if (nnz(~isnan(dz_values)) / nnz(land_data)) < 0.25
                        % Don't use offset where less than 25% of blob's land area is used for calc
                        avg_offset = 0;
                    elseif (nnz(~isnan(dz_values)) * zr_dx^2) < 1000^2
                        % Don't use offset where fewer than 1 sqkm of pixels are used for calc
                        avg_offset = 0;
                    else
                        avg_offset = median(z_values - zr_values, 'all', 'omitnan');
                        if isnan(avg_offset) || abs(avg_offset) > max_offset_meters
                            avg_offset = 0;
                        end
                    end
                end
            end
            avg_offset_array(blob_idxlist) = avg_offset;
        end
        clear z_data CC blob_idxlist z_values zr_values;

        extend_dist_px = ceil(50 / zr_dx);
        avg_offset_array = movmean(movmean(avg_offset_array, extend_dist_px*2+1, 1, 'omitnan'), extend_dist_px*2+1, 2, 'omitnan');
        
        avg_offset_array_zr_res = avg_offset_array;
        avg_offset_array_z_res = interp2(I_ref.x,I_ref.y(:),avg_offset_array,x,y(:),'*nearest');

        % Apply offset to returned z
        fprintf('Applying per-blob offset to z\n');
        z_orig_nans = isnan(z);
        z = z - avg_offset_array_z_res;
        z_reg_nans = isnan(z);
        if ~isequal(z_orig_nans, z_reg_nans)
            error('Developer error: Failed to keep NaN areas unchanged through registerBlobs');
        end
        clear avg_offset_array avg_offset_array_z_res z_orig_nans z_reg_nans;
    end
    clear z_data CC;
end

clear z_masked zr_masked;

if ~(reportOffsetOnly || registerBlobs)
    % Apply offset to returned z
    fprintf('Applying whole-tile offset to z\n');
    z = z - tile_avg_offset;
end

if ~reportOffsetOnly
    % Remake returned z_at_zr_res from adjusted z
    [~,~,~,z_at_zr_res,~,~] = loadSlopefiltWaterfillArrays(z,x,y,refDemTif,[],'skipRefDemLoad');
end
