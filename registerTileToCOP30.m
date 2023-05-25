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
%    registerBlobsSkipregShp = [];
    registerBlobsSkipregShp = "/mnt/pgc/data/elev/dem/setsm/ArcticDEM/mosaic/v4.1/results/output_tiles/arcticdem_may19_mosaic_100m_v4.1_dem_debug-reg_blobs_offset_gt6.shp";
end

n = find(strcmpi('reportOffsetOnly',varargin));
if ~isempty(n)
    reportOffsetOnly = true;
else
    reportOffsetOnly = false;
end

n = find(strcmpi('writeDebugTifs',varargin));
if ~isempty(n)
    writeDebugTifs = true;
else
%    writeDebugTifs = false;
    writeDebugTifs = true;
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

if writeDebugTifs && ~reportOffsetOnly
    z_at_zr_res_unreg = z_at_zr_res;
else
    z_at_zr_res_unreg;
end

z_masked = z_at_zr_res;
zr_masked = I_ref.z;

% Apply water mask
M = (C.z == 80);
watermask_at_zr_res = M;
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
M = slopeDifferenceFilter(I_ref.x,I_ref.y,z_at_zr_res,I_ref.z,watermask_at_zr_res);
z_masked(~M) = nan;
clear M;

clear watermask_at_zr_res;

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

if writeDebugTifs
    avg_offset_array_zr_res_tile = tile_avg_offset * ones(size(z_at_zr_res));
    if isempty(avg_offset_array_zr_res)
        avg_offset_array_zr_res = avg_offset_array_zr_res_tile;
    end
end

if ~reportOffsetOnly
    % Remake returned z_at_zr_res from adjusted z
    [~,~,~,z_at_zr_res,~,~] = loadSlopefiltWaterfillArrays(z,x,y,refDemTif,[],'skipRefDemLoad');
end

%% Write out debug tifs
if writeDebugTifs && (isa(matfile_or_z, 'char') || isa(matfile_or_z, 'string'))
    outNameBase = matfile_or_z;
    outNameBase = strrep(outNameBase,'_reg.mat','.mat');
    outNameBase = strrep(outNameBase,'_unreg.mat.bak','.mat');

    tif_format = 'COG';
%    tif_format = 'GTiff';
    projstr = debugTifProjstr;
    co_predictor = '3';
%    co_predictor = '1';
    
    db_dx = 32;

    outNameBase = strrep(outNameBase,'_2m.mat',sprintf('_%dm.mat', db_dx));
    
    if db_dx == zr_dx
        db_x = I_ref.x;
        db_y = I_ref.y;
        z_at_db_res_unreg = z_at_zr_res_unreg;
        z_at_db_res = z_at_zr_res;
        zr_at_db_res = I_ref.z;
        avg_offset_array_db_res = avg_offset_array_zr_res;
        avg_offset_array_db_res_tile = avg_offset_array_zr_res_tile;
    else
        db_x0 = floor(I_ref.x(1)  /db_dx) * db_dx;
        db_x1 = ceil( I_ref.x(end)/db_dx) * db_dx;
        db_y0 = floor(I_ref.y(end)/db_dx) * db_dx;
        db_y1 = ceil( I_ref.y(1)  /db_dx) * db_dx;
        db_x = db_x0:db_dx:db_x1;
        db_y = db_y1:-db_dx:db_y0;
        
        if ~isempty(z_at_zr_res_unreg)
            z_at_db_res_unreg = interp2(I_ref.x,I_ref.y(:),z_at_zr_res_unreg,db_x,db_y(:),'*bilinear');
        end
        z_at_db_res = interp2(I_ref.x,I_ref.y(:),z_at_zr_res,db_x,db_y(:),'*bilinear');
        zr_at_db_res = interp2(I_ref.x,I_ref.y(:),I_ref.z,db_x,db_y(:),'*bilinear');
        if ~isempty(avg_offset_array_zr_res)
            avg_offset_array_db_res = interp2(I_ref.x,I_ref.y(:),avg_offset_array_zr_res,db_x,db_y(:),'*bilinear');
        end
        if ~isempty(avg_offset_array_zr_res_tile)
            avg_offset_array_db_res_tile = interp2(I_ref.x,I_ref.y(:),avg_offset_array_zr_res_tile,db_x,db_y(:),'*bilinear');
        end
    end

    if reportOffsetOnly
        debug_tif_suffix = 'reportOffsetOnly';
    elseif registerBlobs
        debug_tif_suffix = 'blobs';
    else
        debug_tif_suffix = 'tiles';
    end

    suffix_list = {debug_tif_suffix};
    avg_offset_list = {avg_offset_array_db_res};
    z_at_db_res_list = {z_at_db_res};
%    if ~reportOffsetOnly
%        suffix_list{end+1} = 'reportOffsetOnly';
%        avg_offset_list{end+1} = avg_offset_array_db_res_tile;
%        z_at_db_res_list{end+1} = z_at_db_res_unreg;
%    end

    for i = 1:length(suffix_list)
        debug_tif_suffix = suffix_list{i};
        avg_offset_array = avg_offset_list{i};
        z_at_db_res = z_at_db_res_list{i};

        if ~isempty(avg_offset_array)
            if strcmp(debug_tif_suffix, 'reportOffsetOnly')
                offset_suffix = 'tiles';
            else
                offset_suffix = debug_tif_suffix;
            end
            outNameTif = strrep(outNameBase, '.mat', sprintf('_dem_debug-reg_%s_offset.tif', offset_suffix));
            fprintf('Writing %s\n', outNameTif)
            if exist(outNameTif,'file')
                delete(outNameTif);
            end
            if true
                % Round float values to 1/128 meters to greatly improve compression effectiveness
                avg_offset_array=round(avg_offset_array*128.0)/128.0;
                avg_offset_array(isnan(avg_offset_array)) = -9999;
                writeGeotiff(outNameTif,db_x,db_y,avg_offset_array,4,-9999,projstr,'out_format',tif_format,'co_predictor',co_predictor,'cog_overview_resampling','BILINEAR')
                clear avg_offset_array;
            end
        end

        if ~isempty(z_at_db_res)

            outNameTif = strrep(outNameBase, '.mat', sprintf('_dem_debug-reg_%s.tif', debug_tif_suffix));
            fprintf('Writing %s\n', outNameTif)
            if exist(outNameTif,'file')
                delete(outNameTif);
            end
            if true
                % Round float values to 1/128 meters to greatly improve compression effectiveness
                z_at_db_res_inst=round(z_at_db_res*128.0)/128.0;
                z_at_db_res_inst(isnan(z_at_db_res_inst)) = -9999;
                writeGeotiff(outNameTif,db_x,db_y,z_at_db_res_inst,4,-9999,projstr,'out_format',tif_format,'co_predictor',co_predictor,'cog_overview_resampling','BILINEAR')
                clear z_at_db_res_inst;
            end

            dz = z_at_db_res - zr_at_db_res;

            outNameTif = strrep(outNameBase, '.mat', sprintf('_demdiff_debug-reg_%s.tif', debug_tif_suffix));
            fprintf('Writing %s\n', outNameTif)
            if exist(outNameTif,'file')
                delete(outNameTif);
            end
            if true
                % Round float values to 1/128 meters to greatly improve compression effectiveness
                dz=round(dz*128.0)/128.0;
                dz(isnan(dz)) = -9999;
                writeGeotiff(outNameTif,db_x,db_y,dz,4,-9999,projstr,'out_format',tif_format,'co_predictor',co_predictor,'cog_overview_resampling','BILINEAR')
            end
            clear dz;

        end
    end
end
