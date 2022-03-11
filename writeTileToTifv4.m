function writeTileToTifv4(tilef,projstr,varargin)
% Write 2m or 10m dem matfiles to tif
%   compatible with setsm_postprocesing v4 branch

n = find(strcmpi('outRasterType',varargin));
if ~isempty(n)
    outRasterType = varargin{n+1};
else
    outRasterType = 'GTiff';
end
outRasterType_choices = {'browse', 'GTiff', 'COG'};
if ~any(strcmp(outRasterType, outRasterType_choices))
    error("'outRasterType' must be one of the following: {'%s'}", strjoin(outRasterType_choices, "', '"))
end

if strcmp(outRasterType, 'COG')
    outRasterType_is_cog = true;
    tif_format = 'COG';
else
    outRasterType_is_cog = false;
    tif_format = 'GTiff';
end

fprintf('Source: %s\n',tilef);

% output tile buffer in pixels
buffer = 10;

% load m file and get coordinate vectors
m=matfile(tilef);
x=m.x;
y=m.y;

% find data/buffer boundaries
% get posting distance
dx = x(2)-x(1);

%find buffer size
if dx == 2 % 2m posting, using quarter tile (50km) boundaries
    nx = find(mod(x,50000) == 0);
    ny = find(mod(y,50000) == 0);
elseif dx == 10 % 10m posting, using full tile (100km) boundaries
    nx = find(mod(x,100000) == 0);
    ny = find(mod(y,100000) == 0);
    % EarthDEM UTM mosaic 100km tile edges fall on 50km intervals
    if length(nx) == 1
        nx = find(mod(x,50000) == 0);
        if length(nx) == 3
            nx = [nx(1) nx(end)];
        end
        ny = find(mod(y,50000) == 0);
        if length(ny) == 3
            ny = [ny(1) ny(end)];
        end
    end
else
    error('not compatible with a tile grid size of %dm',dx)
end

% if no index values returned, assume all tile data lies within tile
% boundary - no buffers;
if isempty(nx)
    nx = [1 length(x)+1];
end
if isempty(ny)
    ny = [1 length(y)+1];
end

% if one index value returned, tile does not extend to one boundary.
% Determine which boundary this is by which side its closer to. This may
% fail if the buffer width is larger than the width of data inside tile.
if length(nx) == 1
    % if boundary point is closer to y(1), assume buffer is 1:ny-1 and set
    % data range to the length of the x coordinate , plus one for the crop.
    if nx/length(x) < 0.5
        nx(2) = length(x)+1;
    else
        % if boundary point is closer to y(end), assume buffer is ny-1:end
        % and set minimum range to 1
        nx = [1,nx];
    end
end
if length(ny) == 1
    % if boundary point is closer to y(1), assume buffer is 1:ny-1 
    if ny/length(y) < 0.5
        ny(2) = length(y)+1;
    else
        % if boundary point is closer to y(end), assume buffer is ny-1:end
        ny = [1,ny];
    end
end

% crop last row & col off, since tile begins at x0 and y0, so no overlap
% with next tile
nx(end) = nx(end)-1;
ny(end) = ny(end)-1;

% add standard tile buffer
if buffer > 0
    nx(1) = max(nx(1)-buffer, 1);
    nx(2) = min(nx(2)+buffer, length(x));
    ny(1) = max(ny(1)-buffer, 1);
    ny(2) = min(ny(2)+buffer, length(y));
end

%crop coordinate vectors
x=x(nx(1):nx(2));
y=y(ny(1):ny(2));

%get name without "_reg"
outNameBase = strrep(tilef,'_reg.mat','.mat');


fprintf('Writing DEM\n')
outNameDem = strrep(outNameBase,'.mat','_dem.tif');
if exist(outNameDem,'file')
    fprintf('%s exists, skipping\n',outNameDem);
else
    z=m.z(ny(1):ny(end),nx(1):nx(end));
    % Round DEM values to 1/128 meters to greatly improve compression effectiveness
    z=round(z*128.0)/128.0;
    z(isnan(z)) = -9999;
    if outRasterType_is_cog
        co_predictor = 'FLOATING_POINT';
    else
        co_predictor = '3';
    end
    writeGeotiff(outNameDem,x,y,z,4,-9999,projstr,'out_format',tif_format,'co_predictor',co_predictor,'cog_overview_resampling','BILINEAR')
    clear z
end

flds=fields(m);

if ~strcmp(outRasterType, 'browse')
    if contains('z_mad',flds)
        fprintf('Writing mad\n')
        outNameTif = strrep(outNameBase,'.mat','_mad.tif');
        if exist(outNameTif,'file')
            fprintf('%s exists, skipping\n',outNameTif);
        else
            z_mad=m.z_mad(ny(1):ny(end),nx(1):nx(end));
            % Round MAD values to 1/128 meters to greatly improve compression effectiveness
            z_mad=round(z_mad*128.0)/128.0;
            z_mad(isnan(z_mad)) = -9999;
            if outRasterType_is_cog
                co_predictor = 'FLOATING_POINT';
            else
                co_predictor = '3';
            end
            writeGeotiff(outNameTif,x,y,z_mad,4,-9999,projstr,'out_format',tif_format,'co_predictor',co_predictor,'cog_overview_resampling','BILINEAR')
            clear z_mad
        end
    end

    % data count
    if contains('N',flds)
        fprintf('Writing N\n')
        outNameTif = strrep(outNameBase,'.mat','_count.tif');
        if exist(outNameTif,'file')
            fprintf('%s exists, skipping\n',outNameTif);
        else
            N=m.N(ny(1):ny(end),nx(1):nx(end));
            if outRasterType_is_cog
                co_predictor = 'NO';
%                co_predictor = 'STANDARD';
            else
                co_predictor = '1';
%                co_predictor = '2';
            end
            writeGeotiff(outNameTif,x,y,N,1,0,projstr,'out_format',tif_format,'co_predictor',co_predictor,'cog_overview_resampling','NEAREST')
            clear N
        end
    end

    % matchtag count
    if contains('Nmt',flds)
        fprintf('Writing Nmt\n')
        outNameTif = strrep(outNameBase,'.mat','_countmt.tif');
        if exist(outNameTif,'file')
            fprintf('%s exists, skipping\n',outNameTif);
        else
            Nmt=m.Nmt(ny(1):ny(end),nx(1):nx(end));
            if outRasterType_is_cog
                co_predictor = 'NO';
%                co_predictor = 'STANDARD';
            else
                co_predictor = '1';
%                co_predictor = '2';
            end
            writeGeotiff(outNameTif,x,y,Nmt,1,0,projstr,'out_format',tif_format,'co_predictor',co_predictor,'cog_overview_resampling','NEAREST')
            clear Nmt
        end
    end

    % Maximum date
    if contains('tmax',flds)
        fprintf('Writing tmax\n')
        outNameTif = strrep(outNameBase,'.mat','_maxdate.tif');
        if exist(outNameTif,'file')
            fprintf('%s exists, skipping\n',outNameTif);
        else
            tmax=m.tmax(ny(1):ny(end),nx(1):nx(end));
            if outRasterType_is_cog
                co_predictor = 'NO';
%                co_predictor = 'STANDARD';
            else
                co_predictor = '1';
%                co_predictor = '2';
            end
            writeGeotiff(outNameTif,x,y,tmax,2,0,projstr,'out_format',tif_format,'co_predictor',co_predictor,'cog_overview_resampling','NEAREST')
            clear tmax
        end
    end

    % Minimum date
    if contains('tmin',flds)
        fprintf('Writing tmin\n')
        outNameTif = strrep(outNameBase,'.mat','_mindate.tif');
        if exist(outNameTif,'file')
            fprintf('%s exists, skipping\n',outNameTif);
        else
            tmin=m.tmin(ny(1):ny(end),nx(1):nx(end));
            if outRasterType_is_cog
                co_predictor = 'NO';
%                co_predictor = 'STANDARD';
            else
                co_predictor = '1';
%                co_predictor = '2';
            end
            writeGeotiff(outNameTif,x,y,tmin,2,0,projstr,'out_format',tif_format,'co_predictor',co_predictor,'cog_overview_resampling','NEAREST')
            clear tmin
        end
    end
end


%% 10m browse hillshade
fprintf('Writing browse\n')

% if 2m posting, first downsample to 10m
outNameBrowse = strrep(outNameBase,'.mat','_browse.tif');
if exist(outNameBrowse,'file')
    fprintf('%s exists, skipping\n',outNameBrowse);
else
    [epsg,~,~] = getProjstrInfo(projstr);

    if dx == 2
        outNameTemp = strrep(tilef,'.mat','_temp.tif');
        gdal_cmd = [...
            'gdal_translate',...
            ' -of GTiff',...
            ' -a_srs ', sprintf('EPSG:%d', epsg),...
            ' -tr 10 10 -r bilinear',...
            ' -co TILED=YES -co BIGTIFF=YES',...
            ' -co COMPRESS=LZW -co PREDICTOR=3',...
            ' -a_nodata -9999',...
            ' ', outNameDem, ' ', outNameTemp];
        fprintf('%s\n', gdal_cmd);
        [status, cmdout] = system(gdal_cmd);
        if ~isempty(cmdout)
            disp(cmdout)
        end
        if status ~= 0
            error('Non-zero exit status (%d) from gdal_translate',status)
        end
    elseif dx == 10
        outNameTemp = outNameDem;
    end

    if outRasterType_is_cog
        co_predictor = 'STANDARD';
        co_overviews_args = '-co OVERVIEWS=IGNORE_EXISTING -co RESAMPLING=CUBIC';
    else
        co_predictor = '2';
        co_overviews_args = '';
    end

    % convert to hillshade
    gdal_cmd = [...
        'gdaldem hillshade',...
        ' -of ', tif_format,...
        ' -z 3 -compute_edges',...
        ' -co TILED=YES -co BIGTIFF=YES',...
        ' -co COMPRESS=LZW -co PREDICTOR=', co_predictor,...
        ' ', co_overviews_args,...
        ' ', outNameTemp, ' ', outNameBrowse];
    fprintf('%s\n', gdal_cmd);
    [status, cmdout] = system(gdal_cmd);
    if ~isempty(cmdout)
        disp(cmdout)
    end
    if status ~= 0
        error('Non-zero exit status (%d) from gdaldem hillshade',status)
    end

    if dx == 2
        delete(outNameTemp);
    end
end


if strcmp(outRasterType, 'browse')
    delete(outNameDem)
end







