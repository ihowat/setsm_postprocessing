function writeTileToTifv4(tilef,projstr)
% Write 2m or 10m dem matfiles to tif
%   compatible with setsm_postprocesing v4 branch

fprintf('Source: %s\n',tilef);

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

%crop coordinate vectors
x=x(nx(1):nx(2));
y=y(ny(1):ny(2));

fprintf('Writing DEM\n')
outNameDem = strrep(tilef,'.mat','_dem.tif');
if exist(outNameDem,'file')
    fprintf('%s exists, skipping\n',outNameDem);
else
    z=m.z(ny(1):ny(2),nx(1):nx(2));
    z(isnan(z)) = -9999;
    writeGeotiff(outNameDem,x,y,z,4,-9999,projstr)
    clear z
end

flds=fields(m);

if contains('z_mad',flds)
    fprintf('Writing mad\n')
    outNameTif = strrep(tilef,'.mat','_mad.tif');
    if exist(outNameTif,'file')
        fprintf('%s exists, skipping\n',outNameTif);
    else
        z_mad=m.z_mad(ny(1):ny(end),nx(1):nx(end));
        z_mad(isnan(z_mad)) = -9999;
        writeGeotiff(outNameTif,x,y,z_mad,4,-9999,projstr)
        clear z_mad
    end
end

% data count
if contains('N',flds)
    fprintf('Writing N\n')
    outNameTif = strrep(tilef,'.mat','_count.tif');
    if exist(outNameTif,'file')
        fprintf('%s exists, skipping\n',outNameTif);
    else
        N=m.N(ny(1):ny(end),nx(1):nx(end));
        writeGeotiff(outNameTif,x,y,N,1,0,projstr)
        clear N
    end
end

% matchtag count
if contains('Nmt',flds)
    fprintf('Writing Nmt\n')
    outNameTif = strrep(tilef,'.mat','_countmt.tif');
    if exist(outNameTif,'file')
        fprintf('%s exists, skipping\n',outNameTif);
    else
        Nmt=m.Nmt(ny(1):ny(end),nx(1):nx(end));
        writeGeotiff(outNameTif,x,y,Nmt,1,0,projstr)
        clear Nmt
    end
end

% Maximum date
if contains('tmax',flds)
    fprintf('Writing tmax\n')
    outNameTif = strrep(tilef,'.mat','_maxdate.tif');
    if exist(outNameTif,'file')
        fprintf('%s exists, skipping\n',outNameTif);
    else
        tmax=m.tmax(ny(1):ny(end),nx(1):nx(end));
        writeGeotiff(outNameTif,x,y,tmax,2,0,projstr)
        clear tmax
    end
end
    
    % Minimum date
if contains('tmin',flds)
    fprintf('Writing tmin\n')
    outNameTif = strrep(tilef,'.mat','_mindate.tif');
    if exist(outNameTif,'file')
        fprintf('%s exists, skipping\n',outNameTif);
    else
        tmin=m.tmin(ny(1):ny(end),nx(1):nx(end));
       writeGeotiff(outNameTif,x,y,tmin,2,0,projstr)
        clear tmin
    end
end


%% 10m browse hillshade
fprintf('Writing browse\n')

% if 2m posting, first downsample to 10m
outNameBrowse = strrep(tilef,'.mat','_browse.tif');
if exist(outNameBrowse,'file')
    fprintf('%s exists, skipping\n',outNameBrowse);
else
    if dx == 2
        outNameTemp = strrep(tilef,'.mat','_temp.tif');
        system(['gdal_translate -q -tr 10 10 -r bilinear -co bigtiff=if_safer -a_nodata -9999 ',...
            outNameDem,' ', outNameTemp]);
    elseif dx == 10
        outNameTemp = outNameDem;
    end

    % convert to hillshade
    system(['gdaldem hillshade -q -z 3 -compute_edges -of GTiff -co TILED=YES -co BIGTIFF=IF_SAFER -co COMPRESS=LZW ',...
        outNameTemp,' ', outNameBrowse]);
    
    if dx == 2
        delete(outNameTemp);
    end
end






