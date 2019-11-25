function reproject_database(dbasefile, tilefile)

gdalpath =[]; %set to the path of the gdal binary if not in system path.
if ismac
    gdalpath = '/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/';
elseif ispc
    gdalpath = 'C:/OSGeo4W64/bin/';
end

tempdir = getTempDir();

tiles = load(tilefile);
epsg = tiles.epsg;
a = load(dbasefile);

outname = sprintf('%s_%d.mat', strrep(dbasefile, '.mat', ''), epsg);
projstr_fwd = ['EPSG:',num2str(epsg)];

fprintf('Projection set from tilefile: "%s"\n', projstr_fwd);
fprintf('Output database name: %s\n', outname);

random_rasterFile = strrep(a.f{1}, 'meta.txt', 'dem.tif');

[projinfo_fwd, unit, unitsPerMeter] = getProjInfo(projstr_fwd, random_rasterFile);
tiles.unit = unit;
tiles.unitsPerMeter = unitsPerMeter;
save(tilefile,'-struct','tiles','-v7.3');
clear tiles;

projstr_fwd_is_wgs84 = false;
projstr_bwd_is_wgs84 = false;
cmd = sprintf('python proj_issame.py "%s" "%s" ', 'EPSG:4326', projstr_fwd);
[status, cmdout] = system(cmd);
if ~isempty(cmdout)
    fprintf([cmdout,'\n']);
end
if status == 0
    projstr_fwd_is_wgs84 = true;
end
use_gdaltransform = false;

projstr_bwd = '';
projinfo_bwd = [];
num_strips = length(a.projstr);
for i = 1:num_strips
    
    if ~strcmp(a.projstr{i}, projstr_bwd)
        projstr_bwd = a.projstr{i};
        
        cmd = sprintf('python proj_issame.py "%s" "%s" ', projstr_bwd, projstr_fwd);
        [status, cmdout] = system(cmd);
        if ~isempty(cmdout)
            fprintf([cmdout,'\n']);
        end
        if status == 0
            projinfo_bwd = [];
        else
            use_gdaltransform = false;
            projtif_bwd = fullfile(tempdir, [strrep(strrep(strip(projstr_bwd), ':', '_'), ' ', '_'), '.tif']);
            if exist(projtif_bwd, 'file') ~= 2
                cmd = sprintf('%s -overwrite -q -ts 1 1 -t_srs "%s" "%s" "%s" ', fullfile(gdalpath, 'gdalwarp'), projstr_bwd, random_rasterFile, projtif_bwd);
                [status, cmdout] = system(cmd);
                if ~isempty(cmdout)
                    fprintf([cmdout,'\n']);
                end
            end
            projinfo_bwd = geotiffinfo(projtif_bwd);
            cmd = sprintf('python proj_issame.py "%s" "%s" ', 'EPSG:4326', projstr_bwd);
            [status, cmdout] = system(cmd);
            if ~isempty(cmdout)
                fprintf([cmdout,'\n']);
            end
            if status == 0
                projstr_bwd_is_wgs84 = true;
            else
                projstr_bwd_is_wgs84 = false;
            end
%             projinfo_bwd = 1;
        end
    end
    
    if isempty(projinfo_bwd)
        fprintf('(%d/%d) strip projection matches\n', i, num_strips);
    else
        fprintf('(%d/%d) reprojecting Strip Footprint Vertices\n', i, num_strips);
        
        if ~use_gdaltransform
            try
                if projstr_bwd_is_wgs84
                    lon = a.x{i};
                    lat = a.y{i};
                else
                    [lat, lon] = projinv(projinfo_bwd, a.x{i}, a.y{i});
                end
                if projstr_fwd_is_wgs84
                    x = lon;
                    y = lat;
                else
                    [x, y] = projfwd(projinfo_fwd, lat, lon);
                end
            catch ME
                fprintf('%s\n', ME.identifier);
                fprintf('%s\n', ME.message);
                fprintf('Switching to gdaltransform method for handling strip projection: %s\n', projstr_fwd);
                use_gdaltransform = true;
            end
        end
        
        if use_gdaltransform
            coords_in = strrep(mat2str([a.x{i}(:), a.y{i}(:)]), ';', '\n');
            coords_in = coords_in(2:end-1);
            cmd = sprintf('echo -e "%s" | %s -s_srs "%s" -t_srs "%s" -output_xy', coords_in, fullfile(gdalpath, 'gdaltransform'), projstr_bwd, projstr_fwd);
            [status, cmdout] = system(cmd);
            if status ~= 0
                if ~isempty(cmdout)
                    fprintf([cmdout,'\n']);
                end
                error('Failed to reproject strip footprint vertices using gdaltransform');
            end
            coords_out = str2num(cmdout);
            x = coords_out(:, 1)';
            y = coords_out(:, 2)';
        end
        
        a.projstr{i} = projstr_fwd;
        
        a.x{i} = x;
        a.xmin(i) = min(x);
        a.xmax(i) = max(x);
        a.y{i} = y;
        a.ymin(i) = min(y);
        a.ymax(i) = max(y);
        
    end
    
end

save(outname,'-struct','a','-v7.3');
fprintf('%s saved\n', outname);

end