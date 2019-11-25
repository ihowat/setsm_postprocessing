function reproject_database(dbasefile, tilefile)

tempdir = [char(java.lang.System.getProperty('user.home')),'/setsm_postprocessing_temp'];
if exist(tempdir, 'dir') ~= 7
    mkdir(tempdir);
end

tiles = load(tilefile);
epsg = tiles.epsg;
clear tiles;
a = load(dbasefile);

outname = sprintf('%s_%d.mat', strrep(dbasefile, '.mat', ''), epsg);
projstr_fwd = ['EPSG:',num2str(epsg)];
rasterFile = strrep(a.f{1}, 'meta.txt', 'dem.tif');
rasterFile_proj = sprintf('%s_reprojtemp_%d', tempname(tempdir), epsg);

fprintf('Projection set from tilefile: "%s"\n', projstr_fwd);
fprintf('Output database name: %s\n', outname);

fprintf('Creating dummy raster for reference geotiff header info in tiling scheme projection\n');
cmd = sprintf('gdalwarp -overwrite -q -ts 1 1 -t_srs "%s" "%s" "%s" ', projstr_fwd, rasterFile, rasterFile_proj);
[status, cmdout] = system(cmd);
if ~isempty(cmdout)
    fprintf([cmdout,'\n']);
end
if ~exist(rasterFile_proj, 'file')
    error('gdalwarp call failed: %s', cmd);
end
projinfo_fwd = geotiffinfo(rasterFile_proj);
delete(rasterFile_proj);

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
            cmd = sprintf('gdalwarp -overwrite -q -ts 1 1 -t_srs "%s" "%s" "%s" ', projstr_bwd, rasterFile, rasterFile_proj);
            [status, cmdout] = system(cmd);
            if ~isempty(cmdout)
                fprintf([cmdout,'\n']);
            end
            projinfo_bwd = geotiffinfo(rasterFile_proj);
            delete(rasterFile_proj);
        end
    end
    
    if isempty(projinfo_bwd)
        fprintf('(%d/%d) strip projection matches\n', i, num_strips);
    else
        fprintf('(%d/%d) reprojecting Strip Footprint Vertices\n', i, num_strips);
        
        [lat, lon] = projinv(projinfo_bwd, a.x{i}, a.y{i});
        [x, y] = projfwd(projinfo_fwd, lat, lon);
        
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