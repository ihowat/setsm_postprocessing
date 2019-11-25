function [projinfo, unit, unitsPerMeter] = getProjInfo(projstr, any_rasterFile)

suppress_warning_id = 'MATLAB:structOnObject';

gdalpath =[]; %set to the path of the gdal binary if not in system path.
if ismac
    gdalpath = '/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/';
elseif ispc
    gdalpath = 'C:/OSGeo4W64/bin/';
end

tempdir = getTempDir();
projtif = fullfile(tempdir, [strrep(projstr, ':', '_'), '.tif']);

if exist(projtif, 'file') ~= 2
    fprintf('Creating dummy raster to access target projection info: %s\n', projtif);
    cmd = sprintf('%s -overwrite -q -ts 1 1 -t_srs "%s" "%s" "%s" ', fullfile(gdalpath, 'gdalwarp'), projstr, any_rasterFile, projtif);
    fprintf('%s\n', cmd);
    [status, cmdout] = system(cmd);
    if ~isempty(cmdout)
        fprintf([cmdout,'\n']);
    end
    if ~exist(projtif, 'file')
        error('gdalwarp call failed: %s', cmd);
    end
end

warning('off', suppress_warning_id)
projinfo = geotiffinfo(projtif);
if isfield(struct(projinfo.SpatialRef), 'AngleUnit')
    unit = projinfo.SpatialRef.AngleUnit;
else
    unit = projinfo.UOMLength;
end
warning('on', suppress_warning_id)

if strcmpi(unit, 'm') || strcmpi(unit, 'meter') || strcmpi(unit, 'meters') || strcmpi(unit, 'metre') || strcmpi(unit, 'metres')
    unitsPerMeter = 1;
elseif strcmpi(unit, 'deg') || strcmpi(unit, 'degree') || strcmpi(unit, 'degrees')
    unitsPerMeter = 0.5/16/3600;
else
    error('Cannot convert target coordinate system unit "%s" to meters', unit);
end