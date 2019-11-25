function I= readGeotiff(rasterFile,varargin)
% readGeotiff: read geotiff using imread and assign map info from infinfo.
%
% I = readGeotiff('filename'); reads entire .tif into map structure I
% I = GEOTIFF_READ(...,'pixel_subset', [r0 r1 c0 c1]) reads the subset of
% the image defined by row and column range.
% I = GEOTIFF_READ(...,'map_subset', [x0 x1 y0 y1]); reads the subset of
% the image defined by the map coordinate range.
% I = GEOTIFF_READ(...,'target_projstring', PROJSTRING); projects the image
% into the coordinate system specified by PROJSTRING (may be a PROJ4 string
% or an 'EPSG:XXXX' EPSG string) before reading it in. 
% output:
% I.z, image data
% I.x, x coordinate in map
% I.y, y coordinate in map
% I.info, misc. info
% imshow(I.z, 'xdata', I.x, 'ydata', I.y);
% shows image with map coordinate
% Version by Yushin Ahn, ahn.74@osu.edu
% Glacier Dynamics Laboratory, 
% Byrd Polar Resear Center, Ohio State University 
% Referenced enviread.m (Ian Howat)

gdalpath =[]; %set to the path of the gdal binary if not in system path.
if ismac
    gdalpath = '/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/';
elseif ispc
    gdalpath = 'C:/OSGeo4W64/bin/';
end

Tinfo = imfinfo(rasterFile);
if size(Tinfo, 1) > 1
    Tinfo = Tinfo(1);
end
info.cols   = Tinfo.Width;
info.rows   = Tinfo.Height;
info.imsize = Tinfo.Offset;
info.bands  = Tinfo.SamplesPerPixel;

info.map_info.dx    = Tinfo.ModelPixelScaleTag(1);
info.map_info.dy    = Tinfo.ModelPixelScaleTag(2);
info.map_info.mapx  = Tinfo.ModelTiepointTag(4);
info.map_info.mapy  = Tinfo.ModelTiepointTag(5);

subrows = [1 info.rows];
subcols = [1 info.cols];

minx = info.map_info.mapx;
maxy = info.map_info.mapy;

x = minx + ((0:info.cols-1).*info.map_info.dx);
y = maxy - ((0:info.rows  -1).*info.map_info.dy);

mapInfoOnlyFlag=0;

% get subset bounds if specified
for i=1:numel(varargin)
    if ischar(varargin{i})

        if strcmpi(varargin{i},'pixel_subset')
            subrows = varargin{i+1}(3:4);
            subcols = varargin{i+1}(1:2);

        elseif strcmpi(varargin{i},'map_subset')

            map_subset  = varargin{i+1};
            subcols = (map_subset(1:2)-info.map_info.mapx)./...
                info.map_info.dx+1;
            subrows = (info.map_info.mapy - map_subset([4,3]))./...
                info.map_info.dy+1;
            subcols= round(subcols);
            subrows = round(subrows);

            subcols(subcols < 1) = 1;
            subrows(subrows < 1) = 1;
            subcols(subcols > info.cols) = info.cols;
            subrows(subrows > info.rows) = info.rows;

        elseif strcmpi(varargin{i},'mapinfoonly')
            mapInfoOnlyFlag=1;

        elseif strcmpi(varargin{i},'target_projstr')
            target_projstr = varargin{i+1};
            target_projstr_varargindex = i;
            
            save_cwd = pwd;
            cd(fileparts(mfilename('fullpath')));
            cmd = sprintf('python proj_issame.py "%s" "%s" ', rasterFile, target_projstr);
            [status, cmdout] = system(cmd);
            if ~isempty(cmdout)
                fprintf([cmdout,'\n']);
            end
            if status == 0
                clear target_projstr;
            end
            cd(save_cwd);
        end
    end
end

if exist('target_projstr', 'var')
    fprintf('Reprojecting raster on-the-fly to target projection: %s\n', target_projstr);
    rasterFile_local = [];
    tempdir = getTempDir();
    [~, rasterFname, rasterFext] = fileparts(rasterFile);
    
    rasterFile_local = fullfile(tempdir, [rasterFname, rasterFext]);
    if strcmpi(rasterFile_local, rasterFile)
        rasterFile_local = [];
    else
        copyfile(rasterFile, rasterFile_local);
        rasterFile = rasterFile_local;
    end
    
    [~,~, unitsPerMeter] = getProjInfo(target_projstr, rasterFile);
    
    rasterFile_proj = fullfile(tempdir, sprintf('%s_reproj%s', rasterFname, rasterFext));
    
    if endsWith(rasterFile, '_dem.tif') || endsWith(rasterFile, '_dem_smooth.tif')
        interp_str = 'bilinear';
    elseif endsWith(rasterFile, '_matchtag.tif') || endsWith(rasterFile, '_matchtag_mt.tif')
        interp_str = 'near';
    elseif endsWith(rasterFile, '_ortho.tif')
        interp_str = 'cubic';
    elseif endsWith(rasterFile, 'mask.tif')
        interp_str = 'near';
    else
        error(['Unable to determine interpolation method ', ...
               'for reprojection of %s'], rasterFile);
    end
    
    cmd = [fullfile(gdalpath, 'gdalwarp'), ' -q -co tiled=yes -co compress=lzw -co bigtiff=if_safer '];
    cmd = [cmd, sprintf('-t_srs "%s" -r %s ', target_projstr, interp_str)];
    cmd = [cmd, sprintf('-tap -tr %d %d ', info.map_info.dx*unitsPerMeter, info.map_info.dy*unitsPerMeter)];
%     if exist('map_subset', 'var')
%         % To be written.
%         ;
%     end
    cmd = [cmd, sprintf('"%s" "%s" ', rasterFile, rasterFile_proj)];
    
    [status, cmdout] = system(cmd);
    if ~isempty(cmdout)
        fprintf([cmdout,'\n']);
    end
    
    varargin(target_projstr_varargindex:target_projstr_varargindex+1) = [];
    I = readGeotiff(rasterFile_proj, varargin);
    
    % Create ENVI header for the target projection and keep it for later use when writing output rasters.
    tempdir = getTempDir();
    target_proj_envi = fullfile(tempdir, [strrep(target_projstr, ':', '_'), '.envi']);
    if exist(target_proj_envi, 'file') ~= 2
        fprintf('Creating dummy raster for reference envi header info in later writeGeotiff operations\n');
        cmd = sprintf('gdalwarp -overwrite -q -of ENVI -ts 1 1 "%s" "%s" ', rasterFile_proj, target_proj_envi);
        [status, cmdout] = system(cmd);
        if ~isempty(cmdout)
            fprintf([cmdout,'\n']);
        end
        if ~exist(target_proj_envi, 'file')
            error('gdalwarp call failed: %s', cmd);
        end
    end
    
    delete(rasterFile_proj);
    
    if ~isempty(rasterFile_local)
        delete(rasterFile_local);
    end
    
    return;
end

% define map vectors
I.x = x(subcols(1):subcols(2));
I.y = y(subrows(1):subrows(2));
I.z = [];

if ~mapInfoOnlyFlag
    % read image
    I.z = imread(rasterFile,'PixelRegion',{subrows,subcols});
end


I.info = info;
I.Tinfo = Tinfo;





