function writeGeotiff(OutFileName,x,y,z,fmt,nodata,projstr,varargin)
% writeGeotiff writes raster to geotiff using GDAL
%
% writeGeotiff(OutFileName,x,y,z,fmt,nodata,projstr) writes the map data in
% x, y, z using the format, nodata and projection tags. 'projstr' can be 
% 'polar stereo north', 'polar stereo south' or, for utm, a string the zone
% number and hemisphere (e.g. '12 North'). This function actually writes
% the raster first to a flat binary in envi format using  enviwrite and
% then calls gdal to perform the conversion. If not in the system path,
% gdalpath needs to be set below.
%
% See enviwrite help for format codes and other options.

% Ian Howat, ihowat@gmail.com, Ohio State


gdalpath =[]; %set to the path of the gdal binary if not in system path.
if ismac
    %gdalpath = '/Library/Frameworks/GDAL.framework/Versions/Current/Programs/';
    gdalpath = '/opt/local/bin/';
end


% See GDAL docs for info on GTiff vs. COG
% https://gdal.org/drivers/raster/gtiff.html
% https://gdal.org/drivers/raster/cog.html
n = find(strcmpi('out_format',varargin));
if ~isempty(n)
    out_format = varargin{n+1};
else
    out_format = 'GTiff';
end
out_format_choices = {'GTiff', 'COG'};
if ~any(strcmp(out_format, out_format_choices))
    error("'out_format' must be one of the following, but was '%s': {'%s'}", out_format, strjoin(out_format_choices, "', '"))
end

n = find(strcmpi('co_predictor',varargin));
if ~isempty(n)
    co_predictor = varargin{n+1};
    if isnumeric(co_predictor)
        co_predictor = num2str(co_predictor);
    end
else
    if strcmp(out_format, 'GTiff')
        co_predictor = '1';
    elseif strcmp(out_format, 'COG')
        co_predictor = 'YES';
    end
end
if strcmp(out_format, 'GTiff')
    % See https://gdal.org/drivers/raster/gtiff.html#creation-options
    co_predictor_choices = {'1', '2', '3'};
elseif strcmp(out_format, 'COG')
    % See https://gdal.org/drivers/raster/cog.html#general-creation-options
    co_predictor_choices = {'NO', 'STANDARD', 'FLOATING_POINT', 'YES'};
end
if ~any(strcmp(out_format, out_format_choices))
    error("'co_predictor' must be one of the following, but was '%s': {'%s'}", out_format, strjoin(co_predictor_choices, "', '"))
end

if strcmp(out_format, 'COG')
    co_tiled_arg = '';
    co_overviews_args = '-co OVERVIEWS=IGNORE_EXISTING';

    n = find(strcmpi('cog_overview_resampling',varargin));
    if ~isempty(n)
        cog_overview_resampling = varargin{n+1};
        co_overviews_args = sprintf('%s -co RESAMPLING=%s', co_overviews_args, cog_overview_resampling);
    else
        error("'cog_overview_resampling' varargin option must be provided when 'out_format' is 'COG'")
    end
else
    co_tiled_arg = '-co TILED=YES';
    co_overviews_args = '';
end


%if ismac
%	tempfile =  [tempname(tempdir),'.envi'];
%else
%	 tempfile =  [tempname('/scratch/tmp'),'.envi'];
%end

outdir=OutFileName(1:find(OutFileName=='/',1,'last'));
tempfile =  [tempname(outdir),'.envi'];

[epsg,utm_hemi,utm_zone] = getProjstrInfo(projstr);

if ~isempty(utm_hemi)
    enviwrite(tempfile,x,y(:),z,'format',fmt,'proj','UTM','hemi',utm_hemi,'zone',utm_zone);
else
    enviwrite(tempfile,x,y(:),z,'format',fmt,'proj',projstr);
end
%enviwrite(tempfile,x,y(:),z,'format',fmt,'proj','polar stereo south');
%enviwrite(tempfile,x,y(:),z,'format',fmt,'proj','UTM','hemi','south','zone',23);


% confirm file write
if ~(exist(tempfile,'file') && exist([tempfile,'.hdr'],'file'))
    disp('Warning, envi file not written in getImage, stopping')
    keyboard
end

gdal_cmd = [...
    gdalpath, 'gdal_translate',...
    ' -of ', out_format,...
    ' -a_srs ', sprintf('EPSG:%d', epsg),...
    ' -a_nodata ', num2str(nodata),...
    ' ', co_tiled_arg, ' -co BIGTIFF=YES',...
    ' -co COMPRESS=LZW -co PREDICTOR=', co_predictor,...
    ' ', co_overviews_args,...
    ' ', tempfile, ' ', OutFileName];
fprintf('%s\n', gdal_cmd);
[status, cmdout] = system(gdal_cmd);
if ~isempty(cmdout)
    disp(cmdout)
end
if status ~= 0
    error('Non-zero exit status (%d) from gdal_translate',status)
end

delete([tempfile,'*'])
