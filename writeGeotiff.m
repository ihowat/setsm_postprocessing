function writeGeotiff(OutFileName,x,y,z,fmt,nodata,projstr)
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
    gdalpath = '/Library/Frameworks/GDAL.framework/Versions/Current/Programs/';
end

%if ismac
%	tempfile =  [tempname(tempdir),'.envi'];
%else
%	 tempfile =  [tempname('/scratch/tmp'),'.envi'];
%end

outdir=OutFileName(1:find(OutFileName=='/',1,'last'));
tempfile =  [tempname(outdir),'.envi'];

if strcmpi('polar stereo south',projstr) || strcmpi('polar stereo north',projstr)
    enviwrite(tempfile,x,y(:),z,'format',fmt,'proj',projstr);
else
    if ~isempty(findstr(projstr,'North'))
        zone = str2num(strrep(projstr,'North',''));
        hemi = 'north';
    else
        zone = str2num(strrep(projstr,'South',''));
        hemi = 'south';
    end
    enviwrite(tempfile,x,y(:),z,'format',fmt,'proj','UTM','hemi',hemi,'zone',zone);
end
%enviwrite(tempfile,x,y(:),z,'format',fmt,'proj','polar stereo south');
%enviwrite(tempfile,x,y(:),z,'format',fmt,'proj','UTM','hemi','south','zone',23);


% confirm file write
if ~(exist(tempfile,'file') && exist([tempfile,'.hdr'],'file'))
    disp('Warning, envi file not written in getImage, stopping')
    keyboard
end

system([gdalpath ,'gdal_translate -co bigtiff=if_safer -co compress=lzw -co tiled=yes -a_nodata ',...
    num2str(nodata),' ',tempfile,' ', OutFileName]);

delete([tempfile,'*'])
