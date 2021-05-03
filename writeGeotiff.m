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

if contains(projstr,'UTM','IgnoreCase',true) ...
    || ~isempty(regexp(projstr, '^\s*(?:UTM)?\s*\d{1,2}\s*(?:North|South|N|S)\s*$', 'ignorecase')) ...
    || ~isempty(regexp(projstr, '^\s*(?:UTM)?\s*(?:North|South|N|S)\s*\d{1,2}\s*$', 'ignorecase'))

    projstr_trim = regexprep(projstr,'utm','','ignorecase');

    if contains(projstr,'North','IgnoreCase',true)
        hemi = 'north';
        projstr_trim = regexprep(projstr_trim,'north','','ignorecase');
    elseif contains(projstr,'South','IgnoreCase',true)
        hemi = 'south';
        projstr_trim = regexprep(projstr_trim,'south','','ignorecase');
    elseif contains(projstr,'N','IgnoreCase',true)
        hemi = 'north';
        projstr_trim = regexprep(projstr_trim,'n','','ignorecase');
    elseif contains(projstr,'S','IgnoreCase',true)
        hemi = 'south';
        projstr_trim = regexprep(projstr_trim,'s','','ignorecase');
    else
        error('Cannot parse hemisphere information from UTM ''projstr'': %s', projstr);
    end

    zone = str2num(projstr_trim);
    if isempty(zone)
        error('Cannot parse zone number from UTM ''projstr'': %s', projstr);
    end

    enviwrite(tempfile,x,y(:),z,'format',fmt,'proj','UTM','hemi',hemi,'zone',zone);
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

system([gdalpath ,'$BWPY_PREFIX gdal_translate -co bigtiff=if_safer -co compress=lzw -co tiled=yes -a_nodata ',...
    num2str(nodata),' ',tempfile,' ', OutFileName]);

delete([tempfile,'*'])
