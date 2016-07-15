function writeGeotiff(OutFileName,x,y,z,fmt,nodata,projstr)

gdalpath =[];
if ismac
    gdalpath = '/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/';
end

%if ismac
%	tempfile =  [tempname(tempdir),'.envi'];
%else
%	 tempfile =  [tempname('/scratch/tmp'),'.envi'];
%end

outdir=OutFileName(1:find(OutFileName=='/',1,'last'));
tempfile =  [tempname(outdir),'.envi'];



enviwrite(tempfile,x,y(:),z,'format',fmt,'proj',projstr);
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
