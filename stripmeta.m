function stripmeta(f)
% writes meta data file for strip mosaics of SETSM-generated scenes
%
% stripmeta(f) where f is the file path/*_dem.tif generated a *_meta.txt
% file The file structure must be
% the following:
% 1. strip file located in /path/strips/*dem.tif'
% 2. mosaicking alignment statistics file is /path/region/tif_results/strips/*trans.mat'
% 3. Individual scene file meta data in ../*_meta.txt
%
% Uses gdal to get projection string

% test strip
%f='/data2/ArcticDEM/Greenland/region_03_greenland_southwest/tif_results/strips/WV01_20120619_102001001B432E00_102001001A8B4300_seg1_8m_dem.tif';

gdalpath=[];
if ismac
    gdalpath='/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/';
end


[~,projstr]=system([gdalpath,'gdalsrsinfo -o proj4 ',f]);

projstr=strrep(projstr,'+x_0=1','+x_0=0');
projstr=strrep(projstr,'+y_0=1','+y_0=0');

outfile=strrep(f,'dem.tif','meta.txt');
transfile=strrep(f,'dem.tif','trans.mat');

%outfile='test_meta.txt';

a = readGeotiff(f);
m = a.z ~= -9999;

B = bwboundaries(m, 8, 'noholes'); % find data coverage boundaries

B = cell2mat(B); % merge all the data clusters

if ~isempty(B)

    k = convhull(B(:,2),B(:,1)); % find outer data boundary - the convex hull
                              % gives the straight line trace of the outer 
                              % edge.
                              
    B=DecimatePoly(B(k,:),[1 1]);
                                                    
    x=a.x(B(:,2));
    y=a.y(B(:,1));
else
    x=NaN;
    y=NaN;
end

fid=fopen(outfile,'w');

fprintf(fid,'Strip Metadata \n');
fprintf(fid,'Creation Date: %s\n',datestr(now));
stripdate=dir(transfile); stripdate=stripdate.date;
fprintf(fid,'Strip creation date: %s\n',stripdate);
fprintf(fid,'Strip projection (proj4): %s\n',projstr);
fprintf(fid,'Strip Footprint Vertices \n');
fprintf(fid,'X: '); fprintf(fid,'%d ',x); fprintf(fid,'\n');
fprintf(fid,'Y: '); fprintf(fid,'%d ',y); fprintf(fid,'\n');
fprintf(fid,'\n');

load(transfile);
fprintf(fid,'Mosaicking Alignment Statistics (meters) \n');
fprintf(fid,'scene, rmse, dz, dx, dy\n');
i=1;
N=length(scene);
for i=1:N
    fprintf(fid,'%s %.2f %.4f %.4f %.4f\n',scene{i},rmse(i),trans(:,i)'); 
end

fprintf(fid,'\n');
fprintf(fid,'Scene Metadata \n');
fprintf(fid,'\n');
fclose(fid);

i=1;
N=length(scene);
for i=1:N

	sceneMetaFile=[strrep(fileparts(f),'strips','tif_results'),'/',strrep(scene{i},'dem.tif','meta.txt')];

    system(['echo scene ',num2str(i),' name=',scene{i},' >> ',outfile]);
    if exist(sceneMetaFile,'file')
        system(['cat ',sceneMetaFile,' >> ',outfile]);
    else
         system(['echo ',sceneMetaFile,' not found >> ',outfile]);
    end
    system(['echo " " >> ',outfile]);
end   
    
    
    
    
