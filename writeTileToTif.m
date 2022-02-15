function writeTileToTif(tilef,res,projstr)
% Write dem or reg_dem matfiles to tif

%% Crop Buffers and Write Tiles To Geotiffs
fprintf('writing tif %d of %d\n',i,length(tilef))

if exist(strrep(tilef,'dem.mat','reg_dem.mat'),'file')
    fi=strrep(tilef,'dem.mat','reg_dem.mat');
else
    fi=tilef;
end

fprintf('source: %s\n',fi);

% calc buffer to remove
buffer = floor(200 / res);

load(fi,'x','y');
% crop buffer tile
x=x(buffer+1:end-buffer);
y=y(buffer+1:end-buffer);

OutDemName = strrep(fi,'.mat','.tif');
if ~exist(OutDemName,'file')
    load(fi,'z');
    z=z(buffer+1:end-buffer,buffer+1:end-buffer);
    z(isnan(z)) = -9999;
    writeGeotiff(OutDemName,x,y,z,4,-9999,projstr)
    clear z
end

OutMatchtagName = strrep(fi,'dem.mat','matchtag.tif');
if ~exist(OutMatchtagName,'file')
    load(fi,'mt');
    mt =mt(buffer+1:end-buffer,buffer+1:end-buffer);
    writeGeotiff(OutMatchtagName,x,y,mt,1,0,projstr)
    clear mt
end

% OutOrthoName = strrep(fi,'dem.mat','ortho.tif');
% if ~exist(OutOrthoName ,'file');
%     load(fi,'or');
%     or =or(buffer+1:end-buffer,buffer+1:end-buffer);
%     writeGeotiff(OutOrthoName ,x,y,or,2,0,projstr)
%     clear or
% end;

% OutDaynumName = strrep(fi,'dem.mat','daynum.tif');
% if ~exist(OutDaynumName,'file');
%     load(fi,'dy');
%     dy =dy(buffer+1:end-buffer,buffer+1:end-buffer);
%     writeGeotiff(OutDaynumName,x,y,dy,2,0,projstr)
%     clear dy
% end

hillshade=strrep(OutDemName,'dem.tif','dem_shade.tif');
if ~exist(hillshade,'file');
    system(['gdaldem hillshade -compute_edges -b 1 -of GTiff -co tiled=yes -co compress=lzw -co bigtiff=yes ',...
        OutDemName,' ',hillshade]);
end

clear x y
