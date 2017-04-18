function writeTileToTif(tilef,projstr)
% Write dem or reg_dem matfiles to tif

%% Crop Buffers and Write Tiles To Geotiffs 
fprintf('writing tif %d of %d\n',i,length(tilef))

if exist(strrep(tilef,'dem.mat','reg_dem.mat'),'file')
    fi=strrep(tilef,'dem.mat','reg_dem.mat');
else
    fi=tilef;
end

fprintf('source: %s\n',fi);

load(fi,'x','y');
% crop buffer tile
x=x(101:end-100);
y=y(101:end-100);

OutDemName = strrep(fi,'.mat','.tif');
if ~exist(OutDemName,'file');
    
    load(fi,'z');
    z=z(101:end-100,101:end-100);
    z(isnan(z)) = -9999;
    
    writeGeotiff(OutDemName,x,y,z,4,-9999,projstr)
    clear z
end

OutMatchtagName = strrep(fi,'dem.mat','matchtag.tif');
if ~exist(OutMatchtagName,'file');
    load(fi,'mt');
    mt =mt(101:end-100,101:end-100);
    writeGeotiff(OutMatchtagName,x,y,mt,1,0,projstr)
    clear mt
end

% OutOrthoName = strrep(fi,'dem.mat','ortho.tif');
% if ~exist(OutOrthoName ,'file');
%     load(fi,'or');
%     or =or(101:end-100,101:end-100);
%     writeGeotiff(OutOrthoName ,x,y,or,2,0,projstr)
%     clear or
% end;

% OutDaynumName = strrep(fi,'dem.mat','daynum.tif');
% if ~exist(OutDaynumName,'file');
%     load(fi,'dy');
%     dy =dy(101:end-100,101:end-100);
%     writeGeotiff(OutDaynumName,x,y,dy,2,0,projstr)
%     clear dy
% end

clear x y