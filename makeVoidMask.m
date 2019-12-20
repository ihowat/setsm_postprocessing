function makeVoidMask(tileName)
% makeVoidMask
% load a dem file and make a void mask by drawing polygons on the hillshade
% image. Outputs the voidMask and polygons to a .mat file. The output
% will then be used to fill gaps in the dem using "fillGaps.m". Currently
% set up to read a 2m tile and rescale to 10m for editing.

demDir = 'V:\pgc\data\scratch\claire\pgc\arcticdem\mosaic\2m_v4';
demFiles = dir([demDir,'/',tileName,'/',tileName,'*_2m_dem.tif']);
demFiles = cellfun(@(x) [demDir,'/',tileName,'/',x], {demFiles.name}, 'uniformOutput', false);
%demFile='V:\pgc\data\scratch\claire\pgc\arcticdem\mosaic\2m_v4\44_18\44_18_1_1_2m_dem.tif';

fprintf('Loading tile %s water mask\n',tileName)
waterMask='~/49_10_water.tif';
waterTileDir='V:\pgc\data\scratch\claire\pgc\arcticdem\coastline\global_surface_water\tiles_v2';
waterMask = dir([waterTileDir,'/',tileName,'_water.tif']);
waterMask = cellfun(@(x) [waterTileDir,'/',x], {waterMask.name}, 'uniformOutput',false);
if isempty(waterMask)
    fprintf('Tile %s water mask file does not exist in %s\n',tileName,waterTileDir)
    return
end
waterMask = waterMask{1};

for i=1:length(demFiles)
    demFile = demFiles{i};
    voidMaskFile=strrep(demFile,'.tif','_voidMask.mat');

    %% load data and resize to 10m
    fprintf('loading dem file: %s\n',demFile)
    dem = readGeotiff(demFile);
    sz = size(dem.z);
    
    % convert nodata to nans
    dem.z(dem.z == -9999) = NaN;

    %downsize dem to 10m
    dem10.x = imresize(dem.x,.2);
    dem10.y = imresize(dem.y,.2);
    dem10.z = imresize(dem.z,.2);
    
    clear dem

    % build hillshade
    h10 = hillshade(dem10.z,dem10.x,dem10.y);

    %% build land10 mask
    land = readGeotiff(waterMask);
    land.z = land.z==0;

    land10=interp2(land.x,land.y(:),single(land.z),...
        dem10.x,dem10.y(:),'*nearest');

    land10(isnan(land10)) = 0;
    land10 = logical(land10);

    clear land

%     % make polyshape of this tile with buffer to ensure coverage of border
%     % cells
%     res=100;
%     x0  = min(dem10.x);
%     x1  = max(dem10.x);
%     y0  = min(dem10.y);
%     y1  = max(dem10.y);
%     tilePoly = polyshape([x0-res;x0-res;x1+res;x1+res],[y0-res;y1+res;y1+res;y0-res]);
% 
%     land10 =  false(length(dem10.y),length(dem10.x)); % land10 mask
% 
%     for i=1:length(coastlinePoly)
%         if overlaps(tilePoly,coastlinePoly(i))
% 
%             landPoly = intersect(tilePoly,coastlinePoly(i));
% 
%             NR = [0;find(isnan(landPoly.Vertices(:,1)));...
%                 length(landPoly.Vertices(:,1))+1];
% 
%             for nr = 1:length(NR)-1
%                 xp = landPoly.Vertices(NR(nr)+1:NR(nr+1)-1,1);
%                 yp = landPoly.Vertices(NR(nr)+1:NR(nr+1)-1,2);
%                 land10(roipoly(dem10.x,dem10.y,land10,xp,yp))=true;
%             end
% 
%         end
%     end

    %% build void mask with manual editing

    % make mask of missing data that's not water
    voidMask = ~isnan(dem10.z) | ~land10;

    [voidMask,voidPolys] = manualEdit(dem10.x,dem10.y,h10,voidMask,land10);
    close

    % resize to 2m
    voidMask = imresize(voidMask,sz,'nearest');

    save(voidMaskFile,'voidMask','voidPolys');
    
end
