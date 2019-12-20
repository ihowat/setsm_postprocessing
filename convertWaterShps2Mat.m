
waterTileDir = 'E:\pgc_projects\arcticdem\coastline\water_tiles';
% waterTiles = dir([waterTileDir,'\*coast.shp']);
% waterTiles = cellfun(@(x) [waterTileDir,'\',x], {waterTiles.name}, 'uniformOutput',false);
% 
% for i=1:numel(waterTiles)
%     fprintf('Converting tile %s\n',waterTiles{i});
%     shp = shaperead(waterTiles{i});
%     A= cellfun(@(x,y) polyshape(x,y),{shp.X},{shp.Y},'uniformoutput',0);
%     coastlinePoly = A{1};
%     j=2;
%     for j=1:length(A)
%         coastlinePoly(j) = A{j};
%     end
%     clear A shp
%     mFile = strrep(waterTiles{i},'.shp','.mat'); 
%     save(mFile,'coastlinePoly');
%     clear coastlinePoly
% end

clear waterTiles
waterTiles = dir([waterTileDir,'\*lakes_gte2000.shp']);
waterTiles = cellfun(@(x) [waterTileDir,'\',x], {waterTiles.name}, 'uniformOutput',false);

for i=1:numel(waterTiles)
    fprintf('Converting tile %s\n',waterTiles{i});
    shp = shaperead(waterTiles{i});
    A= cellfun(@(x,y) polyshape(x,y),{shp.X},{shp.Y},'uniformoutput',0);
    lakePoly = A{1};
    j=2;
    for j=1:length(A)
        lakePoly(j) = A{j};
    end
    clear A shp
    mFile = strrep(waterTiles{i},'_gte2000.shp','.mat'); 
    save(mFile,'lakePoly');
    clear lakePoly
end
