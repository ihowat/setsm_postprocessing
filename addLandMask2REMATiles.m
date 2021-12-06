function addLandMask2REMATiles(fileNames)

% fileNames=dir('*10m_reg.mat');
% fileNames={fileNames.name};

if ismac
    tileDefFile = '/Users/ihowat/data4/REMA/rema_tile_definitions.mat';
else
    tileDefFile = '/fs/byo/howat-data4/REMA/rema_tile_definitions.mat';
end

% load tile definition file into structure
tileDefs=load(tileDefFile);

if isfield(tileDefs,'I') && ~isfield(tileDefs,'tileName')
    tileDefs.tileName=tileDefs.I;
end

i=1;
for i=1:length(fileNames)
    
    fprintf('adding land mask to %s...',fileNames{i})
    
    tileName = fileNames{i}(1:5);
    
    % find index of tile in tile def database
    tileInd = find(strcmp(tileDefs.tileName,tileName));
    
    if isempty(tileInd)
        error('tile name %s not found in %s',tileName,tileDefFile)
    end
    
    m=matfile(fileNames{i});
    
    x=m.x;
    y=m.y;
    
    res = x(2)-x(1);
    x0=min(x);
    y0=min(y);
    x1=max(x);
    y1=max(y);
    
    
    m.Properties.Writable = true;
    
    % make land/water(true/false) mask
    %if isfield(tileDefs,'coastlinePolyshape')
    m.land = polyshape2maskmap(...
        tileDefs.coastlinePolyshape(tileInd),[x0 x1 y0 y1],res);
    %else
    %    land= waterTileDir;
    %end
     fprintf('done\n')
end



