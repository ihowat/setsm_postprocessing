function addLandMask2REMATile(fileName,tileDefs)

% if ismac
%     tileDefFile = '/Users/ihowat/data4/REMA/rema_tile_definitions.mat';
% else
%     tileDefFile = '/fs/byo/howat-data4/REMA/rema_tile_definitions.mat';
% end

% load tile definition file into structure
if ~isstruct(tileDefs)
    if exist(tileDefs,'file')
        tileDefs=load(tileDefFile);
        if isfield(tileDefs,'I') && ~isfield(tileDefs,'tileName')
            tileDefs.tileName=tileDefs.I;
        end
    else
        error('tileDefs must be structure or valid file name')
    end
end

m=matfile(fileName);

if any(strcmp(fields(m),'land'))
    fprintf('land variable already exists, skipping\n')
    return
end

[~,tileName] = fileparts(fileName);
if strcmp(tileName(6), 's')
    tileName = tileName(1:6);
else
    tileName = tileName(1:5);
end

% find index of tile in tile def database
tileInd = find(strcmp(tileDefs.tileName,tileName));

if isempty(tileInd)
    error('tile name %s not found in %s',tileName,tileDefFile)
end

x=m.x;
y=m.y;

res = x(2)-x(1);
x0=min(x);
y0=min(y);
x1=max(x);
y1=max(y);


m.Properties.Writable = true;

m.land = polyshape2maskmap(...
    tileDefs.coastlinePolyshape(tileInd),[x0 x1 y0 y1],res);


