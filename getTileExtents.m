function [x0,x1,y0,y1]=getTileExtents(varargin)
% getTileExtents
%
%[x0,x1,y0,y1]=getTileExtents(tileName,tileDefFile)
%[x0,x1,y0,y1]=getArcticDemTileExtents(...,'buffer',buffer,'quadrant',quadrant_string)
% Input arguments:
%tileName=string tile x,y name (e.g. '47_13')
%tileDefFile = matfile list of tile names and ranges (e.g. '/Users/ihowat/unity-home/earthdem/PGC_Imagery_Mosaic_Tiles_Arctic.mat');
%buffer =  size of tile/subtile boundary buffer in meters (100 m default)

buffer=100; % default size of tile/subtile boundary buffer in meters
quadrant = '';

tileName=varargin{1};
tileDefFile = varargin{2};

n = find(strcmpi('buffer',varargin));
if ~isempty(n)
    buffer = varargin{n+1};
end
n = find(strcmpi('quadrant',varargin));
if ~isempty(n)
    quadrant = varargin{n+1};
end

tileDefs=load(tileDefFile);

% find index of tile in tile def database
if startsWith(tileName,'utm')
    sl = split(tileName,'_');
    tilePrefix = [sl{1},'_'];
    tileName_in_tileDef = strjoin(sl(2:3),'_');
else
    tilePrefix = '';
    tileName_in_tileDef = tileName;
end
tileInd = find(strcmp(tileDefs.I,tileName_in_tileDef));

% get tile boundaries with buffer
x0=tileDefs.x0(tileInd);
y0=tileDefs.y0(tileInd);
x1=tileDefs.x1(tileInd);
y1=tileDefs.y1(tileInd);

switch quadrant
    case lower('1_1') %lower left quadrant
      
        x1 = x0 + (x1-x0)./2 + 1000;
        y1 = y0 + (y1-y0)./2 + 1000;
        
    case lower('2_1') %upper left quadrant
        
        x1 = x0 + (x1-x0)./2 + 1000;
        y0 = y0 + (y1-y0)./2 - 1000;
        
    case lower('1_2') %lower right quadrant
        
        x0 = x0 + (x1-x0)./2 - 1000;
        y1 = y0 + (y1-y0)./2 + 1000;
        
    case lower('2_2') %upper right quadrant
        
        x0 = x0 + (x1-x0)./2 - 1000;
        y0 = y0 + (y1-y0)./2 - 1000;
end

% add buffer
x0=x0-buffer;
y0=y0-buffer;
x1=x1+buffer;
y1=y1+buffer;
