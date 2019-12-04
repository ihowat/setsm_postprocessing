function N = stripSearch(varargin)
% stripSearch filters strip metastruct to strips overlapping tile bounds
%
% n = stripSearch(x,y,tilex0,tilex1,tiley0,tiley1) return the index of
% where x and y are cell arrays containing the x and y footprint vertices.
%
% Ian Howat, ihowat@gmail.com, Ohio State
%03-Dec-2019 12:31:01: 1. takes a polyshape as search region
%                      2. uses polyshape function "overlaps", no longer uses intersection

x=varargin{1};
y=varargin{2};

if nargin == 3 % polyShape provided
    
    tilePoly = varargin{3};
    
    tilex0 = nanmin(tilePoly.Vertices(:,1));
    tilex1 = nanmax(tilePoly.Vertices(:,1));
    tiley0 = nanmin(tilePoly.Vertices(:,2));
    tiley1 = nanmax(tilePoly.Vertices(:,2));
    
elseif  nargin == 4  % x,y vertices provided
    %     tilevx = varargin{3};
    %     tilevy = varargin{4};
    
    tilePoly = polyShape(varargin{3},varargin{4});
    
    tilex0 = min(tilevx);
    tilex1 = max(tilevx);
    tiley0 = min(tilevy);
    tiley1 = max(tilevy);
    
elseif nargin == 6  % x,y range provided
    
    tilex0 = varargin{3};
    tilex1 = varargin{4};
    tiley0 = varargin{5};
    tiley1 = varargin{6};
    
    % make tile boundary polygon
    %     tilevx = [tilex0;tilex0;tilex1;tilex1;tilex0];
    %     tilevy = [tiley0;tiley1;tiley1;tiley0;tiley0];
    
    tilePoly = polyShape([tilex0;tilex0;tilex1;tilex1;tilex0],...
        [tiley0;tiley1;tiley1;tiley0;tiley0]);
    
else
    
    error('nargin must be 1,2 or 4\n')
end


% quick search: find strips within range of this tile. This does not
% account for background area of around strips but just pairs them down to
% speed the poly intersection loop

xmax=cellfun(@max,x);
xmin=cellfun(@min,x);
ymax=cellfun(@max,y);
ymin=cellfun(@min,y);


n = xmax > tilex0 & xmin < tilex1 & ...
    ymax > tiley0 & ymin < tiley1;


N=find(n);

% if no overlap, return
if isempty(N); return; end

x=x(n);
y=y(n);

% search for all strip footprints overlapping this tile
n=false(size(x));

for i=1:length(n)
    
    stripPoly = polyshape(x{i},y{i});
    
    n(i) = overlaps(stripPoly,tilePoly);
    
    
    %     n(i) = any(inpolygon(x{i},y{i},tilevx,tilevy)) | ...
    %         any(inpolygon(tilevx,tilevy,x{i},y{i}));
    %
    %     if ~n(i)
    %
    %         [xtest,~] = intersections(x{i},y{i},tilevx,tilevy);
    %
    %         if ~isempty(xtest)
    %             n(i) = true;
    %         end
    %
    %     end
    %
    
end

N=N(n);