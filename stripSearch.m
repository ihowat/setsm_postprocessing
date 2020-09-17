function N = stripSearch(varargin)
% stripSearch filters strip metastruct to strips overlapping roi bounds
%
%n = stripSearch(x,y,roiPolyShape) return the index of polygons that
%overlap roi poly shape, where x and y are cell arrays containing the x and
%y footprint vertices.
%
%n = stripSearch(x,y,roix,roiy) where roix and roiy are cell arrays 
% containing the x and y footprint vertices of the search roi.
%
%n = stripSearch(x,y,roix0,roix1,roiy0,roiy1) where the rois are defined
%by the min/max ranges in x and y.
%
% stripSearch(...,'minOverlapFraction',F) specifies the minimum fraction of
% overlap to accept.
%
%Ian Howat, ihowat@gmail.com, Ohio State
%03-Dec-2019 12:31:01: 1. takes a polyshape as search region
%                      2. uses polyshape function "overlaps", no longer 
%                         uses intersection.
%10-Aug-2020 15:50:31: 1. Added minOverlapFraction option.
%                      2. Restructured varargin parsing


x=varargin{1};
y=varargin{2};

if isa(varargin{3},'polyshape') % polyShape provided
    
    roiPoly = varargin{3};
    
elseif length(varargin{3}) > 1 % x,y vertices provided
    
    roiPoly = polyshape(varargin{3},varargin{4});
    
elseif length(varargin{3}) == 1    % x,y range provided
    
    roiPoly = polyshape([varargin{3};varargin{3}; varargin{4};...
        varargin{4};varargin{3}],[varargin{5};varargin{6};varargin{6};...
        varargin{5};varargin{5}]);    
else
    
    error('input must be a polyshape, vertices or x,y ranges \n')
end

%check minimum overlap fraction argument
minOverlapFraction=0;
n = find(strcmpi(varargin,'minOverlapFraction'));
if n
    minOverlapFraction=varargin{n+1};
end


% quick search: find strips within range of this roi. This does not
% account for background area of around strips but just pairs them down to
% speed the poly intersection loop

xmax=cellfun(@max,x);
xmin=cellfun(@min,x);
ymax=cellfun(@max,y);
ymin=cellfun(@min,y);

roix0 = nanmin(roiPoly.Vertices(:,1));
roix1 = nanmax(roiPoly.Vertices(:,1));
roiy0 = nanmin(roiPoly.Vertices(:,2));
roiy1 = nanmax(roiPoly.Vertices(:,2));

n = xmax > roix0 & xmin < roix1 & ...
    ymax > roiy0 & ymin < roiy1;

N=find(n);

% if no overlap, return
if isempty(N); return; end

x=x(n);
y=y(n);

% search for all strip footprints overlapping this roi
n=false(size(x));

for i=1:length(n)
    
    stripPoly = polyshape(x{i},y{i});
    
    n(i) = overlaps(stripPoly,roiPoly);
    
end

N=N(n);

% if a minimum % overlap is supplied, check for it
if minOverlapFraction > 0
    x=x(n);
    y=y(n);
    
    
    % search for all strip footprints overlapping this roi
    n=false(size(x));
    roiArea = area(roiPoly);
    for i=1:length(n)
        
        stripPoly = polyshape(x{i},y{i});
        
        a = intersect(stripPoly,roiPoly);
        
        n(i) = area(a)/roiArea > minOverlapFraction;
        
        
    end
    
    N=N(n);
    
end