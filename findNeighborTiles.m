function [ntop,nright,ntopRight,nbottomRight] = findNeighborTiles(fileNames)
% findNeighborTiles find neighboring tiles based on map coordinates
%
% [ntop,nright,ntopRight,nbottomRight] = findNeighborTiles(fileNames) uses
% the x y coordinate vectors in each tile matfile to find the neighbors to
% the top, right, top-right and bottom-right.Tile boundaries must overlap.
% Returns n by 2 arrays of indices of the neighbors in the list of file
% names, so that fileNames(ntop) returns a list of top/bottom nighboring 
% files.

fileNames=fileNames(:);

% loop through each tile and get extents
x0 = nan(size(fileNames));
x1= nan(size(fileNames));
y0= nan(size(fileNames));
y1= nan(size(fileNames));
i=1;
for i=1:length(fileNames)
    m = matfile(fileNames{i});
    x0(i)  = min(m.x);
    x1(i) = max(m.x);
    y0(i)  = min(m.y);
    y1(i) = max(m.y);
end
    
% using ranges, find neighbors on each side
ntop         = nan(length(fileNames),2);
nright       = ntop;
ntopRight    = ntop;
nbottomRight = ntop;

i=1;
for i=1:length(fileNames)
    
    mnx=(x0(i)+x1(i))./2;
    mny=(y0(i)+y1(i))./2;
    
    %top
    n = find(y0(i) < y0 & y1(i) >= y0 & mnx > x0 & mnx < x1);
    if ~isempty(n)
        ntop(i,:) = [i,n];
    end

    %right
    n = find(x0(i) < x0 & x1(i) >= x0 & mny > y0 & mny < y1);
    if ~isempty(n)
        nright(i,:) = [i,n];
    end
    
    %top-right
    n = find(y0(i) < y0 & y1(i) >= y0 & x0(i) < x0 & x1(i) >= x0);
    if ~isempty(n)
        ntopRight(i,:) = [i,n];
    end

    %bottom-right
    n = find(y0(i) > y0 & y0(i) <= y1 & x0(i) < x0 & x1(i) >= x0);
    if ~isempty(n)
        nbottomRight(i,:) = [i,n];
    end  
end

ntop(isnan(ntop(:,1)),:) = [];
nright(isnan(nright(:,1)),:) = [];
ntopRight(isnan(ntopRight(:,1)),:) = [];
nbottomRight(isnan(nbottomRight(:,1)),:) = [];


    