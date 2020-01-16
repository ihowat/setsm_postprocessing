function batchMergeTileBuffer(f)
% batchMergeTileBuffer: mergeTileBuffer to all neighboring files
%
% batchMergeTileBuffer(f) where f is a cellstr of files to be merged.


f=f(:);

% loop through each tile and get extents
x0 = nan(size(f));
x1= nan(size(f));
y0= nan(size(f));
y1= nan(size(f));
i=1;
for i=1:length(f)
    m = matfile(f{i});
    x0(i)  = min(m.x);
    x1(i) = max(m.x);
    y0(i)  = min(m.y);
    y1(i) = max(m.y);
end
    
% using ranges, find neighbors on each side
nright= nan(length(f),2);
ntop= nan(length(f),2);
i=1;
for i=1:length(f)

    %right
    n = find(x0(i) < x0 & x1(i) > x0 & (y0(i)+y1(i))./2 > y0 & (y0(i)+y1(i))./2 < y1);
    if ~isempty(n)
        nright(i,:) = [i,n];
    end
    
    %top
    n = find(y0(i) < y0 & y1(i) > y0 & (x0(i)+x1(i))./2 > x0 & (x0(i)+x1(i))./2 < x1);
    if ~isempty(n)
        ntop(i,:) = [i,n];
    end
    
end

n0 = [nright(~isnan(nright(:,1)),1);ntop(~isnan(ntop(:,1)),1)];
n1 = [nright(~isnan(nright(:,1)),2);ntop(~isnan(ntop(:,1)),2)];
    

% run mergeTileBuffer on each pair in list
for i=1:length(n0)
    mergeTileBuffer(f{n0(i)},f{n1(i)}); 
end
    
    
    




