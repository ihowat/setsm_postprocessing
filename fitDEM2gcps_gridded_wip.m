function [dzfit,sf]=fitDEM2gcps(x,y,z,px,py,pz,varargin)
% fitDEM2gcps fit surface to offsets between DEM and ground control points
%
% dzfit=fitDEM2gcps(x,y,z,px,py,pz) returns an array dzfit that is a
% quadratic surface fit to the offsets between DEM x,y,z and points
% px,py,pz. The offset can be applied as z-dzfit.
%
% dzfit=fitDEM2gcps(...,'resizeFactor',val) returns an offset array that is
% size * resizeFactor, saving computation time/memory. The result can be
% resized to the orginal DEM size with imresize(dzfit,size(z)).

minPoints = 1000;
%minPoints = 10000;
%minLand = -1;
minLand = 0.95;

%resizeFactor=0.01;
% standard resizeFactor is for 10m, scale for 2m and other res
resizeFactor_10m = 0.01;
res = x(1,2) - x(1,1);
resizeFactor = min(1.0, resizeFactor_10m * (res/10));

% get varargins
n = find(strcmpi(varargin,'resizefactor'));
if ~isempty(n)
    resizeFactor=varargin{n+1};
end

n = find(strcmpi(varargin,'landMask'));
if ~isempty(n)
    land=varargin{n+1};
end

% interpolate dem to control point coordinates
zi = interp2(x,y,z,px,py,'*linear');

if exist('land','var')
    if nnz(land)/numel(land) < minLand
        warning('not enough land in tile to create accurate surface, returning')
        dzfit = [];
        sf=[];
        return
    end

    if any(~land(:))
        land = interp2(x,y,land,px,py,'*nearest',0);
        zi(land ~= 1) = NaN;
        clear land
    end
end

% subtract to get point offsets
dz = zi - pz;

% remove outliers and nan points
[~,n] = rmoutliers(dz);
%[~,n] = rmoutliers(dz,'quartiles');
n = n | isnan(dz);

if sum(~n) < minPoints
    warning('too few control points (%d), returning',sum(~n))
    dzfit = [];
    sf=[];
    return
end

px(n) = [];
py(n) = [];
dz(n) = [];

% resize the DEM coordinates to calculate the surface
x = imresize(x,resizeFactor);
y = imresize(y,resizeFactor);
grid_res = x(1,2) - x(1,1);
x_min = min(x);
x_max = max(x);
y_min = min(y);
y_max = max(y);

% build solution grid
[x,y] = meshgrid(x,y(:));

% interpolate point offsets to solution grid
warning off
F = scatteredInterpolant(double(px),double(py),double(dz));
warning on
x_flat = x(:);
y_flat = y(:);
dzi = F(x_flat,y_flat);

for i=1:length(dzi)
    grid_x = x_flat(i);
    grid_y = y_flat(i);
    if any((abs(px - grid_x) < (grid_res/2)) & (abs(py - grid_y) < (grid_res/2)))
        ;
    elseif grid_x == x_min || grid_x == x_max
        if mod((grid_y - y_min) / grid_res, 2) == 0
            dzi(i) = 0;
        end
    elseif grid_y == y_min || grid_y == y_max
        if mod((grid_x - x_min) / grid_res, 2) == 0
            dzi(i) = 0;
        end
    else
        dzi(i) = NaN;
    end
end

px = x_flat;
py = y_flat;
dz = dzi;

n = isnan(dz);
px(n) = [];
py(n) = [];
dz(n) = [];

% fit a quadratic surface to point offsets
warning off % get precision and fit warnings, just ignore them
sf = fit([px, py],dz,'poly22');
warning on

% calculate surface
dzfit = sf(x,y);










