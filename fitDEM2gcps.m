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

% get varargins
resizeFactor=0.01;
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

% fit a quadratic surface to point offsets
warning off % get precision and fit warnings, just ignore them
sf = fit([px, py],dz,'poly22');
warning on

% resize the DEM coordinates to calculate the surface
x = imresize(x,resizeFactor);
y = imresize(y,resizeFactor);

% build solution grid
[x,y] = meshgrid(x,y(:));

% calculate surface
dzfit = sf(x,y);










