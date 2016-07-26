function z0 = edgeWarp(z0,z1,varargin)
% EDGEWARP Merge two DEMS, forcing vertical edge alignment
%
%W = edgeWarp(z0,z1,buff) where z0 and z1 are equally-sized arrays
%to be merged with z0 taking presidence. The arrays are merged by applying
%a linearly-increasing offset from 0, buff pixels away from the edge of
%ovlap, to the full offset at the edge of z1. The larger buff, the wider
%the region of warping and the smoother the transition will be for larger
%offsets. Treats NaNs as nodata.

% set default buffer size or set to varargin
buff=51;
if nargin==3
    buff=varargin{1};
end

% make BW array of region of overlap
A= ~isnan(z0) & ~isnan(z1);

% if no overlap, insert z1 data and return
if ~any(A(:))
        z0(~isnan(z1))=z1(~isnan(z1));
        return
end

% define BW array of warp region
warpRegion=imdilate(A,ones(buff)) & ~A & ~isnan(z1);

% define BW of boundary pixels around warp region for interpolation
B = imdilate(warpRegion,ones(11)) & ~warpRegion;

% set pixel indices
nWarp=find(warpRegion); % warp-region pixels
ndz=find(A & B); % difference region boundary pixels
nz1=find(~isnan(z1) & ~A & B); % z1 (zero) difference boundary pixels

% row/column subscripts of boundaries
[R,C]=ind2sub(size(z0),[ndz;nz1]); 

% difference overlap pixels
dz = z0(ndz)-z1(ndz);

% build interpolant using boundary pixels
warning off
F = scatteredInterpolant(C,R,[double(dz);zeros(size(nz1))]);
warning on

% row/column subscripts of warp region
[R,C]=ind2sub(size(z0),nWarp);

% interpolate warp region, add difference to z1 and insert into z0;
z0(nWarp)=z1(nWarp)+F(C,R);

% insert remaining z1 into empty z0
z0(isnan(z0))=z1(isnan(z0));

