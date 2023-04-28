function [z_filled,water_mask] = fillWater(tileMatFile,maskFile,refDemFile,varargin)
%input args
% tileMatFile
% maskefile
% refDemFile

% output5
% z_filled

% %externtal functs
% readGeotiff
% tile
% retile
% inpaint_nans

z = [];
z_mad = [];
N = [];
x = [];
y = [];
if length(varargin) >= 1
    z = varargin{1};
end
if length(varargin) >= 3
    z_mad = varargin{2};
    N = varargin{3};
end
if length(varargin) >= 5
    x = varargin{4};
    y = varargin{5};
end

if any(cellfun(@isempty,{z, z_mad, N, x, y}))
    % read tile DEM mat file
    m = matfile(tileMatFile);
    if isempty(z)
        z = m.z;
    end
    if isempty(z_mad)
        z_mad = m.z_mad;
    end
    if isempty(N)
        N = m.N;
    end
    if isempty(x)
        x = m.x;
    end
    if isempty(y)
        y = m.y;
    end
end

% load mask
C = readGeotiff(maskFile);

% set no data (which should be ocean) to water
C.z(C.z == 0) = 80;

% make sure mask is same grid as dem (just do nn interp)
C = interp2(C.x,C.y(:),C.z,x,y(:),'*nearest');

% define water mask using classification and MAD and repeat criteria
M = ( C == 80 | C == 0 ) & ...
    ( (z_mad > 0.3 & N > 4 ) | N <= 4 );

% remove small non water clusters
M = ~bwareaopen(~M,5);

% read z from matfile and apply water and no data (ocean) mask to z
z_masked = z;
z_masked(M) = NaN;
water_mask = M;
clear z;

% read reference DEM
R = readGeotiff(refDemFile);
R.z(R.z < -200) = NaN;

% make sure reference DEM is same grid as dem
R = interp2(R.x,R.y(:),R.z,x,y(:),'*bilinear');

% Fill the DEM with reference DEM using the Delta method

% difference map between z_masked and reference DEM with NaNs in voids
dz = z_masked - R;

%interpolate differences across void, breaking dz into tiles to speed up
A = tile(dz,5,5,200);

% initiate filled cell array of tiles
B = cell(size(A));

% loop through cells to interpolate voids
i=1;
for i=1:numel(A)
    fprintf('%d of %d\n',i, numel(A));
    if ~any(~isnan(A{i}(:))) % no non-nans (maybe ocean?), fill with zeros
        B{i}= zeros(size(A{i}),'single');
    elseif ~any(isnan(A{i}(:))) % no void, skip
         B{i}=A{i};
    else  % mix of data and nans, interpolate
        B{i}=single(inpaint_nans(double(A{i}),2));
    end
end
clear A

% mosiac filled dz tiles 
 dz_filled = retile(B,200,'linear');

% calculate filled DEM
z_filled = R + dz_filled;

