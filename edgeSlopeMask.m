function [M,dx,dy] = edgeSlopeMask(varargin)
% ESDGESLOPEMASK mask artificats on edges of DEM using high slope values
%
% M = edgeSlopeMask(x,y,z) masks bad edge values in the DEM with coordinate
% vectors x any y.
%
% [M,dx,dy] = edgeSlopeMask(x,y,z) returns the gradient arrays
%
% Ian Howat, ihowat@gmail.com
% 24-Jul-2017 15:39:07
% 04-Oct-2017 15:19:47: Added option/value argument support

%default parameters (designed for data resampled to 8m from 2m)
avg_kernel_size = 21;%size of kernel for calculating mean slope, use 5 for 8m source data
dilate_bad = 13;%dialates masked pixels by this number of surrounding pixels
min_data_cluster = 1000;%isolated clusters of data smaller than this number will be removed

% parse inputs
x = varargin{1};
y = varargin{2};
z = varargin{3};

n = find(strcmpi(varargin,'avg_kernel_size'));
if any(n); avg_kernel_size = varargin{n+1}; end

n = find(strcmpi(varargin,'dilate_bad'));
if any(n); dilate_bad = varargin{n+1}; end

n = find(strcmpi(varargin,'min_data_cluster'));
if any(n); min_data_cluster = varargin{n+1}; end

% slope of z
[dx,dy] = gradient(z,x,y);

% aspect and grade (aspect is used in the cloud filter)
[~,r] = cart2pol(dx,dy);

% mean slope
avg_kernel =ones(avg_kernel_size)/(avg_kernel_size.^2);
Mr = conv2(r,avg_kernel,'same');

% mask mean slopes greater than 1
M = Mr < 1;

% dilate high mean slope pixels and set to false
M(imdilate(Mr > 1,ones(dilate_bad))) = false;

clear Mr

% no data check
if ~any(M(:)); fprintf('boundary filter removed all data\n'); return; end;

% fill interior holes since we're just looking for edges here.
M = imfill(M,'holes');

% get rid of isolated little clusters of data
M = bwareaopen(M, min_data_cluster);

% no data check
if ~any(M(:)); fprintf('boundary filter removed all data\n'); return; end;

% find data coverage boundaries
B = bwboundaries(M, 8, 'noholes');

% merge all the data clusters
B = cell2mat(B);

% find outer data boundary
k = boundary(B(:,2),B(:,1),0.5); %allows for some inward "bending" to conform to the data.
%k = convhull(B(:,2),B(:,1));% straight line trace of the outer edge.

% build edge mask
M = poly2mask(B(k,2),B(k,1), size(M,1),size(M,2));
