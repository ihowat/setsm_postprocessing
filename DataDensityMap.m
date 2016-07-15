function P = DataDensityMap(varargin)
% DataDensityMap density of data points within an array
% 
% P = DataDensityMap(m) returns a size(m) array of the fraction of nodes
% containing ones within an 11 by 11 kernel, where m is a mask array of 0
% for no data and 1 for data.
%
% P = DataDensityMap(m,n) uses an n-sized kernal (default = 11)

n = 11;

m = varargin{1};
if nargin == 2;
    n = varargin{2};
end

P =conv2(double(m), ones(n),'same');
P = P./n.^2;