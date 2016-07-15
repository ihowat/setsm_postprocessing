function M = edgeMask(varargin)
% edgeMask returns mask for bad edges using the matchag field
%
% m1 = edgeMask(m0) where m0 is the matchtag array returns a binary mask of 
% size(m0) designed to filter data bad edges using match point density
m = varargin{1};

%defaults
%2m
n =21; % data density kernel
Pmin=.8; % data density threshold for masking
Amin=1000; % minimum data cluster area
cf = 0.5; %boundary curvature factor (0= point boundary, 1 =conv hull)
crop=n;

%for 2m SETSM 2.xxxx
% n=101;
% Pmin=.99;



% %8m
% n =5; % data density kernel
% Pmin=.8; % data density threshold for masking
% Amin=250; % minimum data cluster area
% cf = 0.5; %boundary curvature factor (0= point boundary, 1 =conv hull)
% crop = n;

% parse inputs
for i = 2:2:length(varargin)
    switch lower(varargin{i})
        case 'n'
            n=varargin{i+1};
        case 'pmin'
            Pmin=varargin{i+1};
        case 'amin'
            Amin=varargin{i+1};
        case 'crop'
            crop=varargin{i+1};
        case 'cf'
            cf=varargin{i+1};
    end
end


P = DataDensityMap(m,n); % data density map is the fraction of pixels in
                         % the kernel containing data points.

M=P>=Pmin; % data density threshold for masking

if ~any(M(:)); return; end

M = imfill(M,'holes'); % fill interior holes since we're just looking for
                       % edges here.
M = bwareaopen(M, Amin); % get rid of isolated little clusters of data

if ~any(M(:)); return; end

B = bwboundaries(M, 8, 'noholes'); % find data coverage boundaries

B = cell2mat(B); % merge all the data clusters

%k = convhull(B(:,2),B(:,1)); % find outer data boundary - the convex hull
                              % gives the straight line trace of the outer 
                              % edge.
                              
k = boundary(B(:,2),B(:,1),cf); % find outer data boundary - this function
                                % is like convhull but allows for some
                                % inward "bending" to conform to the data.

M = poly2mask(B(k,2),B(k,1), size(m,1),size(m,2)); % build edge mask

M = imerode(M,ones(crop)); % apply crop









