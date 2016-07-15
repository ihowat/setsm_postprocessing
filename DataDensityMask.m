function M = DataDensityMask(varargin)
% edgeMask returns mask for bad edges using the matchag field
%
% m1 = edgeMask(m0) where m0 is the matchtag array returns a binary mask of 
% size(m0) designed to filter data bad edges using match point density and
% defualt parameters
%
% m1 = edgeMask(m0,Pmin,Amax) where Pmin is the minimum data density 
% threshold for masking and Amax maximum data gap area to leave filled

m = varargin{1};

% Defaults
n =21; % data density kernel
Pmin=.3; % data density threshold for masking
Amin=1000; % minimum data cluster area
Amax=10000; % maximum data gap area to leave filled

% parse inputs
for i = 2:2:length(varargin)
    switch lower(varargin{i})
        case 'n'
            n=varargin{i+1};
        case 'pmin'
            Pmin=varargin{i+1};
        case 'amin'
            Amin=varargin{i+1};
        case 'amax'
            Amax=varargin{i+1};
    end
end

P = DataDensityMap(m,n); % data density map is the fraction of pixels in
                         % the kernel containing data points.

% apply match density threshold                         
M=P>=Pmin;

if ~any(M(:)); return; end

% remove small data clusters                       
M = bwareaopen(M, Amin); % get rid of isolated little clusters of data

% remove small data gaps
M = ~bwareaopen(~M, Amax);





