function [dtrans,rmse] = coregister2tile(metaFile,m,N,varargin)
%coregrister2tile merge strip into mosaic file
% 
%   Subfunctions: readStripInTile, coregistedems,

%
%   Ian Howat, Ohio State University
%   version 

% Set Defaults

mask=cell(1,2);
minOverlapPixels=100;

dtrans=[NaN;NaN;NaN];
rmse=NaN;

% parse input args
for i=1:2:length(varargin)
    
    switch lower(varargin{i})
        
        
        
        case 'mask'
            
            mask=varargin{i+1};
            
            if ~iscell(mask) || size(mask,2) ~= 2
                error('input variable "mask" must be a cell array with two columns (x,y)')
            end
            
        case 'minoverlappixels'
            
            minOverlapPixels=varargin{i+1};
            
        otherwise
            
            error('Unknown input arguments')
    end
end


% read this strip data
[x,y,z,mt,~,c,r] = readStripInTile(metaFile,m.x,m.y,'mask',mask);

% check for competely masked file
if isempty(x)
    fprintf('no data\n',sum(c)*sum(r));
    return;
end

% convert mosaic subset area to index
c=find(c);
r=find(r);

if length(c)*length(r) < minOverlapPixels
    fprintf('%d pixel overlap is too small, skipping\n',sum(c)*sum(r));
    return;
end

%% Check for new coverage

%subset count grid
Nsub=N(r,c);


A= Nsub ~= 0 & ~isnan(z); % overlapping coverage mask


%% Registration
if sum(A(:)) < minOverlapPixels
    fprintf('%d pixel overlap is too small, skipping\n',sum(A(:)));
    return;
end

% crop to new coverage
co= find(sum(A) ~= 0,1,'first'):find(sum(A) ~= 0,1,'last');
ro= find(sum(A,2) ~= 0,1,'first'):find(sum(A,2) ~= 0,1,'last');

% co-register strip to mosaic subset
try
[~,dtrans,rmse] = coregisterdems(...
    m.x(1,c(co)),m.y(r(ro),1),m.z(r(ro),c(co)),...
    x(co),y(ro),z(ro,co),m.mt(r(ro),c(co)),mt(ro,co));
catch
    
    fprintf('coregistration failure, skipping\n');
    return;
end

% coregisterdems returns dtrans values as positive offsets, whereas
% the gcp registration are negative, so need to reverse the sign:
dtrans=-dtrans;





