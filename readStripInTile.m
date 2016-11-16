function [xi,yi,z,mt,or,c,r] = readStripInTile(metaFile,x,y,varargin)

% parse input args
for i=1:2:length(varargin)
    
    switch lower(varargin{i})    
        case 'mask'
            
            mask=varargin{i+1};
            
            if ~iscell(mask) || size(mask,2) ~= 2
                error('input variable "mask" must be a cell array with two columns (x,y)')
            end

        otherwise
            
            error('Unknown input arguments')
    end
end


% readStrip read the strip and crop to a defined grid
buff=0;

% data file names
demFile= strrep(metaFile,'meta.txt','dem.tif');
matchFile= strrep(metaFile,'meta.txt','matchtag.tif');
orthoFile= strrep(metaFile,'meta.txt','ortho.tif');

% read dem first
z=readGeotiff(demFile,'map_subset',[min(x),max(x),min(y),max(y)]);

xi=z.x;
yi=z.y(:);
z=z.z;

% reset noData to NaNs
z(z < -100 | z == 0 | z == -NaN ) = NaN;

% check for blank DEM
M = ~isnan(z);


%% Apply mask if exists
if any(M(:)) && exist('mask','var')
    
    fprintf('applying mask\n')
    for i=1:size(mask,1)
        if ~isempty(mask{i,1}) && ~isempty(mask{i,2})
            M(roipoly(xi,yi,M,mask{i,1},mask{i,2}))=0;
        end
    end
    z(~M) = nan;
end

if ~any(M(:))
    xi=[];
    yi=[];
    z=[];
    mt=[];
    or=[];
    c = [];
    r = [];
    
    
    fprintf('no non-NaN data\n');
    return; 
end

% get boundaries

rowsum = sum(M) ~= 0;
colsum = sum(M,2) ~= 0;

% crop the nans
r = [];
c = [];

c(1) = find(rowsum,1,'first')-buff;
c(2) = find(rowsum,1,'last')+buff;

r(1) = find(colsum,1,'first')-buff;
r(2) = find(colsum,1,'last')+buff;

if c(1) < 1; c(1)=1; end
if r(1) < 1; r(1)=1; end

sz=size(z);

if c(2) > sz(2); c(2)=sz(2); end
if r(2) > sz(1); r(2)=sz(1); end

z = z(r(1):r(2),c(1):c(2));
xi = xi(c(1):c(2));
yi = yi(r(1):r(2));

mt=readGeotiff(matchFile,'map_subset',[min(x),max(x),min(y),max(y)]);
mt = mt.z(r(1):r(2),c(1):c(2));

or=readGeotiff(orthoFile,'map_subset',[min(x),max(x),min(y),max(y)]);
or = or.z(r(1):r(2),c(1):c(2));

% crop mosic to this DEM's boundary
c = x >= min(xi) & x <= max(xi);
r = y >= min(yi) & y <= max(yi);

if diff(x(1:2)) ~= diff(xi(1:2));
    
    z=interp2(xi,yi,z,x(c),y(r),'*linear');
    
    mt=single(mt); mt(mt == 0)  =NaN;
    mt=interp2(xi,yi,mt,x(c),y(r),'*nearest');
    mt(isnan(mt)) = 0;
    mt=logical(mt);
    
    or=single(or); or(or == 0)  =NaN;
    or=interp2(xi,yi,or,x(c),y(r),'*nearest');
    or(isnan(or)) = 0;
    or=int16(or);
    
    xi = x(c); yi=y(r);
    
end







