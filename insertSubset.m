function insertSubset(baseFile,insertFile,xr,yr,p,outName)
% insertSubset: overlay a subset from one tile onto another, including
% metadata
%
% insertSubset(baseFile,insertFile,xr,yr,p,outName) places the subset defined by xr =
% [min(x) max(x)] and yr = [min(y) max(y)] and polygon p from the insertFile into
% baseFile. Removed median offset around subset boundaries and applies edge
% blending. Currently supports mat files only. Appends results to outName. 
% If p=[], p is set to the polygon [xr,yr].

source_dzmed_dist_meters = 30;
feather_dist_meters = 500;

%% if no p given, make p vertices the xr yr rectangle boundary
%if isempty(p)
%    p = [[xr(1);xr(1);xr(2);xr(2);xr(1)],...
%         [yr(1);yr(2);yr(2);yr(1);yr(1)]];
%end

%% load insert file into mat object
%m1=matfile('~/Desktop/arcticdem_subset_mosaics/arcticdem_subset_mosaics_2020_m05/18_39a_10m.mat');
m1=matfile(insertFile);

dx = abs(m1.x(1,2)-m1.x(1,1));

% % extract coordinate arrays
% x1=m1.x;
% y1=m1.y;

% cols and rows of the insert DEM in the range
subcols = find(m1.x >= min(xr) & m1.x <= max(xr));
subrows = find(m1.y >= min(yr) & m1.y <= max(yr));
%
% % load subset of insert DEM
% z1 = m1.z(subrows,subcols);
x1=m1.x(1,subcols);
y1=m1.y(1,subrows);
y1=y1(:); % make vertical for interp
%
% %% crop all nan cols/rows
% nancols =  ~any(~isnan(z1));
% z1(:,nancols) = [];
% x1(nancols)=[];
%
% nanrows =  ~any(~isnan(z1'));
% z1(nanrows,:) = [];
% y1(nanrows)=[];

%% Load full mosaic as mat obj

if strcmp(baseFile,outName)
    m = matfile(baseFile ,'Writable',true);
    overwriteFlag=true;
else
    m=matfile(baseFile);
    overwriteFlag=false;
end

% load coordintate vectors
x=m.x;
y=m.y;
y=y(:);

% find extent of subset in mosaic
cols = find(x1(1) == m.x) :  find(x1(end) == m.x);
rows = find(y1(1) == m.y) :  find(y1(end) == m.y);

% below needed if insert/base not on same grid
% cols = x >= x1(1) &  x <= x1(end);
% rows = y >= y1(end) &  y <= y1(1);

%% Register and prepare insert
% Difference insert and base over subset
z1 = m1.z(subrows,subcols);
dz = z1 - m.z(rows,cols);

buff_px = ceil(source_dzmed_dist_meters / dx);
if ~isempty(p)
    % apply polygon mask
    BW0 = roipoly(x1,y1,dz,p(:,1),p(:,2));
    BW1 = imdilate(BW0,ones(buff_px*2+1));
    BW1(BW0) = false;
else
    BW0 = ones(size(dz),'logical');
    BW1 = zeros(size(dz),'logical');
    BW1(1,  :) = 1;
    BW1(end,:) = 1;
    BW1(:,  1) = 1;
    BW1(:,end) = 1;
    BW1 = imdilate(BW1,ones(buff_px*2+1));
end

% Median of differences along edges
%dzmed = nanmedian([dz(1,:),dz(end,:),dz(:,1)',dz(:,end)']);
dzmed = nanmedian(dz(BW1));

% subtract dzed offset from insert
z1 = z1 - dzmed;

% set NaNs in insert to zero
z1(isnan(z1)) = 0;

%% make an edge feathering array

% % this code can efficiently handle cases of no NaN values around border.
% A = ones(size(dz));
% C = (1:100)./100;
% A(:,1:length(C)) = repmat(C,size(dz,1),1);
% A(:,end-length(C)+1:end) = repmat(fliplr(C),size(dz,1),1);
%
% B = ones(size(dz));
% C = C(:);
% B(1:length(C),:) = repmat(C,1,size(dz,2));
% B(end-length(C)+1:end,:) = repmat(flipud(C),1,size(dz,2));
%
% A(B < A) = B(B<A);

% This code can handle nans around border but needs inpaint_nan function
% A = ~isnan(dz);
% A = padarray(A,[1,1],0,'both');
% B= imerode(A,ones(100));
% A = double(A);
% A(A~=B) = NaN;
% clear B
% A=inpaint_nans(A,2);
% A=A(2:end-1,2:end-1);


% This code can handle nans around border but needs inpaint_nan function
buff_px = ceil(feather_dist_meters / dx);
A = ~isnan(dz) & BW0;
A = padarray(A,[1,1],0,'both');
B= imerode(A,ones(buff_px*2+1));
A = double(A);
A(A~=B) = NaN;
clear B
A=inpaint_nans(A,2);
A=A(2:end-1,2:end-1);

A(A < 0) = 0;
A(A > 1) = 1;
A_interp = (A > 0) & (A < 1);

%% Load base dem and merge insert
if overwriteFlag
    m.z(rows,cols) = z1.*A + ...
        m.z(rows,cols).*(1-A);

    clear z1

    % merge z_mad
    m.z_mad(rows,cols) = m1.z_mad(subrows,subcols).*A + ...
        m.z_mad(rows,cols).*(1-A);

    arr1 = m1.N(subrows,subcols);
    arr2 = m.N(rows,cols);
    I = uint8(single(arr1).*A + single(arr2).*(1-A));
    I_interp = max(arr1, arr2);
    I(A_interp) = I_interp(A_interp);
    I(arr1 == 0) = arr2(arr1 == 0);
    I(arr2 == 0) = arr1(arr2 == 0);
    m.N(rows,cols) = I;

    arr1 = m1.Nmt(subrows,subcols);
    arr2 = m.Nmt(rows,cols);
    I = uint8(single(arr1).*A + single(arr2).*(1-A));
    I_interp = max(arr1, arr2);
    I(A_interp) = I_interp(A_interp);
    I(arr1 == 0) = arr2(arr1 == 0);
    I(arr2 == 0) = arr1(arr2 == 0);
    m.Nmt(rows,cols) = I;

    arr1 = m1.tmax(subrows,subcols);
    arr2 = m.tmax(rows,cols);
    I = uint16(single(arr1).*A + single(arr2).*(1-A));
    I_interp = max(arr1, arr2);
    I(A_interp) = I_interp(A_interp);
    I(arr1 == 0) = arr2(arr1 == 0);
    I(arr2 == 0) = arr1(arr2 == 0);
    m.tmax(rows,cols) = I;

    arr1 = m1.tmin(subrows,subcols);
    arr2 = m.tmin(rows,cols);
    I = uint16(single(arr1).*A + single(arr2).*(1-A));
    I_interp = min(arr1, arr2);
    I(A_interp) = I_interp(A_interp);
    I(arr1 == 0) = arr2(arr1 == 0);
    I(arr2 == 0) = arr1(arr2 == 0);
    m.tmin(rows,cols) = I;

    clear I;
    
    
else
    
    z = m.z;
    
    % merge with edge blending
    z(rows,cols) = z1.*A + ...
        z(rows,cols).*(1-A);
    
    save(outName,'-append','x','y','z');
    clear z z1
    
    % merge z_mad
    z_mad1=m1.z_mad(rows,cols); % load subset insert
    z_mad = m.z_mad;
    z_mad(rows,cols) = z_mad1.*A + ...
        z_mad(rows,cols).*(1-A);
    
    save(outName,'-append','z_mad');
    clear z_mad z_mad1
    
    % merge N
    N1=m1.N(rows,cols);
    N=m.N;
    N(rows,cols) = N1;
    
    save(outName,'-append','N');
    clear N1 N
    
    % merge Nmt
    Nmt1=m1.Nmt(rows,cols);
    Nmt=m.Nmt;
    Nmt(rows,cols) = Nmt1;
    
    save(outName,'-append','Nmt');
    clear Nmt1 Nmt
    
    % merge tmax and tmin
    tmax1=m1.tmax(rows,cols);
    tmax=m.tmax;
    tmax(rows,cols) = tmax1;
    
    save(outName,'-append','tmax');
    clear tmax1 tmax
    
    tmin1=m1.tmin(rows,cols);
    tmin=m.tmin;
    tmin(rows,cols) = tmin1;
    
    save(outName,'-append','tmin');
    clear tmin1 tmin
end


