function [dz_med,dz_mad,N] = diffDEMFromGCPs(m,gcp,waterMask)
% diffDEMFromGCPs vertical difference between DEM and gcp's

maxGcps=1000; % maximum number of gcps to use. Will sample gcp vector evenly to give the closest integer N > maxGcps.
footprintRadius = 30;

% initialize output
dz_med =NaN;
dz_mad =NaN;
N =0;

%subset gcps over this tile
n =find(gcp.x > min(m.x) & gcp.x < max(m.x) & gcp.y > min(m.y) & gcp.y < max(m.y));

if isempty(n)
    fprintf('no gcp''s overlap this tile, retuning NaN\n')
    return
end

% crop gcp structure
gcp = structfun( @(x) x(n),gcp,'uniformoutput',0);

% sample gcp's over water mask
landFlag = interp2(m.x,m.y,waterMask,gcp.x,gcp.y,'*nearest');

n = find(landFlag == 1);

if isempty(n)
    fprintf('all gcp''s over water, retuning NaN\n')
    return
end

%subsample GCPs to near maximum # for speed/memory
if length(n) > maxGcps
    n=n(1:floor(length(n)/ maxGcps):length(n));
end

gcp = structfun( @(x) x(n),gcp,'uniformoutput',0);

% initialize coordinate subset cells
zsub=nan(size(gcp.z));


% subsetting loop
j=1;
for j=1:length(zsub)
    minCol=find(m.x >=  gcp.x(j) - footprintRadius,1,'first');
    maxCol=find(m.x <=  gcp.x(j) + footprintRadius,1,'last');
    
    minRow=find(m.y <=  gcp.y(j) + footprintRadius,1,'first');
    maxRow=find(m.y >=  gcp.y(j) - footprintRadius,1,'last');
    
    if any(any(~waterMask(minRow:maxRow,minCol:maxCol)))
        continue
    end
    
    z = m.z(minRow:maxRow,minCol:maxCol);
    
    zsub(j)=nanmean(z(:));
end

dz = zsub - gcp.z;

dz_med = nanmedian(dz);
dz_mad = mad(dz,1);
N = length(~isnan(dz));
