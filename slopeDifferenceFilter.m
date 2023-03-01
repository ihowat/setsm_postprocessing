function M = slopeDifferenceFilter(x,y,z,zr)
% slopeDifferenceFilter: mask using the fractional difference in roughness
%

% kernel over which to calculate slope standard deviation
kernel = 5;

% threshold in the fractional slope standard deviation
FThreshGood  = 10; % max for good data
FThreshBad = 50;  % min for bad data

% threshold in absolute height difference from reference
dzThreshGood = 10;% max for good data
dzThreshBad = 100;% min for bad data

% calculate SD of slope (roughness) over kernel
[zsx,zsy] = gradient(z,x,y);
zslope = sqrt(zsx.^2 + zsy.^2);
zslopeSDF = stdfilt(zslope,ones(kernel));

[zrsx,zrsy] = gradient(zr,x,y);
zrslope = sqrt(zrsx.^2 + zrsy.^2);
zrslopeSDF = stdfilt(zrslope,ones(kernel));

% fractional change in roughness
F = (zslopeSDF - zrslopeSDF)./zrslopeSDF;

% zero out very low slope, slope sdf areas
F(zslopeSDF < 0.02 & zrslopeSDF < 0.02) = 0;

dz = z-zr;
dz = dz - median(dz(:),'omitnan');

% make output mask
M = nan(size(F),'single');

% define point with high likliood of being good data
M(F < FThreshGood & abs(dz) < dzThreshGood) = 1;

% define points with highly liklihood of being bad data
M(F >= FThreshBad | abs(dz) >= dzThreshBad) = 0;

% interpolate voids, so that voids with no bad data will be ones and voids
% with bad data will have values < 1. There's got to be a better/faster way
% of doing this - maybe using bwconnect
M = inpaint_nans(double(M),2);

% make binary mask, set all M's less than 1 (by 2 sigfigs) as bad
M = M > 0.99;

%M = F < FThresh & abs(z-zr) < dzThresh;

% remove borders
M = imopen(M,ones(kernel));

% remove small clusters of ones and zeroes
M = bwareaopen(M,1000);
M = ~bwareaopen(~M,100);

