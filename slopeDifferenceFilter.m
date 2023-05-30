function M = slopeDifferenceFilter(x,y,z,zr,C,avoidFilteringWaterFlag)
% slopeDifferenceFilter: mask using the fractional difference in roughness
%

if ~isempty(C)
    water_mask = (C.z == 0 | C.z == 80);
    snowice_mask = (C.z == 70);
else
    water_mask = [];
    snowice_mask = [];
end

dx = x(2)-x(1);

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

if ~isempty(water_mask) && avoidFilteringWaterFlag
    % Define a buffer zone near the edge of water in the water mask.
    water_edge_px = ceil(250 / dx);
    water_edge_zone = xor(...
        imdilate(water_mask, circleMaskSE(water_edge_px)),...
        imerode( water_mask, circleMaskSE(water_edge_px))...
    );

    % Don't use fracional slope stdev as an indicator close to the
    % edge of water where the reference surface is very smooth
    % (copernicus has a very smooth water fill).
    F((water_mask | water_edge_zone) & zrslopeSDF < 0.02) = 0;
end

dz = z-zr;
dz = dz - median(dz(:),'omitnan');

% make output mask
M = nan(size(F),'single');

% define point with high likliood of being good data
M(F < FThreshGood & abs(dz) < dzThreshGood) = 1;

% define points with highly liklihood of being bad data
M(F >= FThreshBad | abs(dz) >= dzThreshBad) = 0;

if ~isempty(snowice_mask)
    % Don't apply filter over snow/ice
    buffer_px = ceil(250 / dx);
    M(imdilate(snowice_mask, circleMaskSE(buffer_px))) = 1;
end

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
if ~isempty(water_mask)
    % Don't allow removal of small chunks of ones ("good" data)
    % where those chunks are in small islands (as determined from water mask).
    islands = imdilate(xor(water_mask, ~bwareaopen(~water_mask, 2000)), circleMaskSE(10));
    remove_good_data_chunks = xor(M, bwareaopen(M,1000));
    remove_good_data_chunks(islands) = 0;
    M(remove_good_data_chunks) = 0;
    M = ~bwareaopen(~M,100);
else
    M = bwareaopen(M,1000);
    M = ~bwareaopen(~M,100);
end

if ~isempty(water_mask) && avoidFilteringWaterFlag
    % Only keep pixels marked as "bad" near the edge of water
    % that are close enough to pixels marked as bad further inland.
    M = ~M;
    M_old = M;
    M_new = M;
    M_new(water_edge_zone | water_mask) = 0;
    dilate_px = ceil(700 / dx);
    M_new = M_old & imdilate(M_new, circleMaskSE(dilate_px));
    M = ~M_new;
end
