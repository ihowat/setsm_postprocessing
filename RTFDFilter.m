function M = RTFDFilter(z,zr)
% Residual Topography Fractional Difference Filter
%
% M = RTFDFilter(z,zr) Uses a relative of deviaton in the residual
% topography (RT,the difference between the DEM and a smoothed DEM) between
% the DEM z and a reference DEM zr to filter noise and create binary mask
% M. DEM arrays z and zr are assumed coregistered and must have equal rows
% and columms. z can be a 3-D array, in which case each z(:,:,i) is
% fitlered independently and M is the same size as z.
% 
% Requires inpaint_nans


% kernel size (NxN)
N = 5;
% connected groups of pixels less than this value will be removed
minPixClusterSize=1e4; 

%% Residual topography from reference dem

% kernel average convolution of reference
B = conv2(zr,ones(N),'same')./N^2;

% apply standard deviation filter to RT
B = stdfilt(zr-B,ones(N));

% standard deviation filter returns zeros where there are Nans, so set the
% reference std dev filter results to nan to avoid artifacts in boundary
% interpolation
B(B==0) = NaN;

% set very small reference RT values to 0.1 to prevent small denominator in
% difference ratio
B(B < 0.01) = 0.01;

%% Make filter mask for each z(:,:,i)

% intitiate output mask stack
M = true(size(z));

% loop through stack and make masks
i=1;
for i=1:size(z,3)
    
    % kernel average convolution of strip
    A = conv2(z(:,:,i),ones(N),'same')./N^2;
    
    % standard deviation filter of residual topography
    A = stdfilt(z(:,:,i)-A,ones(N));
    
    % difference of residual topography
    dAB = A-B;
    
    % select bad data as fractional difference in residual topography and
    % residual topography above a threshold
    bad = dAB./B >= 10 & dAB > 0.5;
    
    % steps below address boundary of std filter by setting out-of-bounds 
    % to zero and the interpolating across bounary pixels, setting values
    % gt 0 to bad
    
    % covert to doubles for inserting nans and inpaint_nans only takes
    % doubles
    bad = double(bad);
    
    % set filtering boundary to NaN
    bad(A ==0)=NaN;
    
    % set strip boundary to zero
    bad(isnan(z(:,:,i)))=0;
    
    % interpolate between M=1 (bad) and boundary
    bad = inpaint_nans(bad,2);

    % set logical mask as M gt zero (accounting for interp) and strip 
    % boundary
    M(:,:,i) = (bad < 0.01 | isnan(zr)) & ~isnan(z(:,:,i));
    
    % remove small clusters of data
    M = bwareaopen(M,minPixClusterSize);
end
