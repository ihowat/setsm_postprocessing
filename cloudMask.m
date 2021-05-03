function M = cloudMask(z,or,P)
% cloudMask mask bad surfaces on DEM based on slope and radiance

% M = cloudMask(z,or)masks bad edge values in the DEM with coordinate
% vectors x any y.
%
% Ian Howat, ihowat@gmail.com
% 24-Jul-2017 15:39:07

%parameters
avg_kernel_size = 21;%size of kernel for calculating mean slope
dilate_bad = 21;%dialates masked pixels by this number of surrounding pixels
min_data_cluster = 1e4;%isolated clusters of data smaller than this number will be removed
PMinCloud = 0.9;
min_nodata_cluster=1000;

% make sure sufficient non NaN pixels exist, otherwise cut to the chase
if sum(~isnan(z(:))) < 2.*min_nodata_cluster
    M=true(size(z));
    return
end

% calculate standard deveiation of elevation
avg_kernel = ones(avg_kernel_size)/(avg_kernel_size .^2);
Mz =conv2(z,avg_kernel,'same');
Sz = sqrt( conv2(z.^2,avg_kernel,'same') - Mz.^2);
Sz = real(Sz);

% calculate elevation percentile difference
dprctile = prctile(z(:),80) - prctile(z(:),20);

% set standard deviation difference based on percentile difference
if dprctile <= 40
    SzThresh = 10.5;
elseif dprctile > 40 && dprctile  <= 50
    SzThresh = 15;
elseif dprctile > 50 && dprctile  <= 75
    SzThresh = 19;
elseif dprctile > 75 && dprctile  <= 100
    SzThresh = 27;
elseif dprctile > 100
    SzThresh = 50;
end

fprintf('20/80 prctile elevation difference: %.1f, sigma-z threshold: %.1f\n',dprctile, SzThresh)

% apply mask conditions

% 3e version - refined SzThresh's and lowered Pmin to 0.9
 M = ~isnan(z) & ( ...
    ( or > 70  & P < PMinCloud ) |...
    P  < 0.6 |... 
    Sz > SzThresh );

% fill holes in masked clusters
M = imfill(M,'holes');

% remove small masked clusters
M = bwareaopen(M,min_nodata_cluster);

% remove thin borders cauased by cliffs/ridges
Me = imerode(M,ones(31));
Me = imdilate(Me,ones(61));

M = M & Me;

% dilate no data
M = imdilate(M,ones(dilate_bad));

% remove small clusters of unfiltered data
M = ~bwareaopen(~M,min_data_cluster);
