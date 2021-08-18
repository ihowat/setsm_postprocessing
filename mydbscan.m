function [idx,corepts] = mydbscan(X,epsilon,minpts)
%DBSCAN Density-Based algorithm for clustering
%   IDX = DBSCAN(X, EPSILON, MINPTS) partitions the points in the N-by-P
%   data matrix X into clusters based on parameters EPSILON and MINPTS.
%   EPSILON is a threshold for a neighborhood search query. MINPTS is a
%   positive integer used as a threshold to determine whether a point is a
%   core point. IDX is an N-by-1 vector containing cluster indices. An
%   index equal to '-1' implies a noise point.
%
%   IDX = DBSCAN(D, EPSILON, MINPTS, 'DISTANCE', 'PRECOMPUTED') is an
%   alternative syntax that accepts distances D between pairs of
%   observations instead of raw data. D may be a vector or matrix as
%   computed by PDIST or PDIST2, or a more general dissimilarity vector or
%   matrix conforming to the output format of PDIST or PDIST2.
%
%   [IDX, COREPTS] = DBSCAN(...) returns a logical vector COREPTS
%   indicating indices of core-points as identified by DBSCAN.
%
%   IDX = DBSCAN(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional
%   parameter name/value pairs to control the algorithm used by DBSCAN.
%   Parameters are:
%
%   'Distance'      -   a distance metric which can be any of the distance 
%                       measures accepted by the PDIST2 function. The 
%                       default is 'euclidean'. For more information on 
%                       PDIST2 and available distances, type HELP PDIST2. 
%                       An additional choice is:
%     'precomputed' -   Needs to be specified when a custom distance matrix
%                       is passed in
%
%    'P'            -   A positive scalar indicating the exponent of Minkowski
%                       distance. This argument is only valid when 'Distance'
%                       is 'minkowski'. Default is 2.
%  
%    'Cov'          -   A positive definite matrix indicating the covariance
%                       matrix when computing the Mahalanobis distance. This
%                       argument is only valid when 'Distance' is
%                       'mahalanobis'. Default is NANCOV(X).
%  
%    'Scale'        -   A vector S containing non-negative values, with length
%                       equal to the number of columns in X. Each coordinate
%                       difference between X and a query point is scaled by the
%                       corresponding element of S. This argument is only valid
%                       when 'Distance' is 'seuclidean'. Default is NANSTD(X).
%
%   Example:
%      % Find clusters in data X, using the default distance metric 
%      % 'euclidean'.
%      X = [rand(20,2)+2; rand(20,2)];
%      idx = dbscan(X,0.5,2);
%
%   See also KMEANS, KMEDOIDS, PDIST2, PDIST.

%   Copyright 2018-2020 The MathWorks, Inc.

% Parse X
validateattributes(X,{'single','double','logical'},{'2d','nonsparse','real',...
    'nonempty'},'','X');

% Parse epsilon
validateattributes(epsilon,{'single','double'},{'real',...
    'nonsparse','nonnegative'},'','epsilon');

% Epsilon must be [] when X is logical and must be non-empty when X is
% numeric
if islogical(X) && ~isempty(epsilon)
    error(message('stats:dbscan:EpsilonEmpty')); % Epsilon must be empty
elseif isnumeric(X) && isempty(epsilon)
    error(message('stats:dbscan:EpsilonNonEmpty')); % Epsilon must not be empty
elseif ~isscalar(epsilon) && ~isempty(epsilon)
    error(message('stats:dbscan:EpsilonScalar'));
end

% Parse minpts
validateattributes(minpts,{'numeric'},{'scalar','real','nonempty',...
    'nonsparse','integer','positive'},'','minpts');

% Store size of X
[nobs,ndims] = size(X);

idx = internal.stats.dbscan(X',epsilon,minpts,'euc',2,false,nobs,ndims,[],false);
end
