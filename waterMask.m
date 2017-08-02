function M = waterMask(or,meanSunElevation,P)

Jthresh=0.2;
orThresh = 20;

if meanSunElevation < 30; orThresh = 5; end

% subtraction image
or_subtraction =  movmax(or,5) - movmin(or,5);

% entropy image
J = entropyfilt(uint8(or_subtraction),ones(5));

% set edge-effected values to zero
J(or == 0) = 0;

% mask data with entopy less than threshold.
M1 =  or~=0 & J < Jthresh;

% remove isolated clusters of masked pixels
M1 = bwareaopen(M1,500);

% dialate masked pixels
M1 = imdilate(M1, ones(7));

% mask data with low radiance and matchpoint density
M2 =  or~=0 & or < orThresh & P < .98;

M2 = bwareaopen(M2,500);

M = ~M1 & ~M2 & or~=0;

% remove isolated clusters of data
M = bwareaopen(M,500);
M = ~bwareaopen(~M,500);
