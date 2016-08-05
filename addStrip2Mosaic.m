function N = addStrip2Mosaic(metaFile,m,dy0,N,dtrans,varargin)
%addStrip2Meta merge strip into mosaic file
%   N = addStrip2Mosaic(metaFile,m,dy0,N,dtrans) adds the strip data
%   specified in the *_meta.txt file into the mat file handle m. N is the
%   mosaic data/no data byte mask. dtrans are [dz;dx;dy] registration
%   offsets to be applied to the strip before merging (where z1=z0+dz). If
%   dtrans=[0;0;0], no registration will be applied. If dtrans=[], the
%   strip will be corgistered to existing data in the mosaic if exists, or
%   inserted without registration if not. 
%
%   N = addStrip2Mosaic(...,mergeMethodString) where mergeMethodString can
%   be one of:
%       'feather'   Apply a linear weighting between the edges of
%                   overlapping strips. DEFAULT
% 
%       'warp'      Force the edge of the new data to align with the edge
%                   of the mosaic and throw away new data that overlaps.
%
%       'underprint' Only add new data where there is no existing data
%
%   Subfunctions: readStripInTile, edgeFeather, edgeWarp, coregistedems, 
%                 interpolate2grid
%
%   Ian Howat, Ohio State University
%   version 1; 28-Jul-2016 09:50:32


% option to disable using control registration files
mergeMethod='feather';
if ~isempty(varargin); mergeMethod=varargin{1}; end

% read this strip data
[x,y,z,mt,or,c,r] = readStripInTile(metaFile,m.x,m.y);

% check for competely masked file
if isempty(x);
    m.dtrans=[m.dtrans,[NaN;NaN;NaN]];
    m.rmse=[m.rmse,NaN];
    return;
end

% convert mosaic subset area to index
c=find(c);
r=find(r);

if length(c)*length(r) < 1000;
    fprintf('%d pixel overlap is too small, skipping\n',sum(c)*sum(r));
    m.dtrans=[m.dtrans,[NaN;NaN;NaN]];
    m.rmse=[m.rmse,NaN];
    return;
end

%% Check for new coverage

%subset count grid
Nsub=N(r,c);

%quit if too little new pixel added
if sum(Nsub(:)==0 & ~isnan(z(:))) < 1000;
    disp('redundant coverage, skipping');
    m.dtrans=[m.dtrans,[NaN;NaN;NaN]];
    m.rmse=[m.rmse,NaN];
    return;
end

A= Nsub ~= 0 & ~isnan(z); % overlapping coverage mask

%% Registration

% preset rmse
rmse=0;

% if no dtrans given, but there's no overlap to coregister, then set 0
if ~any(A(:)) && isempty(dtrans); dtrans=[0;0;0]; end
 
% if dtrans is empty and there's overlap, then apply coregistration
if isempty(dtrans)
    
    % crop to new coverage
    co= find(sum(A) ~= 0,1,'first'):find(sum(A) ~= 0,1,'last');
    ro= find(sum(A,2) ~= 0,1,'first'):find(sum(A,2) ~= 0,1,'last');
    
    % co-register strip to mosaic subset
    [~,dtrans,rmse] = coregisterdems(...
        m.x(1,c(co)),m.y(r(ro),1),m.z(r(ro),c(co)),...
        x(co),y(ro),z(ro,co),m.mt(r(ro),c(co)),mt(ro,co));
    
    % coregisterdems returns dtrans values as positive offsets, whereas
    % the gcp registration are negative, so need to reverse the sign:
    dtrans=-dtrans;
    
    % check for coregistration failure
    if isnan(rmse); return; end;
    
end

% apply dtrans and interpolate back to grid (but only dtrans is not zeros)
if any(dtrans ~= 0)
    
    x = x + dtrans(2);
    y = y + dtrans(3);
    z = z + dtrans(1);
    
    % interpolate the first dem to the same grid
    [z,mt,or] = interpolate2grid(x,y,z,mt,or,m.x(1,c),m.y(r,1));
    
end

% put into file
m.dtrans=[m.dtrans,dtrans(:)];
m.rmse=[m.rmse,rmse];


%% Data Merge

% if no overlap, set to underprint
if ~any(A(:));  mergeMethod = 'underprint'; end

% make date grid
dy=~isnan(z).*dy0;

% for the matchtag, just straight combination
m.mt(r,c) = m.mt(r,c) | mt;
clear mt

switch mergeMethod
    
    case 'underprint'

        n = find(Nsub == 0 & ~isnan(z));
        
        zsub = m.z(r,c);
        zsub(n)=z(n);
        m.z(r,c)=zsub;
        clear z zsub

        orsub=m.or(r,c);
        orsub(n)=or(n);
        m.or(r,c) = orsub;
        clear or orsub
        
        dysub=m.dy(r,c);
        dysub(n)=dy0;
        m.dy(r,c) = dysub;
        clear dy dysub        
        
    case 'feather'
        
        
        dx = m.x(1,2)-m.x(1,1);
        buff=10*dx+1;
        
        W = edgeFeather(Nsub~=0,~isnan(z),buff);
        
        % make weighted elevation grid
        zsub=m.z(r,c);
        A=  zsub.*W + z.*(1-W);
        n = isnan(zsub) & ~isnan(z);  A(n)= z(n);
        n = ~isnan(zsub) & isnan(z);  A(n)= zsub(n);
        
        % put strip subset back into full array
        m.z(r,c) = A;
        clear A z
        
        % make weighted ortho grid
        or = single(or);
        or(or ==0) = NaN;
        orsub=single(m.or(r,c));
        orsub(orsub==0) = NaN;
        A=  orsub.*W + or.*(1-W);
        
        A( isnan(orsub) & ~isnan(or))=   or(  isnan(orsub) & ~isnan(or));
        A(~isnan(orsub) &  isnan(or))= orsub(~isnan(orsub) &  isnan(or));
        A(isnan(A)) = 0; % convert back to uint16
        A = uint16(A);
        m.or(r,c) = A;
        clear A or
        
        % make weighted dy grid
        dy = single(dy);
        dy(dy==0) = NaN;
        dysub = single(m.dy(r,c));
        dysub(dysub==0) = NaN;
        
        A=  dysub.*W + dy.*(1-W);
        
        A( isnan(dysub) & ~isnan(dy))=   dy ( isnan(dysub) & ~isnan(dy));
        A(~isnan(dysub) &  isnan(dy))= dysub(~isnan(dysub) &  isnan(dy));
        
        A(isnan(A)) = 0; % convert back to uint16
        A = uint16(A);
        m.dy(r,c) = A;
        clear A dy
        
        clear W
        
    case 'warp'
        
        % merge z by warping edge of z to match zsub
        m.z(r,c) = edgeWarp(m.z(r,c),z);
        
        orsub=m.or(r,c);
        orsub(orsub == 0)=or(orsub == 0);
        m.or(r,c) = orsub;
        clear orsub
        
        dysub=m.dy(r,c);
        dysub(dysub == 0)=dy(dysub == 0);
        m.dy(r,c) = dysub;
        clear dysub
        
end

N(r,c)=~isnan(m.z(r,c));



