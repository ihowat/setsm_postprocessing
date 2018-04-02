function [N,returnFlag] = addStrip2Mosaic(metaFile,m,dy0,N,dtrans,rmse,varargin)
%addStrip2Meta merge strip into mosaic file
%   N = addStrip2Mosaic(metaFile,m,dy0,N,dtrans) adds the strip data
%   specified in the *_meta.txt file into the mat file handle m. N is the
%   mosaic data/no data byte mask. dtrans are [dz;dx;dy] registration
%   offsets to be applied to the strip before merging (where z1=z0+dz). If
%   dtrans=[0;0;0], no registration will be applied. If dtrans=[], the
%   strip will be corgistered to existing data in the mosaic if exists, or
%   inserted without registration if not. 
%
%   N = addStrip2Mosaic(...,'option',value) where options are:
%   'mergeMethod' can be one of:
%       'feather'   Apply a linear weighting between the edges of
%                   overlapping strips. DEFAULT
% 
%       'warp'      Force the edge of the new data to align with the edge
%                   of the mosaic and throw away new data that overlaps.
%
%       'underprint' Only add new data where there is no existing data
%
%   'mask', mask where mask is an nx2 cell array of mask vertices
%   coordinates with columns x and y, and each row a polygon
%
%   'maxrmse', value specifies a maximum rmse, over which the results will
%   be rejected and nan's will be returned
%
%   Subfunctions: readStripInTile, edgeFeather, edgeWarp, coregistedems, 
%                 interpolate2grid
%
%   Ian Howat, Ohio State University
%   version 1.0; 28-Jul-2016 09:50:32
%   version 1.1; 06-Dec-2016 14:29:02
%   - changed to NaNs in dtrans triggers coregistration

% Set Defaults
mergeMethod='feather';
mask=cell(1,2);
maxrmse=inf;
minNewPixels=100;
minOverlapPixels=100;
returnFlag = false;
addErrorFlag=false;
refineRegFlag=false;
resetDate=false;

% parse input args
for i=1:2:length(varargin)
    
    switch lower(varargin{i})
        
        case 'mergemethod'
            
            mergeMethod=varargin{i+1};
            
            if ~strcmp(mergeMethod,'feather') && ~strcmp(mergeMethod,'warp') &&  ~strcmp(mergeMethod,'underprint')
                error('unrecognized merge method string')
            end
            
        case 'mask'
            
            mask=varargin{i+1};
            
            if ~iscell(mask) || size(mask,2) ~= 2
                error('input variable "mask" must be a cell array with two columns (x,y)')
            end
            
        case 'maxrmse'
            
            maxrmse=varargin{i+1};
            
        case 'minnewpixels'
            
            minNewPixels=varargin{i+1};
            
            
        case 'minoverlappixels'
            
            minOverlapPixels=varargin{i+1};
            
        case 'adderrors'
            
            addErrorFlag=varargin{i+1};
            
            if ~islogical(addErrorFlag)
                error('addError argument must be logical')
            end
            
        case 'refinereg'
            
            refineRegFlag=varargin{i+1};
            
            if ~islogical(refineRegFlag)
                error('refineReg argument must be logical')
            end
            
        case 'resetdate'
            
            resetDate=varargin{i+1};
            
            if ~islogical(refineRegFlag)
                error('resetDate argument must be logical')
            end
            
        otherwise
            
            error('Unknown input arguments')
    end
end


% read this strip data
[x,y,z,mt,or,c,r] = readStripInTile(metaFile,m.x,m.y,'mask',mask);

% check for competely masked file
if isempty(x)
    m.dtrans=[m.dtrans,[NaN;NaN;NaN]];
    m.rmse=[m.rmse,NaN];
    returnFlag=true;
    return;
end

% convert mosaic subset area to index
c=find(c);
r=find(r);

if length(c)*length(r) < minOverlapPixels
    fprintf('%d pixel overlap is too small, skipping\n',sum(c)*sum(r));
    return;
end

%% Check for new coverage

%subset count grid
Nsub=N(r,c);

%quit if too few new pixels added
newPixels=sum(Nsub(:)==0 & ~isnan(z(:)));
if newPixels < minNewPixels
    disp('redundant coverage, skipping');
    m.dtrans=[m.dtrans,[NaN;NaN;NaN]];
    m.rmse=[m.rmse,NaN];
    returnFlag=true;
    return;
end

A= Nsub ~= 0 & ~isnan(z); % overlapping coverage mask


%% Registration

% if no dtrans given, but there's no overlap to coregister, then set 0
if ~any(A(:)) && isempty(dtrans); dtrans=[0;0;0]; end
 
% if dtrans is empty or regRefineFlag and there's overlap, then apply 
% coregistration
if isempty(dtrans) || any(isnan(dtrans)) || refineRegFlag 
    
    if sum(A(:)) < minOverlapPixels
        fprintf('%d pixel overlap is too small, skipping\n',sum(A(:)));
        return;
    end

    % crop to new coverage
    co= find(sum(A) ~= 0,1,'first'):find(sum(A) ~= 0,1,'last');
    ro= find(sum(A,2) ~= 0,1,'first'):find(sum(A,2) ~= 0,1,'last');
    
    % set dtrans to zeros if this is a new adjustment
    if ~refineRegFlag
        dtrans0 = zeros(3,1); 
    else
        dtrans0=dtrans;
    end
    
    % co-register strip to mosaic subset
    [~,dtrans,rmse] = coregisterdems(...
        m.x(1,c(co)),m.y(r(ro),1),m.z(r(ro),c(co)),...
        x(co)+dtrans0(2),y(ro)+dtrans0(3),z(ro,co)+dtrans0(1),...
        m.mt(r(ro),c(co)),mt(ro,co));
    
    % coregisterdems returns dtrans values as positive offsets, whereas
    % the gcp registration are negative, so need to reverse the sign:
    dtrans=-dtrans;
    
    dtrans = dtrans(:) + dtrans0(:);
    
    % check for coregistration failure
    if isnan(rmse) || rmse > maxrmse
        
        disp('coregistration failure, skipping');
        m.dtrans=[m.dtrans,[NaN;NaN;NaN]];
        m.rmse=[m.rmse,NaN];
        
        return
    end    
    
    coregFlag=true;
    
end

% apply dtrans and interpolate back to grid (but only dtrans is not zeros)
if any(dtrans ~= 0)
     
        x0=x;
        y0=y;
     
        x = x + dtrans(2);
        y = y + dtrans(3);
        z = z + dtrans(1);

        % check for precision issues that result in a non-uniform grid
        if range(x(2:end)-x(1:end-1)) ~= 0
            fprintf('non-uniform x-spacing detected following dtrans apply, correcting\n')
            % force a uniform spacing by adding a vector of the cumultive sum of the target resolution to the first value of the coordinate vector. 
            x = x(1) + [0,cumsum(ones(1,length(x)-1).*(m.x(1,2)-m.x(1,1)))];  
        end
        if range(y(2:end)-y(1:end-1)) ~= 0
            fprintf('non-uniform y-spacing detected following dtrans apply, correcting\n')
            % same as for x
            y = y(1) + [0,cumsum(ones(1,length(y)-1).*(m.y(2,1)-m.y(1,1)))]';
        end

    try
        % interpolate the first dem to the same grid
        [z,mt,or] = interpolate2grid(x,y,z,mt,or,m.x(1,c),m.y(r,1));
    
    catch
        fprintf('stupid dtrans precision error, taking a dump\n')
   
        save dtrans_precision_error_dump.mat
   end
   
   clear x0 y0
    
end

%% Data Merge

% if no overlap, set to underprint
if ~any(A(:));  mergeMethod = 'underprint'; end

switch mergeMethod
    
    case 'underprint'

        n = find(Nsub == 0 & ~isnan(z));
        
        zsub = m.z(r,c);
        zsub(n)=z(n);
        m.z(r,c)=zsub;
        clear zsub
        
        vars = whos(m);
        if ismember('ze', {vars.name})
            zesub=m.ze(r,c);
            if addErrorFlag
                nn = ~isnan(zesub) & ~isnan(z);
                zesub(n)= sqrt(rmse.^2 + nanmean(zesub(nn)).^2);
            else
                zesub(n)=rmse;
            end
            m.ze(r,c) = zesub;
            clear ze zesub
        end
        clear z

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
   
        W = edgeFeather(Nsub~=0,~isnan(z),'buffer',buff);
        
        % make weighted elevation grid
        zsub=m.z(r,c);
        A=  zsub.*W + z.*(1-W);
        n = isnan(zsub) & ~isnan(z);  A(n)= z(n);
        n = ~isnan(zsub) & isnan(z);  A(n)= zsub(n);
        
        % put strip subset back into full array
        m.z(r,c) = A;
        clear A zsub
        
        % ze
        vars = whos(m);
        if ismember('ze', {vars.name})
            
             zesub=m.ze(r,c);
             
             % make error grid
             if addErrorFlag
                 nn = ~isnan(zesub) & ~isnan(z);
                 ze=~isnan(z).*sqrt(rmse.^2 + nanmean(zesub(nn)).^2);
             else
                 ze=~isnan(z).*rmse;
             end
             ze(isnan(z)) = NaN; 
             
            A=  zesub.*W + ze.*(1-W);
            
            A( isnan(zesub) & ~isnan(ze))=   ze ( isnan(zesub) & ~isnan(ze));
            A(~isnan(zesub) &  isnan(ze))= zesub(~isnan(zesub) &  isnan(ze));
            
            m.ze(r,c) = A;
            clear A ze zesub
        end

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
        clear A or orsub
        
        
        % make weighted dy grid
        dysub = single(m.dy(r,c));
        dysub(dysub==0) = NaN;
        
        if resetDate;  dy0 = nanmean(dysub(~isnan(z))); end
        
        dy=~isnan(z).*dy0;
        dy = single(dy);
        dy(dy==0) = NaN;
 
        
        A=  dysub.*W + dy.*(1-W);
        
        A( isnan(dysub) & ~isnan(dy))=   dy ( isnan(dysub) & ~isnan(dy));
        A(~isnan(dysub) &  isnan(dy))= dysub(~isnan(dysub) &  isnan(dy));
        
        A(isnan(A)) = 0; % convert back to uint16
        A = uint16(A);
        m.dy(r,c) = A;
        clear A dy dysub
        
        
        clear W
        
    case 'warp'
        
        % merge z by warping edge of z to match zsub
         [m.z(r,c),edgeWarpflag] = edgeWarp(m.z(r,c),z);
        if edgeWarpflag
            disp('redundant coverage, skipping');
            m.dtrans=[m.dtrans,[NaN;NaN;NaN]];
            m.rmse=[m.rmse,NaN];
            returnFlag=true;
            return;
        end

        n = find(Nsub == 0 & ~isnan(z));
        vars = whos(m);
        if ismember('ze', {vars.name})
            
            zesub=m.ze(r,c);
            if addErrorFlag
                nn = ~isnan(zesub) & ~isnan(z);
                zesub(n)= sqrt(rmse.^2 + nanmean(zesub(nn)).^2);
            else
                zesub(n)=rmse;
            end
            m.ze(r,c) = zesub;
            clear ze zesub
        end

        clear z
        
        orsub=m.or(r,c);
        orsub(n)=or(n);
        m.or(r,c) = orsub;
        clear or orsub
        
        dysub=m.dy(r,c);
        dysub(n)=dy0;
        m.dy(r,c) = dysub;
        clear dy dysub        
        
end

% for the matchtag, alwats just straight combination
m.mt(r,c) = m.mt(r,c) | mt;
clear mt

% put coregistration data into file
m.dtrans=[m.dtrans,dtrans(:)];
m.rmse=[m.rmse,rmse];


N(r,c)=~isnan(m.z(r,c));
returnFlag=true;
