function [X,Y,Z,M,O,trans,rmse,f]=scenes2strips(varargin)
%SCENES2STRIPS merge scenes into strips
%
%   [x,y,z,m,o,trans,rmse,f]=scenes2strips(demdir,f) merges the
%   scene geotiffs listed in cellstr f within directory demdir after
%   ordering them by position. If a break in coverage is detected between
%   scene n and n+1 only the first 1:n scenes will be merged. The data are
%   coregistered at overlaps using iterative least squares, starting with
%   scene n=1.
%   Outputs are the strip grid coorinates x,y and strip elevation, z, 
%   matchtag, m and orthoimage, o. The 3D translations are given in 3xn 
%   vector trans, along with root-mean-squared of residuals, rmse. The
%   output f gives the list of filenames in the mosaic. If a break is
%   detected, the list of output files will be less than the input.
%
%   [...]=scenes2strips(...,'maskFileSuffix',value) will apply the mask
%   identified as the dem filename with the _dem.tif replaced by
%   _maskFileSuffix.tif
%   [...]=scenes2strips(...,'max_coreg_rmse',value) will set a new maximum
%   coregistration error limit in meters (default=1). Errors above this
%   limit will result in a segment break.
%
% Version 3.1, Ian Howat, Ohio State University, 2015.

max_coreg_rmse = 1;%meters; coregistration error larger than this will
                   %cause a segment break.
                   
         
%% Parse argins
demdir=varargin{1};
f=varargin{2};

n=find(strcmpi(varargin,'max_coreg_rmse')); % set this max coreg limit
if ~isempty(n); max_coreg_rmse=varargin{n+1}; end;

n=find(strcmpi(varargin,'maskFileSuffix')); % set this max coreg limit
if ~isempty(n); maskFileSuffix=varargin{n+1}; end;

%% Order Scenes in north-south or east-west direction by aspect ratio
fprintf('ordering %d scenes\n',length(f))

f = orderPairs(demdir,f);

% intialize output stats
trans=zeros(6,length(f));
rmse=zeros(1,length(f));

% file loop
for i=1:length(f)
    
    % construct filenames
    demFile = [demdir,'/',f{i}];
    matchFile= strrep(demFile,'dem.tif','matchtag.tif');
    orthoFile= strrep(demFile,'dem.tif','ortho.tif');
    maskFile= [];
    if exist('maskFileSuffix','var')
        maskFile= strrep(demFile,'dem.tif',[maskFileSuffix,'.tif']);
    else
        fprintf('No Mask Applied');
    end
    
    fprintf('scene %d of %d: %s\n',i,length(f),demFile)
    
    try
        [x,y,z,o,m,md] = loaddata(demFile,matchFile,orthoFile,maskFile);
    catch
        fprintf('data read error, skipping \n'); 
        continue; 
    end
                         
    % check for no data
    if  ~any(md(:)); fprintf('no data, skipping \n'); continue; end; 
    
    %Apply Masks
    [x,y,z,o,m] =  applyMasks(x,y,z,o,m,md);
     
    dx = x(2)-x(1);
    dy = y(2)-y(1);
    
    % check for non-integer grid
    if rem(x(1)./dx,1) ~= 0 || rem(y(1)./dy,1)  ~= 0
        [x,y,z,m,o] = regrid(x,y,z,m,o);
    end
    
    % if first scene in strip, set as strip and continue to next scene
    if ~exist('X','var') 
        X=x; clear x;
        Y=y; clear y;
        Z=z; clear z;
        M=m; clear m;
        O=o; clear o;
        
        continue;
    end
    
    % pad new arrays to stabilize interpolation
     buff=10*dx+1;
     z=padarray(z,[buff buff],NaN);
     m=padarray(m,[buff buff]);
     o=padarray(o,[buff buff]);
     x=[x(1)-dx.*(buff:-1:1),x,x(end)+dx.*(1:buff)];
     y=[y(1)+dx.*(buff:-1:1),y,y(end)-dx.*(1:buff)];
     
    % expand strip coverage to encompass new scene
    if x(1) < X(1)
        X1=x(1):dx:X(1)-dx;
        X=[X1,X];
        Z=[nan(size(Z,1),length(X1)),Z];
        M=[false(size(M,1),length(X1)),M];
        O=[zeros(size(O,1),length(X1)),O];
        clear X1
    end
    if x(end) > X(end)
        X1=X(end)+dx:dx:x(end);
        X=[X,X1];
        Z=[Z,nan(size(Z,1),length(X1))];
        M=[M,false(size(M,1),length(X1))];
        O=[O,zeros(size(O,1),length(X1))];
        clear X1
    end
    if y(1) > Y(1)
        Y1=y(1):-dx:Y(1)+dx;
        Y=[Y1,Y];
        Z=[nan(length(Y1),size(Z,2));Z];
        M=[false(length(Y1),size(M,2));M];
        O=[zeros(length(Y1),size(O,2));O];
        clear Y1
    end
    if y(end) < Y(end)
        Y1=Y(end)-dx:-dx:y(end);
        Y=[Y,Y1];
        Z=[Z;nan(length(Y1),size(Z,2))];
        M=[M;false(length(Y1),size(M,2))];
        O=[O;zeros(length(Y1),size(O,2))];
        clear Y1;
    end
    
    % map new dem pixels to swath. These must return integers. If not, 
    % interpolation will be required, which is currently not supported.
    c0 = find(x(1) == X);
    c1 = find(x(end) == X);
    r0 = find(y(1) == Y);
    r1 = find(y(end) == Y);
    
    %crop to overlap
    Zsub=Z(r0:r1,c0:c1);
    Xsub=X(c0:c1);
    Ysub=Y(r0:r1);
    Msub=M(r0:r1,c0:c1);
    Osub=O(r0:r1,c0:c1);
    
    %% NEW MOSAICING CODE
    
    % crop to just region of overlap
    A= single( ~isnan(Zsub) & ~isnan(z));
    
    %check for segment break
    if sum(A(:)) <= 1000
        f=f(1:i-1);
        trans = trans(:,1:i-1);
        rmse  = rmse(1:i-1);
        break
    end
    
    A(A==0)=NaN;
    
    [~,r,c] = cropnans(A,buff);
    
    %Make overlap mask removing isolated pixels
    A=single(bwareaopen(...
        isnan(Zsub(r(1):r(2),c(1):c(2))) & ~isnan(z(r(1):r(2),c(1):c(2))),...
        1000)); % nodata in strip and data in scene is a one
    
    % check for redundant scene
    if sum(A(:)) <= 1000; fprintf('redundant scene, skipping \n'); continue; end; 

    A(bwareaopen(...
        ~isnan(Zsub(r(1):r(2),c(1):c(2))) &  isnan(z(r(1):r(2),c(1):c(2))),...
        1000))=2;
         % data in strip and no data in scene is a two
    
%     tic
%     % USING REGIONFILL - Requires matlab 2015a or newer
%     A = regionfill(A,A==0) -1;
%     toc
    
    Ar=imresize(A,.1,'nearest');
    
    [C, R] = meshgrid(1:size(Ar,2),1:size(Ar,1));
    
    % pixles on outside of boundary of overlap region
    B = bwboundaries(Ar~=0, 8, 'noholes');
    B = cell2mat(B);

    n=sub2ind(size(Ar),B(:,1),B(:,2));

    warning off
    F = scatteredInterpolant(C(n),R(n),double(Ar(n)));
    warning on    

    try
        Ar(Ar==0)=F(C(Ar==0),R(Ar==0));
    catch ME
        % segment break
        fprintf('Error: %s\n', ME.message);
        fprintf('Due to this error, breaking segment\n');
        f=f(1:i-1);
        trans = trans(:,1:i-1);
        rmse  = rmse(1:i-1);
        break
    end

    Ar=imresize(Ar,size(A),'bilinear');
    Ar(A==1 & Ar ~=1)=1;
    Ar(A==2 & Ar ~=2)=2;
    A=Ar-1;
    A(A < 0) = 0;
    A(A > 1) = 1;
   
    W=single(~isnan(Zsub));
    W(r(1):r(2),c(1):c(2)) = A;  clear A
    W(isnan(Zsub) & isnan(z)) = NaN;

   
    % shift weights so that more of the reference layer is kept
    f0=.25; % overlap fraction where ref z weight goes to zero
    f1=.55; % overlap fraction where ref z weight goes to one
    
    W=(1/(f1-f0)).*W-f0/(f1-f0);
    W(W > 1) = 1;
    W(W < 0) = 0;
  
    
    % remove <25% edge of coverage from each in pair
    Zsub(W == 0) = NaN;
    Msub(W == 0) = 0;
    Osub(W == 0) = 0;

    z(W >= 1) = NaN;
    m(W >= 1) = 0;
    o(W >= 1) = 0;
 
    %% Coregistration
 
    P0 = DataDensityMap(Msub(r(1):r(2),c(1):c(2))) > 0.9;

    %check for segment break
    if ~any(P0(:))
        f=f(1:i-1);
        trans = trans(:,1:i-1);
        rmse  = rmse(1:i-1);
        break
    end

    P1 = DataDensityMap(m(r(1):r(2),c(1):c(2))) > 0.9;

    % check for redundant scene
    if ~any(P1(:)); fprintf('redundant scene, skipping \n'); continue; end; 
    
    
    % coregister this scene to the strip mosaic, only use areas with > 95%
    [~,p,perr,rmse(i)] = ...
            coregisterdems(Xsub(c(1):c(2)),Ysub(r(1):r(2)),...
            Zsub(r(1):r(2),c(1):c(2)),...
            x(c(1):c(2)),y(r(1):r(2)),...
            z(r(1):r(2),c(1):c(2)),...
            P0,P1);
    trans(:,i) = [p;perr];
    
    %check for segment break
    if isnan(rmse(i)) || (rmse(i) > max_coreg_rmse)
        fprintf('Unable to coregister, breaking segment\n')
        f=f(1:i-1);
        trans = trans(:,1:i-1);
        rmse  = rmse(1:i-1);
        break
    end
    
    % interpolation grid
    xi=x - trans(2,i);
    yi=y - trans(3,i);
    
    %check uniform spacing is maintained (sometimes rounding errors)
    if length(unique(diff(xi))) > 1; xi=round(xi,4); end
    if length(unique(diff(yi))) > 1; yi=round(yi,4); end
    
    % interpolate the floating data to the reference grid
    zi = interp2(xi,yi,z-trans(1,i),Xsub(:)',Ysub(:),'*linear');
    clear z
    
    % interpolate the mask to the same grid
    mi = interp2(xi,yi,single(m),Xsub(:)',Ysub(:),'*nearest');
    mi(isnan(mi)) = 0; % convert back to uint8
    mi = logical(mi);
    clear m
    
    % interpolate ortho to same grid
    oi = single(o);
    oi(oi==0) = NaN; % set border to NaN so wont be interpolated
    oi = interp2(xi,yi,oi,Xsub(:)',Ysub(:),'*cubic');
    clear o 
    
    clear Xsub Ysub
    
    % remove border 0's introduced by nn interpolation
    M3 = ~isnan(zi);
    M3=imerode(M3,ones( 6 )); % border cutline
    
    zi(~M3)=NaN;
    mi(~M3)=0;
    clear M3 
    
    % remove border on orthos seperately
    M4= ~isnan(oi);
    M4=imerode(M4,ones( 6 ));
    oi(~M4)=NaN;
    clear M4
    
    % make weighted elevation grid
    A=  Zsub.*W + zi.*(1-W);
    A( isnan(Zsub) & ~isnan(zi))=  zi( isnan(Zsub) & ~isnan(zi));
    A(~isnan(Zsub) &  isnan(zi))=Zsub(~isnan(Zsub) &  isnan(zi));
    clear zi Zsub
    
    % put strip subset back into full array
    Z(r0:r1,c0:c1) = A;
    clear A
    
    % for the matchtag, just straight combination
    M(r0:r1,c0:c1) = Msub | mi;
    clear Msub Mi
  
     % make weighted ortho grid
     Osub = single(Osub);
     Osub(Osub==0) = NaN;
     A=  Osub.*W + oi.*(1-W);
    
    clear W
    
    A( isnan(Osub) & ~isnan(oi))=   oi( isnan(Osub) & ~isnan(oi));
    A(~isnan(Osub) &  isnan(oi))= Osub(~isnan(Osub) &  isnan(oi));
    clear Osub oi
    
    A(isnan(A)) = 0; % convert back to uint16
    A = uint16(A);
    
    O(r0:r1,c0:c1) = A;
    clear A
   
    
end

%crop to data
if exist('Z','var') && any(~isnan(Z(:)))
    [Z,rcrop,ccrop] = cropnans(Z);
    if ~isempty(rcrop)
        X = X(ccrop(1):ccrop(2));
        Y = Y(rcrop(1):rcrop(2));
        M = M(rcrop(1):rcrop(2),ccrop(1):ccrop(2));
        O = O(rcrop(1):rcrop(2),ccrop(1):ccrop(2));
        Z(isnan(Z)) = -9999;
    end
    
else
    X=[]; Y=[]; Z=[]; M=[]; O=[];
end

function [x,y,z,o,m,md] = loaddata(demFile,matchFile,orthoFile,maskFile)
% loaddata load data files and perform basic conversions

d=readGeotiff(demFile); x=d.x; y=d.y; z=d.z; clear d
sz=size(z);

d=readGeotiff(matchFile);
szd=size(d.z);
if any(szd ~= sz)
    m = interp2(d.x,d.y,single(d.z),...
        x(:)',y(:),'*nearest');
    m(isnan(m)) = 0; % convert back to uint16
    m = logical(m);
else
    m=d.z;
end
clear d

o=zeros(size(z));
if exist(orthoFile,'file')
    d=readGeotiff(orthoFile);
    szd=size(d.z);
    if any(szd ~= sz)
        d.z = single(d.z);
        d.z(d.z==0) = NaN; % set border to NaN so wont be interpolated
        o = interp2(d.x,d.y,d.z,...
            x(:)',y(:),'*cubic');
        o(isnan(o)) = 0; % convert back to uint16
        o = uint16(o);
    else
        o=d.z;
    end
    clear d
end

if ~isempty(maskFile)
    if exist(maskFile,'file')
        d=readGeotiff(maskFile);
        szd=size(d.z);
        if any(szd ~= sz); error('maskFile wrong dimensions'); end
        md=d.z;
    else
        error('No maskFile Found')
    end
else
    md = true(sz);
end

z(z < -100 | z == 0 | z == -NaN ) = NaN;

function [x,y,z,o,m] =  applyMasks(x,y,z,o,m,md)

z(~md) = NaN;
m(~md) = 0;
%o(~me) = 0; Note no longer mask the ortho


if any(~isnan(z(:)))
    [z,rcrop,ccrop] = cropnans(z);
    x = x(ccrop(1):ccrop(2));
    y = y(rcrop(1):rcrop(2));
    m = m(rcrop(1):rcrop(2),ccrop(1):ccrop(2));
    o = o(rcrop(1):rcrop(2),ccrop(1):ccrop(2));
end

function [X,Y,Z,M,O] = regrid(X,Y,Z,M,O)

dx = X(2)-X(1);
dy = Y(2)-Y(1);

Xi = X(1) + (dx-rem(X(1)/dx,1)*dx):dx:X(end);
Yi = Y(1) - (rem(Y(1)/dy,1)*dy):dy:Y(end);

Zi = interp2(X,Y(:),Z,Xi,Yi(:),'*linear');

M = interp2(X,Y(:),single(M),Xi,Yi(:),'*nearest');
M(isnan(M)) = 0; % convert back to uint8
M = logical(M);

% interpolate ortho to same grid
O = single(O);
O(isnan(Z)) = NaN; % set border to NaN so wont be interpolated
O = interp2(X,Y(:),O,Xi,Yi(:),'*cubic');
O(isnan(O)) = 0; % convert back to uint16
O = uint16(O);

Z = Zi;
X = Xi;
Y = Yi;


function [A,r,c] = cropnans(varargin)
% cropnans crop array of bordering nans
%
% [A,r,c] = cropnans(A)

A=varargin{1};
buff=0;
if nargin == 2; buff=varargin{2}; end

r = [];
c = [];


M = ~isnan(A);

if ~any(M(:)); return; end

rowsum = sum(M) ~= 0;
colsum = sum(M,2) ~= 0;

c(1) = find(rowsum,1,'first')-buff;
c(2) = find(rowsum,1,'last')+buff;

r(1) = find(colsum,1,'first')-buff;
r(2) = find(colsum,1,'last')+buff;

if c(1) < 1; c(1)=1; end
if r(1) < 1; r(1)=1; end
if c(2) > size(A,2); c(2)=size(A,2); end
if r(2) > size(A,1); r(2)=size(A,1); end

A = A(r(1):r(2),c(1):c(2));
