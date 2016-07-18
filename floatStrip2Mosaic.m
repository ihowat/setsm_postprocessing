function N = floatStrip2Mosaic(metaFile,m,dy0,N,varargin)

% read this strip data
[x,y,z,mt,or,c,r] = readStripInTile(metaFile,m.x,m.y);

if isempty(x);
    m.dtrans=[m.dtrans,[NaN;NaN;NaN]];
    m.rmse=[m.rmse,NaN];
    return;
end


c=find(c);
r=find(r);

if length(c)*length(r) < 1000;
    fprintf('%d pixel overlap is too small, skipping\n',sum(c)*sum(r));
    m.dtrans=[m.dtrans,[NaN;NaN;NaN]];
    m.rmse=[m.rmse,NaN];
    return;
end

% check for new coverage
Nsub=N(r,c);

if sum(Nsub(:)==0 & ~isnan(z(:))) < 1000;
    disp('redundant coverage, skipping');
    m.dtrans=[m.dtrans,[NaN;NaN;NaN]];
    m.rmse=[m.rmse,NaN];
    return;
end

% crop to new coverage
A = Nsub ~= 0 & ~isnan(z); % overlapping coverage mask
co= find(sum(A) ~= 0,1,'first'):find(sum(A) ~= 0,1,'last');
ro= find(sum(A,2) ~= 0,1,'first'):find(sum(A,2) ~= 0,1,'last');

if isempty(co) || isempty(ro);	
	fprintf('No overlapping coverage\n');
	return;
end


% apply registration if given
if length(varargin) == 1;
   
    m.dtrans=[m.dtrans,varargin{1}'];
    m.rmse=[m.rmse,0];
    
else
    
    % co-register strip to mosaic subset
    [~,m.dtrans(:,size(m.dtrans,2)+1), m.rmse(1,length(m.rmse)+1) ] = ...
        coregisterdems(m.x(1,c(co)),m.y(r(ro),1),m.z(r(ro),c(co)),...
        x(co),y(ro),z(ro,co),m.mt(r(ro),c(co)),mt(ro,co));    
    
    % check for coregistration failure
    if isnan(m.rmse(1,end)); return; end;
    
end

% Apply transform to x,y,z
x=x - m.dtrans(2,end);
y=y - m.dtrans(3,end);
z=z - m.dtrans(1,end);

% interpolate transformed surface to tile grid
[z,mt,or] = interpolate2grid(x,y,z,mt,or,m.x(1,c),m.y(r,1));

% make date grid
dy=~isnan(z).*dy0;

if any(A(:)); % check if any overlapping
    
    %%remove ovelapping area with existing data with buffer
    dx = m.x(1,2)-m.x(1,1);
    buff=10*dx+1;
    %A = imerode(A,ones(buff));
    %z(A)=NaN;
    %mt(A)=0;
    %or(A)=0;
    %clear A;
    
    % New Feathering Code
    W = edgeFeather(Nsub~=0,~isnan(z),buff);
    
    % make weighted elevation grid
    zsub=m.z(r,c);
    A=  zsub.*W + z.*(1-W);
    n = isnan(zsub) & ~isnan(z);  A(n)= z(n);
    n = ~isnan(zsub) & isnan(z);  A(n)= zsub(n);
    
    % put strip subset back into full array
    m.z(r,c) = A;
    clear A 
    
    Nsub = Nsub + uint8(z.*(1-W) > 0);
    N(r,c)=Nsub;
    clear Nsub z
    
    % for the matchtag, just straight combination
    m.mt(r,c) = m.mt(r,c) | mt;
    
    clear mt
    
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
else
    
    
    
    N(r,c)= uint8(~isnan(z));
    m.z(r,c) = z; clear z
    m.mt(r,c) = mt; clear mt
    m.or(r,c) = or; clear or
    m.dy(r,c) = uint16(dy); clear dy
end


% %% Match density based weighting
% % data density masks
% Pm = DataDensityMap(mtsub,21);
% Pi = DataDensityMap(mi,21);
%
% % set density to zero where no data exists
% Pm(isnan(zsub)) =0;
% Pi(isnan(zi)) =0;
%
% % sum of densities
% sumP = Pm + Pi;
%
% % put density-scaled mosaic and strip subs into merge array
% ztemp(:,:,1) = zsub.*Pm./sumP;
% ztemp(:,:,2) = zi.*Pi./sumP;
% a = sum(isnan(ztemp),3) == 2; % keep track of missing data
% ztemp = nansum(ztemp,3); % sum scaled data to get weighted avg
% ztemp(a) = NaN; % put back NaNs
% z(r,c) = ztemp;
% clear ztemp zsub zi a
%
% % put density-scaled date stamps
% dytemp = (double(dysub).*Pm./sumP) + (double(dyi).*Pi./sumP);
% dytemp = int16(dytemp);
% dy(r,c) = dytemp;
% clear dytemp dysub dyi
%
% ortemp = (double(orsub).*Pm./sumP) + (double(oi).*Pi./sumP);
% ortemp =   int16(ortemp);
% or(r,c) = ortemp;
% clear ortemp orsub oi
%
% clear Pm Pi sumP
%
% mt(r,c) = mtsub | mi;
%
% clear mtsub mi

