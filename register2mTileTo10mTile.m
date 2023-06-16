function register2mTileTo10mTile(tile2m,tile10m)
% register 2m quad tile to 10m full tile w/ translation and dz warp
%
% register2to10(2mTileName,10mTileName)

% set matfile handles
m2=matfile(tile2m);
m10=matfile(tile10m);

% set outname
outName=strrep(tile2m,'.mat','_reg.mat');

% read 2m x,y,z
x=m2.x;
y=m2.y;
z=m2.z;

% read 10m, loading only what overlaps 2m quad tile
cols = find(m10.x >= min(m2.x) & m10.x <= max(m2.x));
rows = find(m10.y >= min(m2.y) & m10.y <= max(m2.y));

xr = m10.x(1,cols);
yr = m10.y(1,rows);
zr = m10.z(rows,cols);

% downsample 2m to 10m grid
z_interp = interp2(x,y(:),z,xr,yr(:),'*linear');

% coregister 2m downsampled to 10m reference
[~,reg.p,reg.perr,reg.sd] = coregisterdems(xr,yr,zr,xr,yr,z_interp);

clear z_interp

% apply registration (transform) to 2m quad tile
z_reg = interp2(x - reg.p(2),y(:) - reg.p(3),z - reg.p(1),...
    x,y(:),'*linear');
    
clear z

%smooth the 2m grid to 10m through 5x5 kernel mean
z_reg_smooth = conv2(z_reg,(1/25)*ones(5),'valid');

% interpolate the smoothed 2m grid to the 10m grid
z_reg_smooth_interp = interp2(x(3:end-2),y(3:end-2)',z_reg_smooth,xr,yr(:),'*linear');

clear z_reg_smooth

% difference smoothed & downsampled 2m and 10m reference
dz = z_reg_smooth_interp - zr;

clear z_reg_smooth_interp 

%replace border NaNs with nearest values
N = ~isnan(dz);
% left
n = find(any(N),1,'first');
dz(:,1:n-1) = repmat(dz(:,n),1,n-1);
% right
n = find(any(N),1,'last');
n = size(dz,2) - n;
if n > 0
    dz(:,end-n+1:end) = repmat(dz(:,end-n),1,n);
end
% top
n = find(any(N,2),1,'first');
dz(1:n-1,:) = repmat(dz(n,:),n-1,1);
% bottom
n = find(any(N,2),1,'last');
n = size(dz,1) - n;
if n > 0
    dz(end-n+1:end,:) = repmat(dz(end-n,:),n,1);
end

% smooth the dz map to preserve small scale variability due to resolution
k=51; %square kernel dimension
dz_smooth = conv2(dz,(1/(k^2))*ones(k),'valid');

dz_smooth = padarray(dz_smooth,(size(dz)-size(dz_smooth))./2,'replicate');
    
clear dz

% fill remaining NaNs
dz_smooth = inpaint_nans(double(dz_smooth),2);
dz_smooth = single(dz_smooth);

% upsample difference map 2m grid
dz_smooth_interp = interp2(xr,yr(:),dz_smooth,x,y(:),'*linear');
   
clear dz_smooth

% appy upsampled to map to transformed 2m DEM and add to file
z = z_reg - ;
  
% save to out matfile
save(outName,'reg','x','y','z','dz_smooth_interp','-v7.3');

m2reg=matfile(outName);
m2reg.Properties.Writable = true;

% make coregistration cluster mask - not currently used
C = true(nrows,ncols);

% register z_mad and add to mat file
m2reg.z_mad =  applyRegistration2Var(reg.p,m2,C,'gridvar','z_mad','subsetSize',5000);

% register N and add to mat file
m2reg.N =  applyRegistration2Var(reg.p,m2,C,'gridvar','N','subsetSize',5000);

% register Nmt and add to mat file
m2reg.Nmt =  applyRegistration2Var(reg.p,m2,C,'gridvar','Nmt','subsetSize',5000);

% register tmin and add to mat file
m2reg.tmin =  applyRegistration2Var(reg.p,m2,C,'gridvar','tmin','subsetSize',5000);

% register tmax and add to mat file
m2reg.tmax = applyRegistration2Var(reg.p,m2,C,'gridvar','tmax','subsetSize',5000);

% register land mask to mat file if exists
if any(strcmp(m_fields,'land'))
    m2reg.land = applyRegistration2Var(reg.p,m2,C,'gridvar','land','subsetSize',5000);
end

% additional fields to transfer over
add_fields = {'stripList', 'version'};
for i=1:length(add_fields)
    fld = add_fields{i};
    if any(strcmp(m_fields,fld))
        eval(['m2reg.',fld,' = m2.',fld,';']);
    end
end


function z = applyRegistration2Var(dtrans,m,N,varargin)
% ApplyRegistration2Var applies 3-D transformation to grid. V4 compatible
%   z = applyRegistration2Var(dtrans,m,z,N) dtrans is [x,y,z] offsets to be
%   *subtracted* . m is matfile object. N is coregistration mask.
%   z = applyRegistration2Var(...,'gridVar','varname')
%   z = applyRegistration2Var(...,'subsetSize',length)

%check dtrans
if numel(dtrans) ~= 3
    error('dtrans must be a 3 element vector')
end
if any(isnan(dtrans))
    error('dtrans contains NaNs')
end

% parse/error check gridVar argument
gridVar='z'; %default
subsetSize=5000; % dimension of interpoaltion subset(pixels=subsetSize^2)


if ~isempty(varargin)
  
    n = find(strcmpi(varargin,'gridvar'));
    if ~isempty(n)
        gridVar=varargin{n+1};
    end
    
    n = find(strcmpi(varargin,'subsetsize'));
    if ~isempty(n)
        subsetSize=varargin{n+1};
    end
end

% extract grid vectors for speed rather in pulling on each iter
x=m.x;
y=m.y;
y = y(:);
res = x(2)-x(1);

% to speed up/converve memory, we'll perform gridded interpolation over
% 1000x1000 pixel subsets, overlapping by one pixel. 
% Make a list of min max indices.
sz = [length(y) length(x)];

if subsetSize < 0.75*sz(2)
    
    intsCol=1:subsetSize:sz(2);
    intsCol=[intsCol(1:end-1)',intsCol(2:end)'];
    intsCol(end)=sz(2);
    intsCol(1:end-1,2)=intsCol(1:end-1,2)-1;
    
else
    intsCol =[1 sz(2)];
end


if subsetSize < .75*sz(1)
    
    intsRow=1:subsetSize:sz(1);
    intsRow=[intsRow(1:end-1)',intsRow(2:end)'];
    intsRow(end)=sz(1);
    intsRow(1:end-1,2)=intsRow(1:end-1,2)-1;
    
else
    intsRow =[1 sz(1)];
end

% need to add pixel buffer to int to account for shifts, scaling with dtrans
buff=max(ceil(abs(dtrans(2:3)./res)));
if size(intsCol,1) == 1
    buffsCol = [0,0];
else
    buffsCol=[[0 buff];repmat([-buff buff],size(intsCol,1)-2,1);[-buff 0]];
end

if size(intsRow,1) == 1
    buffsRow=[0,0];
else
    buffsRow=[[0 buff];repmat([-buff buff],size(intsRow,1)-2,1);[-buff 0]];
end

%% nested interpolation loops with interpolation dependent on grid variable

% initialize output
switch gridVar
    case {'z','z_mad'}
        z = nan(sz,'single');
    case {'tmin','tmax'}
        z= zeros(sz,'uint16');
    case {'N','Nmt','land'}
        z = zeros(sz,'uint8');  
end

% subset interpolation loop
for i=1:size(intsRow,1)
    
    % make subset row vector with buffer
    intRowBuff=intsRow(i,1)+buffsRow(i,1):intsRow(i,2)+buffsRow(i,2);
    
    for j=1:size(intsCol,1)
        
        % make subset column vector with buffer
        intColBuff=intsCol(j,1)+buffsCol(j,1):intsCol(j,2)+buffsCol(j,2);
        
        % extract subset from unregistered tile with buffer and add offset
        eval(['z0 = m.',gridVar,'(intRowBuff,intColBuff);']);
        
        % convert to single for interpolation
        z0 = single(z0);
        
        % apply coregistration cluster mask
        z0(~N(intRowBuff,intColBuff)) = NaN;
        
        % check for empty subset to skip
        if ~any(~isnan(z0(:)) & z0(:) ~= 0); continue; end
        
        % make unbuffered row and column vectors to insert into new
        % array
        intRow=intsRow(i,1):intsRow(i,2);
        intCol=intsCol(j,1):intsCol(j,2);
        
        x1 = x(intColBuff)-dtrans(2);
        y1 = y(intRowBuff)-dtrans(3);
        
        % interpolate depending on variable
        switch gridVar
            
            case 'z'
                
                % apply xy,z offsets and interpolate
                z(intRow,intCol) =interp2(x1,y1, z0-dtrans(1),...
                    x(intCol),y(intRow),'*linear',NaN);

            case 'z_mad' % error grid
                
                % apply only x and y offsets and interpolate
                z(intRow,intCol) =interp2(x1,y1, z0,...
                    x(intCol),y(intRow),'*linear',NaN);
                
            case {'N','Nmt','land'} % stack, matchtag numbers and land mask
                
                z0(isnan(z0))=0;
                
                % apply horizontal offsets and interpolate
                z0=interp2(x1,y1, z0, x(intCol),y(intRow),'*nearest',0);
                
                % convert back to uint8 and insert into array
                z0(isnan(z0)) = 0; 
                z(intRow,intCol)=uint8(z0);
                
                clear z0
                
            case 'or'
                
                % apply horizontal offsets and interpolate
                z0=interp2(x1,y1, z0, x(intCol),y(intRow),'*cubic',0);
                
                % convert back to int16 and insert into array
                z0(isnan(z0)) = 0; 
                z(intRow,intCol)=int16(z0);
                
            case {'tmin','tmax'}
                
                % apply horizontal offsets and interpolate
                z0=interp2(x1,y1, z0, x(intCol),y(intRow),'*linear',0);
                
                z0(isnan(z0)) = 0; % convert back to uint8
                z(intRow,intCol)=uint16(z0);
                
            otherwise
                error('grid variable input argument not recognized, must be ''z'',''z_mad'',''tmin'',''tmax'',''N'',''land'' or ''Nmt''')
                  
        end
        
        clear z0
        
    end
end
