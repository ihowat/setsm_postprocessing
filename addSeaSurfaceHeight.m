function [z,seaSurfaceMask]=addSeaSurfaceHeight(x,y,z,land,varargin)
% Add sea surface height above the ellipsoid where land mask is false
%
%  z=addSeaSurfaceHeight(x,y,z,land)
% 
%  z=addSeaSurfaceHeight(...,'epsg',value) speficifies the projection of x
%  and y by the epsg value (3413 is default);
%
%  z=addSeaSurfaceHeight(...,'adaptCoastline')

%epsg = 3413;
epsg = [];

land = logical(land);

seaSurfaceMask = false(size(z));

% parse vargins
n = find(strcmpi(varargin,'epsg'));
if ~isempty(n)
    epsg = varargin{n+1};
end
if isempty(epsg)
    error('epsg varargin must be provided')
end

adaptCoastlineFlag = false;
n = strcmpi(varargin,'adaptCoastline');
if any(n)
    adaptCoastlineFlag = true;
end

landIsQcMaskFlag = false;
n = strcmpi(varargin,'landIsQcMask');
if any(n)
    landIsQcMaskFlag = true;
end

n = find(strcmpi(varargin,'doNotFillBorder'));
if ~isempty(n)
    doNotFillBorder = varargin{n+1};
else
    doNotFillBorder = [];
end

if ~any(~land(:))
    fprintf('land mask has no zero (water) values\n')
    return
end

if ~landIsQcMaskFlag
    landfraction = nnz(land) / numel(land);
    fprintf('tile land fraction: %g\n', landfraction)
    if landfraction > 0.9998
        % This may no longer be needed since the addition of the 'doNotFillBorder' arg,
        % whose purpose is to avoid filling the gaps on tile borders that come from
        % horizontal shifts during registration.
        fprintf('skipping sea surface application for tile that is more than 99.98 percent land\n')
        return
    end
end

% Trying to sample sea level at full (2m) resolution may never return!
% Sample at 100m resolution. EGM96 geoid native resolution is ~15 min anyways.
xi_min = idivide(int64(x(1,1)), int64(100), 'floor') * 100;
xi_max = idivide(int64(x(1,end)), int64(100), 'ceil') * 100;
yi_max = idivide(int64(y(1,1)), int64(100), 'ceil') * 100;
yi_min = idivide(int64(y(1,end)), int64(100), 'floor') * 100;
xi = double(xi_min):100:double(xi_max);
yi = double(yi_max):-100:double(yi_min);

% grid of x and y coords
[X,Y] = meshgrid(single(xi),single(yi));

% convert to lat lon coordinates
if epsg == 3413
    [LAT,LON]=polarstereo_inv(X,Y,[],[],70,-45);
elseif epsg == 3031
    [LAT,LON]=polarstereo_inv(X,Y,[],[],-71,0);
else
    error('epsg not handled: %d',epsg)
end

% egm96 sea level height above
fprintf('retrieving ellipsoid heights\n')
try
    % Note that older matlab versions egm96geoid don't support
    % egm96geoid(LAT,LON) syntax and will error.
    ellipsoidHeight = egm96geoid(LAT,LON);
catch ME
    warning('Caught err in egm96geoid method')
    fprintf(1,'Err identifier: %s\n',ME.identifier);
    fprintf(1,'Err message:\n%s\n',ME.message);
    fprintf(1,'Trying geoidheight method instead\n')

%    ellipsoidHeight = geoidheight(LAT,LON,'egm96');
%    ellipsoidHeight = reshape(ellipsoidHeight, size(LAT));

    % The above 2-line approach is bugged when LON values
    % cross the 180 line (both +179.X and -179.X values exist).
    % Call geoidheight method on those two groups of LON value
    % points separately.
    LON_pos_ind = find(LON >= 0);
    LON_neg_ind = find(LON < 0);
    ellipsoidHeight_pos = geoidheight(LAT(LON_pos_ind),LON(LON_pos_ind),'egm96');
    ellipsoidHeight_neg = geoidheight(LAT(LON_neg_ind),LON(LON_neg_ind),'egm96');
    ellipsoidHeight = NaN(size(LAT),'single');
    ellipsoidHeight(LON_pos_ind) = ellipsoidHeight_pos;
    ellipsoidHeight(LON_neg_ind) = ellipsoidHeight_neg;
    if any(isnan(ellipsoidHeight))
        error('failed to assemble ellipsoid heights')
    end
end
fprintf('got ellipsoid heights\n')

ellipsoidHeight = interp2(xi,yi(:),ellipsoidHeight,x,y(:),'*bilinear');

% set heights below the sea level height to sea level height
M = z < ellipsoidHeight;
if ~isempty(doNotFillBorder)
    overlap = (M & doNotFillBorder);
    if any(overlap(:))
        fprintf('preventing waterfill around border of DEM that is extrapolated NoData from horizontal shifts in registration\n')
        M(doNotFillBorder) = 0;
    end
end
z(M) = ellipsoidHeight(M);
seaSurfaceMask = M;


if adaptCoastlineFlag
    % adapt coastline to ice front change/ice bergs
    fprintf('adapting coastline\n')

    M = imdilate(land,ones(100));
    M = M & ( (z-ellipsoidHeight) > 2 );
    
    L = bwlabel(M);
    
    uniqueL = unique(L(:));
    uniqueLLand  = unique(L(find(land)));
    
    removeTheseLs = setdiff(uniqueL,uniqueLLand);
    
    for i=1:length(removeTheseLs)
        M(L == removeTheseLs(i)) = false;
    end
    
    land = M;
end

fill_ocean = ~land;
if ~isempty(doNotFillBorder)
    overlap = (fill_ocean & doNotFillBorder);
    if any(overlap(:))
        fprintf('preventing waterfill around border of DEM that is extrapolated NoData from horizontal shifts in registration\n')
        fill_ocean(doNotFillBorder) = 0;
    end
end

z(fill_ocean) = ellipsoidHeight(fill_ocean);
seaSurfaceMask(fill_ocean) = 1;

fprintf('sea surface height applied\n')
