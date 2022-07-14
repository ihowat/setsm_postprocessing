function z=addSeaSurfaceHeight(x,y,z,land,varargin)
% Add sea surface height above the ellipsoid where land mask is false
%
%  z=addSeaSurfaceHeight(x,y,z,land)
% 
%  z=addSeaSurfaceHeight(...,'epsg',value) speficifies the projection of x
%  and y by the epsg value (3413 is default);
%
%  z=addSeaSurfaceHeight(...,'adaptCoastline')

% parse vargins
epsg = 3413;
n = find(strcmpi(varargin,'epsg'));
if ~isempty(n)
    epsg = varargin{n+1};
end

adaptCoastlineFlag = false;
n = strcmpi(varargin,'adaptCoastline');
if any(n)
    adaptCoastlineFlag = true;
end

if ~any(~land(:))
    fprintf('land mask has no zero (water) values\n')
    return
end

% grid of x and y coords
[X,Y] = meshgrid(single(x),single(y));

% convert to lat lon coordinates
if epsg == 3412
    [LAT,LON]=polarstereo_inv(X,Y,[],[],70,-45);
elseif epsg == 3031
    [LAT,LON]=polarstereo_inv(X,Y,[],[],-71,0);
end

% egm96 sea level height above 
ellipsoidHeight = egm96geoid(LAT,LON);

% set heights below the sea level height to sea level height
M = z < ellipsoidHeight;
z(M) = ellipsoidHeight(M);


if adaptCoastlineFlag
    % adapt coastline to ice front change/ice bergs
    M = imdilate(land,ones(100));
    M = M & ( (z-ellipsoidHeight) > 2 );
    
    L = bwlabel(M);
    
    uniqueL = unique(L(:));
    uniqueLLand  = unique(L(land(:)));
    
    removeTheseLs = setdiff(uniqueL,uniqueLLand);
    
    for i=1:length(removeTheseLs)
        M(L == removeTheseLs(i)) = false;
    end
    
    land = M;
end

z(~land) = ellipsoidHeight(~land);

fprintf('sea surface height applied\n')