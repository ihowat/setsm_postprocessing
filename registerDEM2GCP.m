function [p,dz,dtrans,sigma_dtrans,p68,p90,N] = registerDEM2GCP(x,y,z,rx,ry,rz,varargin)
% registerDEM2GCP Register raster DEM to LIDAR point cloud
%
% p = registerDEM2GCP(x,y,z,rx,ry,rz) where x,y,z are the DEM
% coordinate vectors and array and rx,ry,rz are the control point vectors.
% The output p is the dz,dx,dy offsets so that to apply the correction is:
%    zn = z - p(1);
%    xn = x - p(2);
%    yn = y - p(3);
%
% The complete outputs are:
% [p,dz,dtrans,sigma_dtrans,fit_stats,p68,p90,N] = registerDEM2GCP(...)
% where:
% dz           =    the difference between the lidar elevation points and
%                   the DEM elevations interpolted to the lidar point cloud
%                   (lidar minus dem)
% dtrans        =   the correction parameters at each iteration
% sigma_dtrans  =   the 1sigma errors of each dtrans
% p68,p90       =   the 68th and 90th percentiles of the absolute vertical
%                   residuals [abs(dz)] for each iteration.
% N             =   index of lidar points used in each iteration (i.e.
%                   points not filtered or not over NaNs on the DEM.
%
% Optional Inputs:
%  [...] = registerDEM2GCP(...,'maxOffset',val)
% where maxOffset is the maximum horizontal offset allowed, above which the
% iterations will terminate and only the vertical and rotation-corrections
% will be returned.
%
% Ian Howat, ihowat@gmail.com
% '09-Dec-2020 10:38:43'

% Default maximum offset allowed
maxp = 15;

% set output variables
dz = []; % vector of vertical differences between DEM and lidar
dtrans = []; % translation vector [dz,dx,dy, dzdx,dzdy]
sigma_dtrans = []; % std dev of translations
p = [0;0;0]; % new iteration translations
p68=[]; % 68th percentile of the absolute value of dz
p90=[];% 90th percentile of the absolute value of dz
N=logical([]); % index of points used in the regression
perr = [0;0;0];

if iscell(z) % cell block mode
    
    % insure input DEM are double
    x = cellfun(@(x) double(x), x, 'UniformOutput',false);
    y = cellfun(@(x) double(x), y, 'UniformOutput',false);
    z = cellfun(@(x) double(x), z, 'UniformOutput',false);
    
    % get slopes iterating through each ceell
    sx=cell(size(x));
    sy=cell(size(y));
    
    i=1;
    for i=1:length(x)
        
        % make sure cell is properly formed
        if length(x{i}) < 3 || length(y{i}) < 3
            fprintf('dem must be at least 3x3, returning\n')
            return
        end
        
        [sx{i},sy{i}]=gradient(z{i},x{i},y{i});
        sx{i}=-sx{i};
        sy{i}=-sy{i};
        
    end
    
else % single array mode
    
    % insure DEM is double
    x = single(x);
    y = single(y);
    z = single(z);
    
    % make array is large enough
    if length(x) < 3 || length(y) < 3
        fprintf('dem must be at least 3x3, returning\n')
        return
    end
    
    % slopes
    [sx,sy] = gradient(z,x,y);
    sx = -sx;
    sy = -sy;
    
end

% parse optional inargs
n=find(strcmpi(varargin,'maxOffset'));
if ~isempty(n)
    maxp=varargin{n+1}(1); % reset max offset
end

% insure inputs are double for interpolation
rx = double(rx);
ry = double(ry);
rz = double(rz);

% set new iteration values
xn = x;
yn = y;
zn = z;

% set iteration loop variables
d0 = inf;
it = 1;

%% Planimetric Correction Iteration Loop
while it
    
    % interpolate DEM to lidar x,y
    if iscell(z) % iterate  on each cell
        
        % create vector same size as lidar data and fill by iteration
        zi = nan(size(rz)); sxi = zi; syi = zi;
        i=1;
        for i=1:length(z)
            zi(i)   = interp2(xn{i},yn{i},zn{i},rx(i),ry(i));
            sxi(i)  = interp2(xn{i},yn{i},sx{i},rx(i),ry(i));
            syi(i)  = interp2(xn{i},yn{i},sy{i},rx(i),ry(i));
        end
    else
        % in array mode, simply interp to lidar x,y
        zi  = interp2(xn,yn,zn,rx,ry,'*linear');
        sxi  = interp2(xn,yn,sx,rx,ry,'*linear');
        syi  = interp2(xn,yn,sy,rx,ry,'*linear');
    end
    
    % difference DEM and LIDAR
    dzn = zi-rz;
    
    % calculate 70th and 90th percentiles of absolute dz residulas
    d1 = prctile(abs(dzn),68);
    d2 = prctile(abs(dzn),90);
    
    ptdzn = prctile(dzn,[5 95]);
    ptsx = prctile(sxi,[5 95]);
    ptsy = prctile(syi,[5 95]);
    
    n = dzn > ptdzn(1) & dzn < ptdzn(2) & sxi > ptsx(1) & sxi < ptsx(2) &...
        syi > ptsy(1) & syi < ptsy(2);

    % if too few non-nan points found, return
    if sum(n) < 4
        fprintf('Too few (%d) non-nan DEM points, returning\n',sum(n));
        return;
    end
    
    % calculate & display rmse
    fprintf('it=%d N=%d dz=%.2f dx=%.2f dy=%.2f median=%.3f LE68=%.3f LE90=%.3f\n',...
        it-1,sum(n),sum([dtrans,p],2)',nanmedian(dzn),d1,d2)
   
    % break if LE68 increases or remains constant, keeping last iteration values
    if d1 >= d0
        fprintf('increased error this iteration, stopping\n');
        break
    end

    
    % if better, then set new values to old
    z   = zn;
    x   = xn;
    y   = yn;
    dz  = dzn;
    
    % add shifts and stats to record arrays
    dtrans = [dtrans,p]; % add offset (5 x # of iterations)
    sigma_dtrans = [sigma_dtrans,perr]; % 1 sigma +/- of shifts
    p68 = [p68,d1]; % 70th percentile of abs value of residuals
    p90 = [p90,d2];% 90th percentile of abs value of residuals
    N = [N,n]; % index of points used
    
    % check for convergence
    if (d0 - d1) <= 1e-3
            fprintf('improvement less than threshold, stopping\n');
        break
    end
    
    d0 = d1;
    
    X = [ones(size(dz(n))),sxi(n),syi(n)];% solve for x,y,z shifts
    
    % solve for parameters
    warning off % turn off rank deficiency messages, will instead asses stats
    p = X\dz(n); % old solver, fast but does not provide errors ests
    %[p,pint,~,~,stats]=regress(dz(n),X,0.32); %soltn with 1sigma errors
    warning on
    
    % calculate p errors
    [~,R,perm] = qr(X,0);
    RI = R\eye(3);
    nu = size(X,1)-size(X,2); % Residual degrees of freedom
    yhat = X*p;                     % Predicted responses at each data point.
    r = dz(n)-yhat;                     % Residuals.
    normr = norm(r);
    
    rmse = normr/sqrt(nu);      % Root mean square error.
    tval = tinv((1-0.32/2),nu);
    
    se = zeros(size(X,2),1);
    se(perm,:) = rmse*sqrt(sum(abs(RI).^2,2));
    perr = tval*se;

    %  check if total horizontal shift exceeds maxp
    if any(abs(sum([dtrans(2:3,:),p(2:3)],2)) > maxp)
        fprintf('maximum horizontal offset reached\n')
        break
    end
    
    % apply offsets to generate new values to test
    if iscell(z) % apply iteratively to cells
        zn = cell(size(z)); xn = zn; yn = zn;
        i=1;
        for i=1:length(z)
            zn{i} = z{i} - p(1);
            xn{i} = x{i} - p(2);
            yn{i} = y{i} - p(3);
        end
    else
        zn = z - p(1);
        xn = x - p(2);
        yn = y - p(3);
    end
    % next interation
    it=it+1;
end

p = sum(dtrans,2);


%       % z.zn = interp2(z.x - p(2), z.y - p(3),z.z - p(1),z.x,z.y,'*linear');
% X = atm.dz(abs(atm.dz) < 10);
% X = atm.dzi(abs(atm.dzi) < 10);
% [N,edges] = histcounts(X, 'Normalization','pdf');
% edges = edges(2:end) - (edges(2)-edges(1))/2;
% plot(edges, N);
