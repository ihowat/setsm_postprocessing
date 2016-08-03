function [p,dz] = registerDEM2LIDAR(x,y,z,rx,ry,rz)
% registerDEM2LIDAR Register raster DEM to LIDAR point cloud
%
% [dtrans,dz] = registerDEM2LIDAR(x,y,z,rx,ry,rz) where x,y,z are the DEM
% coordinate vectors and array and rx,ry,rz are the control point vectors.
% The output p is the dz,dx,dy offets and dz are the residuals, so that the
% corrected dem is x+p(2),y+(3),z+p(1) with the 1-sigma of the fit as 
% nanstd(dz).



if iscell(z)
    
    x = cellfun(@(x) double(x), x, 'UniformOutput',false);
    y = cellfun(@(x) double(x), y, 'UniformOutput',false);
    z = cellfun(@(x) double(x), z, 'UniformOutput',false);
    
    % slopes
    dx = diff(x{1}(1:2));
    [sx,sy] = cellfun(@(x) gradient(x,dx), z, 'UniformOutput',false);
    sx = cellfun(@(x) -x, sx, 'UniformOutput',false);
    
else

    x = double(x);
    y = double(y);
    z = double(z);
    
    
    % slopes
    [sx,sy] = gradient(z,x(2)-x(1));
    sx = -sx;
    
end


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
dz = [];
dtrans = [];
p = [0;0;0];

%% Planimetric Correction Iteration Loop
while it
    
    % interpolate DEM to lidar x,y
    
    if iscell(z)
        
        zi = nan(size(rz)); sxi = zi; syi = zi;
        i=1;
        for i=1:length(z);            
            zi(i)   = interp2(xn{i},yn{i},zn{i},rx(i),ry(i));
            sxi(i)  = interp2(xn{i},yn{i},sx{i},rx(i),ry(i));
            syi(i)  = interp2(xn{i},yn{i},sy{i},rx(i),ry(i));
        end
    else
    
        zi  = interp2(xn,yn,zn,rx,ry);
        sxi  = interp2(xn,yn,sx,rx,ry);
        syi  = interp2(xn,yn,sy,rx,ry);
    end
    
    % difference DEM and LIDAR
    dzn = zi-rz;
    
    % filter NANs for speed
    n = ~isnan(sxi) & ~isnan(syi) & ~isnan(dzn);
    % filter outliers for solution
    n = n & (abs(dzn - median(dzn(n))) < std(dzn(n)));
    
    % if too few non-nan points found, return
     if sum(n) < 4
        fprintf('Too few (%d) non-nan DEM points, returning\n',sum(n));
        return;
    end

    % calculate & display rmse
    d1 = sqrt(mean(dzn(n).^2));
    fprintf('N=%d  dz=%.2f dx=%.2f dy=%.2f rmse=%.3f\n ',...
        sum(n),sum([dtrans,p],2)',d1)
    
    % break if rmse increases, keeping last iteration values
    if d1 >= d0
        fprintf('returning values from iteration n-1\n')
        break
    end
    
    % if better, then set new values to old
    z   = zn;
    x   = xn;
    y   = yn;
    dz  = dzn;
    
    dtrans = [dtrans,p];
    
    % break if rmse decrease is small, keeping these iteration values
    if (d0 - d1)./d0 <= 1e-3
        break
    end
    
    % dependent variable matrix
    X = [ones(size(dz(n))),sxi(n),syi(n)];
    
    % solve for parameters
    %warning off
    p = X\dz(n);
    % warning on
    
    % apply offsets to generate new values to test
    if iscell(z)
        zn = cell(size(z)); xn = zn; yn = zn;
        i=1;
        for i=1:length(z);
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
    it = it+1;
    d0 = d1;
    
end

p = -sum(dtrans,2);

% 
% %% x,y,z - Dependent Bias Correction
% %disp('Ramp and Elevation Bias Correction')
% 
% % filter NaNs and outliers for solution
% n = ~isnan(sxi) & ~isnan(syi) & ~isnan(dz);
% n = n & (abs(dz - median(dz(n))) < 2.*std(dz(n)));
% 
% % dependent variable matrix
% X = [ones(size(dz(n))),rz(n),rx(n),ry(n)];
% 
% % solve for parameters and display offsets
% %p = X\dz(n);
% 
% [b,bint,r,rint,stats] = regress(dz(n),X);
% 
% 
% 
% % grid dem coordinates and apply corrections
% [xgrid,ygrid] = meshgrid(x,y);
% z = z - (p(1) + p(2).*z + p(3).*xgrid + p(4).*ygrid);
% 
% delev = p;
% 
% % calculate final rmse
% zi  = interp2(x,y,z,rx,ry);
% dz = zi-rz;
% n = ~isnan(dz);
% % filter outliers for solution
% n = n & (abs(dz - median(dz(n))) < 2.*std(dz(n)));
% zrms = rms(dz(n));
% zmed = median(dz(n));
% zmn = mean(dz(n));
% fprintf('rmse = %.2f, dz mean = %.2f, dz median = %.2f \n',...
%     zrms,zmed,zmn);
