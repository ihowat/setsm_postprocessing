function registerTile(m,gcp)
% registerTile register tiles to ground control
%
%   registerTile(m,gcp) registers the DEM tile mosaic z in matlab file
%   handle m (with m.z as elevation) to the gcp points in structure gcp,
%   which has fields gcp.x, gcp.y and gcp.z.
% 
%   subfuncs: polyCropGCPs, pointAreaSubsets, registerDEM2LIDAR, 
%   applyRegistration

%% Check for multiple coregistration clusters

% fastest is to look for a zero dtrans's
coregClusters=sum(sum(m.dtrans==0)==3);

% if no un-co-registered data, return
if coregClusters == 0; 
    fprintf('All data is registered, skipping\n');
    return
end

% load coregistration cluster
C=m.C;

sz = whos(m,'z'); sz = sz.size; % image dimensions info
z = nan(sz,'single'); % initialize output

%create output file
outname=strrep(f,'.mat','_reg.mat');

% cluster coregistraton loop
i=1;
for i=1:coregClusters
    
    % make a mask of this cluster
    N = C == i+1;
    
    %% load overlapping GCPs full tile boundaries
    
    % make tile polygon footprint vertices from rectangular coordinate range
    x=m.x;
    y=m.y;
    
    minx=min(x);
    maxx=max(x);
    miny=min(y);
    maxy=max(y);
    
    xv=[minx,minx,maxx,maxx,minx];
    yv=[miny,maxy,maxy,miny,miny];
    
    n=polyCropGCPs(gcp,xv,yv,'rectangle');
    
    % skip if too few
    if length(n) < 4;
        fprintf('%d overlapping points, too few,skipping\n',sum(nn));
        return;
    end
    
    % find gcp's over these pixels
    col = round((gcp.x(n) - minx)/(x(2)-x(1)));
    row = round((gcp.y(n) - maxy)/(y(2)-y(1)));
    col(col < 1)=1; row(row < 1) = 1;
    col(col > length(x)) = length(x); row(row > length(y)) = length(y);
    gcpind = sub2ind([length(y),length(x)],row,col);
    nn = N(gcpind);
    
    n = n(nn);
    
    % skip if too few
    if length(n) < 4;
        fprintf('%d overlapping points, too few,skipping\n',sum(nn));
        return;
    end
    
    %subsample GCPs to near maximum # for speed/memory
    maxGcps=100;
    if length(n) > maxGcps
        n=n(1:floor(length(n)/ maxGcps):length(n));
    end
    
    %% Registration with subsetting
    % For fitting to control, we can either load the whole image or load
    % subsets of the image in the neighborhood of each control point with
    % a size large enough to allow for shifting within the expected image
    % registration error. The former is fastest for larger numbers of GCPs
    % and or smaller image sizes. We can select by calculating the total
    % number of pixels to be loaded by each method and use the less.
    
    % subset size expect +/- 20m maximum image displacement
    res=diff(x(1:2));
    dd=20/res;
    
    [xsub,ysub,zsub]=pointAreaSubsets(gcp.x(n),gcp.y(n),x,y,m,dd,N);
    
    % send to registration fx
    [dtrans,dzall] = ...
        registerDEM2LIDAR(xsub,ysub,zsub,gcp.x(n),gcp.y(n),gcp.z(n));
    
    %% Apply registration
    ztemp = applyRegistration(dtrans,m,z,N);
    
    n=isnan(z) & ~isnan(ztemp);
    z(n) = ztemp(n);
    
    clear ztemp n
    
    m1.dtrans{i} = dtrans;
    m1.dzall{i} = dzall;
    
end



m1 = matfile(outname,'Writable',true);
m1.x=x;
m1.y=y;
m1.z=z;



