function out=coregisterStack(varargin)
% coregisterStack coregister DEM layers in 3D array
%
% offsets=coregisterStack(x,y,z,mask,strip_ind)registered each z(:,:,i) to
% every other, ignoring zero pixels in 2D array mask. Pairs with the same
% value of strip_ind (vector of length size(z,3)) will be skipped.
% x,y,z,mask,strip_ind,t,dzdt

%parse requd argins
x = varargin{1};
y = varargin{2};
z = varargin{3};

% set up/read optional argins
mask = [];
strip_ind = [];
t = [];
dzdt = [];
if nargin > 3; mask = varargin{4}; end
if nargin > 4; strip_ind  = varargin{5}; end
if nargin > 5; t = varargin{6}; end
if nargin > 6; dzdt  = varargin{7}; end

% set missing opt argins
if isempty(mask); mask = true(length(y),length(x)); end
if isempty(strip_ind); strip_ind = 1:size(z,3); end

N = size(z,3); % # of strips

out.i=NaN;
out.j=NaN;
out.dz = NaN;
out.dx = NaN;
out.dy = NaN;

out.dze = NaN;
out.dxe = NaN;
out.dye = NaN;

out.mean_dz_uncoreg=NaN;
out.median_dz_uncoreg=NaN;
out.sigma_dz_uncoreg=NaN;

out.mean_dz_coreg=NaN;
out.median_dz_coreg=NaN;
out.sigma_dz_coreg=NaN;

if N < 2 || length(unique(strip_ind)) < 2
    fprintf('less than two DEMs in the stack, skipping\n')
    return
end

I = nchoosek(1:N,2); % indices of unique pairs

i=I(:,1); % first strip in pair
j=I(:,2); % second strip in pair

clear I

% remove pairs of segments from same strip
n = strip_ind(i) == strip_ind(j);
i(n) = [];
j(n) = [];

Npairs=size(i,1); % # of pairs

dz = nan(Npairs,1); % initialize pair offset vector
dx = nan(Npairs,1); % initialize pair offset vector
dy = nan(Npairs,1); % initialize pair offset vector

dze = nan(Npairs,1); % initialize pair offset sigma vectors
dxe = nan(Npairs,1); % initialize pair offset sigma vectors
dye = nan(Npairs,1); % initialize pair offset sigma vectors

mean_dz_uncoreg =  nan(Npairs,1);
median_dz_uncoreg =  nan(Npairs,1);
sigma_dz_uncoreg =  nan(Npairs,1);

mean_dz_coreg =  nan(Npairs,1);
median_dz_coreg =  nan(Npairs,1);
sigma_dz_coreg =  nan(Npairs,1);

% Pair coregistration loop
for pair_n=1:Npairs
    
 %  fprintf('i:%d, j:%d,pair %d of %d\n',i(pair_n),j(pair_n),pair_n,Npairs)
 
    % apply time seperation dependent ice mask
    pairMask = mask;
    
    if ~isempty(t) && ~isempty(dzdt)
        %time difference
        dt =  t(i(pair_n))-t(j(pair_n));
        
        % predicted displacement
        delz = dt.*dzdt;
        
        % if displacement > 1 m, don't use to coregister
        pairMask(abs(delz) > 1) = false;
        
        if ~any(pairMask(:))
            continue
        end
        
    end
    
    % get overlap stats before coregistration
    p = z(:,:,i(pair_n)) - z(:,:,j(pair_n));
    p(~pairMask) = NaN;
    mean_dz_uncoreg(pair_n) = nanmean(p(:));
    
    if isnan(mean_dz_uncoreg(pair_n))
        continue
    end
    
    median_dz_uncoreg(pair_n) = nanmedian(p(:));
    sigma_dz_uncoreg(pair_n) = nanstd(p(:));
    
    % coregister pair
    [zj,p,perr] = coregisterdems(x,y,z(:,:,i(pair_n)),x,y,z(:,:,j(pair_n)),pairMask);
    
    dz(pair_n) = p(1);
    dx(pair_n) = p(2);
    dy(pair_n) = p(3);
    
    dze(pair_n) = perr(1);
    dxe(pair_n) = perr(2);
    dye(pair_n) = perr(3);

    if ~isnan(dz(pair_n))
        p = z(:,:,i(pair_n)) - zj;
        p(~pairMask) = NaN;
        mean_dz_coreg(pair_n) = nanmean(p(:));
        median_dz_coreg(pair_n) = nanmedian(p(:));
        sigma_dz_coreg(pair_n) = nanstd(p(:));
    end
end

out.i=i;
out.j=j;
out.dz = dz;
out.dx = dx;
out.dy = dy;

out.dze = dze;
out.dxe = dxe;
out.dye = dye;

out.mean_dz_uncoreg=mean_dz_uncoreg;
out.median_dz_uncoreg=median_dz_uncoreg;
out.sigma_dz_uncoreg=sigma_dz_uncoreg;

out.mean_dz_coreg=mean_dz_coreg;
out.median_dz_coreg=median_dz_coreg;
out.sigma_dz_coreg=sigma_dz_coreg;