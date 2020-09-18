function  [dZ,dX,dY,global_ids] = adjustOffsets(offsets,varargin)
%adjustOffsets perform WLSQ adjustment on 3D coregisration offsets

%[dZ,dX,dY,dem_index] = adjustOffsets(offsets) returns the optimal shifts in x,y,z 
% given the offsets. Inarg "offsets" is a structure with fields:
% 
%                     i: [n×1 double] dem 1 index
%                     j: [n×1 double] dem 2 index
%                    dz: [n×1 double] z offset (dem 1 - dem 2)
%                    dx: [n×1 double] x offset (dem 1 - dem 2)
%                    dy: [n×1 double] y offset (dem 1 - dem 2)
%                   dze: [n×1 double] z offset 1-sigma error
%                   dxe: [n×1 double] x offset 1-sigma error
%                   dye: [n×1 double] y offset 1-sigma error
%         mean_dz_coreg: [n×1 double] mean diff in z after corgestration
%       median_dz_coreg: [n×1 double] median diff in z after corgestration
%        sigma_dz_coreg: [n×1 double] std dev of diff in z after corgestration
%
% [...] = adjustOffsets(offsets,'parameter',value) specifies filter values
% for ignoring pairwise offsets. Paremeters and defaults are:
% 'offsetDiffMax',20  (filter vertical offsets larger than value)
% 'offsetErrMax',0.1  (filter offsets with errors larger than value)
% 'min_sigma_dz_coregMax', 4  (filter dems with minimum std dev's more than value)
% 'min_abs_mean_dz_coregMax, 0.1 (filter dems with minimum absolute mean
%                              post-coregistration offsets more than value )
% min_abs_median_dz_coregMax = 1 (filter dems with minimum absolute median
%                             post-coregistration offsets more than value )
%
%

% pairwise coregistration statistics filter threshold defaults
offsetDiffMax= 20;
offsetErrMax = 0.1; % use 2 m for setsm
min_sigma_dz_coregMax=4;
min_abs_mean_dz_coregMax=0.1;
min_abs_median_dz_coregMax = 1;

if length(varargin) > 1
    
    varargin(~cellfun(@isnumeric,varargin))=...
        lower(varargin(~cellfun(@isnumeric,varargin)));
    
    narg=find(strcmpi('offsetdiffmax',varargin));
    if narg
        offsetDiffMax=varargin{narg+1};
    end
    
    narg=find(strcmpi('offseterrmax',varargin));
    if narg
        offsetErrMax=varargin{narg+1};
    end
    
    narg=find(strcmpi('min_sigma_dz_coregmax',varargin));
    if narg
        min_sigma_dz_coregMax=varargin{narg+1};
    end
    
    narg=find(strcmpi('min_abs_mean_dz_coregmax',varargin));
    if narg
        min_abs_mean_dz_coregMax=varargin{narg+1};
    end
    
    narg=find(strcmpi('min_abs_median_dz_coregmax',varargin));
    if narg
        min_abs_median_dz_coregMax=varargin{narg+1};
    end
    
end

% make sure all fields are column vectors
offsets=structfun( @(x) x(:),offsets,'UniformOutput',false);

% covert global i,j to local i,j
% make sorted list of unique ids with index.
[global_ids,~,c] = unique([offsets.i;offsets.j]);
% indices for i and j colimns
ci = c(1:length(offsets.i));
cj = c(length(ci)+1:end);

% local_ids are just a monotonic list
local_ids = (1:length(global_ids))';

% save global indices
offsets.i0 = offsets.i;
offsets.j0 = offsets.j;

% populate i and j lists with local idices
offsets.i= local_ids(ci);
offsets.j = local_ids(cj);

% get the number of DEMs in the stack from the max j of i-j pairs
Ndems=length(global_ids);

% initialize output
dZ = nan(Ndems,1);
dX = nan(Ndems,1);
dY = nan(Ndems,1);


%calculate minumum pairwise errors for each dem
min_sigma_dz_coreg =  accumarray([offsets.i;offsets.j],[offsets.sigma_dz_coreg;offsets.sigma_dz_coreg],[],@min);
min_abs_mean_dz_coreg =  accumarray([offsets.i;offsets.j],abs([offsets.mean_dz_coreg;offsets.mean_dz_coreg]),[],@min);
min_abs_median_dz_coreg =  accumarray([offsets.i;offsets.j],abs([offsets.median_dz_coreg;offsets.median_dz_coreg]),[],@min);

% apply threshold
bad_dems = find(min_sigma_dz_coreg > min_sigma_dz_coregMax | ...
    min_abs_mean_dz_coreg > min_abs_mean_dz_coregMax | ...
    min_abs_median_dz_coreg > min_abs_median_dz_coregMax);

% remove all pairs that include these DEMs
n = ~ismember(offsets.i,bad_dems) & ~ismember(offsets.j,bad_dems);

% remove pairs missing offsets and filter high errors
n = n & ~isnan(offsets.dz);

% filter for maximum offsets and offset errors
nz = n & abs(offsets.dz) < offsetDiffMax & abs(offsets.dze) < offsetErrMax;
nx = n & abs(offsets.dx) < offsetDiffMax & abs(offsets.dxe) < offsetErrMax;
ny = n & abs(offsets.dy) < offsetDiffMax & abs(offsets.dye) < offsetErrMax;

% apply filter to structures
offsetsz=structfun( @(x) x(nz,:), offsets,'uniformoutput',0);
offsetsx=structfun( @(x) x(nx,:), offsets,'uniformoutput',0);
offsetsy=structfun( @(x) x(ny,:), offsets,'uniformoutput',0);

%% Adjust z offsets
Npairs = length(offsetsz.dz);

if Npairs > 0
    % find DEMs with no acceptable pairs
    i_missing = setdiff(local_ids,unique(offsetsz.i));
    j_missing = setdiff(local_ids,unique(offsetsz.j));
    n_missing = intersect(i_missing,j_missing);
    
    % Build design and weight matrices
    A = zeros(Npairs,Ndems); % initialize design matrix
    
    linearInd = sub2ind([Npairs Ndems], (1:Npairs)', offsetsz.i);
    A(linearInd) = 1;
    linearInd = sub2ind([Npairs Ndems], (1:Npairs)', offsetsz.j);
    A(linearInd) = -1;
    
    % remove filtered dems
    A(:,n_missing) = [];
    
    % add delta=0 lines
    A = [A;diag(ones(1,size(A,2)))];
    
    dz = [offsetsz.dz;zeros(size(A,2),1)];
    dze = [offsetsz.dze; ones(size(A,2),1).*4];
    wz = 1./dze.^2;
    n = local_ids;
    n(n_missing) = [];
    dZ(n) = (wz.*A)\(wz.*dz);
end

%% Adjust x offsets
Npairs = length(offsetsx.dx);

if Npairs > 0
    % find DEMs with no acceptable pairs
    i_missing = setdiff(local_ids,unique(offsetsx.i));
    j_missing = setdiff(local_ids,unique(offsetsx.j));
    n_missing = intersect(i_missing,j_missing);
    
    % Build design and weight matrices
    A = zeros(Npairs,Ndems); % initialize design matrix
    
    linearInd = sub2ind([Npairs Ndems], (1:Npairs)', offsetsx.i);
    A(linearInd) = 1;
    linearInd = sub2ind([Npairs Ndems], (1:Npairs)', offsetsx.j);
    A(linearInd) = -1;
    
    % remove filtered dems
    A(:,n_missing) = [];
    
    % add delta=0 lines
    A = [A;diag(ones(1,size(A,2)))];
    dx = [offsetsx.dx;zeros(size(A,2),1)];
    dxe = [offsetsx.dxe; ones(size(A,2),1).*4];
    wx = 1./dxe.^2;
    n = local_ids;
    n(n_missing) = [];
    dX(n) = (wx.*A)\(wx.*dx);
end

%% Adjust y offsets
Npairs = length(offsetsy.dy);

if Npairs > 0
    % find DEMs with no acceptable pairs
    i_missing = setdiff(local_ids,unique(offsetsy.i));
    j_missing = setdiff(local_ids,unique(offsetsy.j));
    n_missing = intersect(i_missing,j_missing);
    
    % Build design and weight matrices
    A = zeros(Npairs,Ndems); % initialize design matrix
    
    linearInd = sub2ind([Npairs Ndems], (1:Npairs)', offsetsy.i);
    A(linearInd) = 1;
    linearInd = sub2ind([Npairs Ndems], (1:Npairs)', offsetsy.j);
    A(linearInd) = -1;
    
    % remove filtered dems
    A(:,n_missing) = [];
    
    % add delta=0 lines
    A = [A;diag(ones(1,size(A,2)))];
    
    dy = [offsetsy.dy;zeros(size(A,2),1)];
    dye = [offsetsy.dye; ones(size(A,2),1).*4];
    wy = 1./dye.^2;
    n = local_ids;
    n(n_missing) = [];
    dY(n) = (wy.*A)\(wy.*dy);
    
end