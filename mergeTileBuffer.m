function mergeTileBuffer(f0,f1,varargin)
% mergeTileBuffer: merge the edge buffers two tiles with linear feathering
%
% mergeTileBuffer(f0,f1) where f0 and f1 are the file names or mat file
% handles of neighboring tiles.

% top = 1, Bottom = 2, left = 3, right = 4;

% Maximum distance from the center of the tile overlap zone toward the
% tile edge that can be feathered, in meters.
max_feather_halfwidth = 100;

% Minimum distance from the center of the tile overlap zone toward the
% tile edge that must be feathered, in meters.
% If feathering criteria are not met at the minimum distance,
% the overlap will not be feathered.
min_feather_halfwidth = 80;

%% first test if input args are either valid filenames or mat file handles
fprintf('Attempting to merge tile %s with %s\n', f0, f1)

if isstr(f0) % it's a string, might be a filename
    if exist(f0,'file') % yup its a file
        m0 = matfile(f0); % load it
    else
        error('File does not exist');
    end
elseif isvalid(f0) % not a string, is it a valid file handle?
    m0 = f0; % yes, so set m to f and get the filename for f
    f0 = m0.Properties.Source;
else
    error('input arg must be a filename or valid matfile handle')
end

if isstr(f1) % it's a string, might be a filename
    if exist(f1,'file') % yup its a file
        m1 = matfile(f1); % load it
    else
        error('File does not exist');
    end
elseif isvalid(f1) % not a string, is it a valid file handle?
    m1 = f1; % yes, so set m to f and get the filename for f
    f1 = m1.Properties.Source;
else
    error('input arg must be a filename or valid matfile handle')
end

overwriteBuffer=false;
if ~isempty(varargin)
    if any(strcmpi(varargin,'overwrite'))
        overwriteBuffer=true;
    end
end


% check if tiles are locked by concurrent merge process
if lock_tiles(f0, f1, 'check')
    error('one or both tiles are locked by concurrent merge process, try again later')
end
% ensure lock tiles created during this process will be removed upon cleanup
cleanup_locks = onCleanup(@()lock_tiles(f0, f1, 'unlock'));
% lock tiles for this merge process
if ~lock_tiles(f0, f1, 'lock')
    error('failed to create tile lock files for this merge process')
end


% make sure writeable
m0.Properties.Writable = true;
m1.Properties.Writable = true;

%% crop tiles to buffer
c0 = m0.x >= min(m1.x) & m0.x <= max(m1.x);
r0 = m0.y >= min(m1.y) & m0.y <= max(m1.y);
c0 = [find(c0,1,'first'),find(c0,1,'last')];
r0 = [find(r0,1,'first'),find(r0,1,'last')];

if isempty(c0) || isempty(r0); error('no overlap between tiles'); end

c1 = m1.x >= min(m0.x) & m1.x <= max(m0.x);
r1 = m1.y >= min(m0.y) & m1.y <= max(m0.y);
c1 = [find(c1,1,'first'),find(c1,1,'last')];
r1 = [find(r1,1,'first'),find(r1,1,'last')];

%% Make weight array based on boundary being merged
info0=whos(m0,'z');
sz0=info0.size;
varlist0 = who(m0);
varlist1 = who(m1);

%if c0(1) > 1 % merging on f0's right boundary [old method can't account
%for deviations in tile sizes]
%if y1m > min(m0.y) && y1m < max(m0.y) && x1m > max(m0.x) [rev1 fails in
%cases where m0 is < 1/2 size of m1
if diff(c0) < diff(r0) && c0(2) == sz0(2)
    disp('Detected merge attempt on 1st tile RIGHT edge & 2nd tile LEFT edge')
    if  any(strcmp(varlist0,'mergedRight')) && any(strcmp(varlist1,'mergedLeft')) && ~overwriteBuffer
        if m0.mergedRight && m1.mergedLeft
            disp('already merged');
            return;
        end
    end

    n0 = 4;
    n1 = 3;

    [w_r0,w_c0,w_r1,w_c1,failure] = adjustFeatherZone(m0,m1,r0,c0,r1,c1,'right','left',n0,n1,max_feather_halfwidth,min_feather_halfwidth);
    if failure
        error('adjustFeatherZone failure, will not merge these edges')
    end

    W0 = [ones(1,w_c0(1)-c0(1)), linspace(1,0,diff(w_c0)+1), zeros(1,c0(2)-w_c0(2))];
    W0 = repmat(W0,diff(r0)+1,1);

    m0.mergedRight=true;
    m1.mergedLeft=true;

%elseif c0(2) < sz0(2) % merging on f0's left boundary
%elseif y1m > min(m0.y) && y1m < max(m0.y) && x1m < min(m0.x)
elseif diff(c0) < diff(r0) && c0(1) == 1
    disp('Detected merge attempt on 1st tile LEFT edge & 2nd tile RIGHT edge')
    if  any(strcmp(varlist0,'mergedLeft')) && any(strcmp(varlist1,'mergedRight')) && ~overwriteBuffer
        if m0.mergedLeft && m1.mergedRight
            disp('already merged');
            return;
        end
    end

    n0 = 3;
    n1 = 4;

    [w_r0,w_c0,w_r1,w_c1,failure] = adjustFeatherZone(m0,m1,r0,c0,r1,c1,'left','right',n0,n1,max_feather_halfwidth,min_feather_halfwidth);
    if failure
        error('adjustFeatherZone failure, will not merge these edges')
    end

    W0 = [zeros(1,w_c0(1)-c0(1)), linspace(0,1,diff(w_c0)+1), ones(1,c0(2)-w_c0(2))];
    W0 = repmat(W0,diff(r0)+1,1);

    m0.mergedLeft=true;
    m1.mergedRight=true;

%elseif r0(2) < sz0(1) % merging on f0's top boundary
%elseif x1m > min(m0.x) && x1m < max(m0.x) && y1m > max(m0.y)
elseif diff(c0) > diff(r0) && r0(1) == 1
    disp('Detected merge attempt on 1st tile TOP edge & 2nd tile BOTTOM edge')
    if  any(strcmp(varlist0,'mergedTop')) && any(strcmp(varlist1,'mergedBottom')) && ~overwriteBuffer
        if m0.mergedTop && m1.mergedBottom
            disp('already merged');
            return;
        end
    end

    n0 = 1;
    n1 = 2;

    [w_r0,w_c0,w_r1,w_c1,failure] = adjustFeatherZone(m0,m1,r0,c0,r1,c1,'top','bottom',n0,n1,max_feather_halfwidth,min_feather_halfwidth);
    if failure
        error('adjustFeatherZone failure, will not merge these edges')
    end

    W0 = [zeros(1,w_r0(1)-r0(1)), linspace(0,1,diff(w_r0)+1), ones(1,r0(2)-w_r0(2))];
    W0 = repmat(W0(:),1,diff(c0)+1);

    m0.mergedTop=true;
    m1.mergedBottom=true;

%elseif r0(1) > 1 % merging on f0's bottom boundary
%elseif x1m > min(m0.x) && x1m < max(m0.x) && y1m < min(m0.y)
elseif diff(c0) > diff(r0) && r0(2) == sz0(1)
    disp('Detected merge attempt on 1st tile BOTTOM edge & 2nd tile TOP edge')
    if  any(strcmp(varlist0,'mergedBottom')) && any(strcmp(varlist1,'mergedTop')) && ~overwriteBuffer
        if m0.mergedBottom && m1.mergedTop
            disp('already merged');
            return;
        end
    end

    n0 = 2;
    n1 = 1;

    [w_r0,w_c0,w_r1,w_c1,failure] = adjustFeatherZone(m0,m1,r0,c0,r1,c1,'bottom','top',n0,n1,max_feather_halfwidth,min_feather_halfwidth);
    if failure
        error('adjustFeatherZone failure, will not merge these edges')
    end

    W0 = [ones(1,w_r0(1)-r0(1)), linspace(1,0,diff(w_r0)+1), zeros(1,r0(2)-w_r0(2))];
    W0 = repmat(W0(:),1,diff(c0)+1);

    m0.mergedBottom=true;
    m1.mergedTop=true;

else
    error('cant determine side being merged')
end

fprintf('merging tiles...\n')

%% merge z
if ~any(strcmp(varlist0,'zbuff'))
	m0.zbuff = cell(4,1);
end

if ~any(strcmp(varlist1,'zbuff'))
    m1.zbuff = cell(4,1);
end

z0=m0.zbuff(n0,1);
z0=z0{1};

if isempty(z0)
	z0= m0.z(r0(1):r0(2),c0(1):c0(2));
	z0 = real(z0);
	m0.zbuff(n0,1) = {z0};
end

z1=m1.zbuff(n1,1);
z1=z1{1};

if isempty(z1)
 	z1= m1.z(r1(1):r1(2),c1(1):c1(2));
	z1 = real(z1);
    m1.zbuff(n1,1) = {z1};
end

z=(z0.*W0)+(z1.*(1-W0));

n=isnan(z0) & ~isnan(z1);
z(n)=z1(n);
n=~isnan(z0) & isnan(z1);
z(n)=z0(n);

m0.z(r0(1):r0(2),c0(1):c0(2))=z;
m1.z(r1(1):r1(2),c1(1):c1(2))=z;
clear z0 z1 z;

%%
if ~isempty(varargin)
    if any(strcmpi(varargin,'zonly'))
        fprintf('z only flag, stopping\n')
        return;
    end
end

%% merge z_mad
if ~any(strcmp(varlist0,'z_madbuff'))
        m0.z_madbuff = cell(4,1);
end

if ~any(strcmp(varlist1,'z_madbuff'))
        m1.z_madbuff = cell(4,1);
end

z_mad0=m0.z_madbuff(n0,1);
z_mad0=z_mad0{1};

if isempty(z_mad0)
        z_mad0= m0.z_mad(r0(1):r0(2),c0(1):c0(2));
        z_mad0 = real(z_mad0);
        m0.z_madbuff(n0,1) = {z_mad0};
end

z_mad1=m1.z_madbuff(n1,1);
z_mad1=z_mad1{1};

if isempty(z_mad1)
        z_mad1= m1.z_mad(r1(1):r1(2),c1(1):c1(2));
        z_mad1 = real(z_mad1);
        m1.z_madbuff(n1,1) = {z_mad1};
end

z_mad=(z_mad0.*W0)+(z_mad1.*(1-W0));

n=isnan(z_mad0) & ~isnan(z_mad1);
z_mad(n)=z_mad1(n);
n=~isnan(z_mad0) & isnan(z_mad1);
z_mad(n)=z_mad0(n);

m0.z_mad(r0(1):r0(2),c0(1):c0(2))=z_mad;
m1.z_mad(r1(1):r1(2),c1(1):c1(2))=z_mad;
clear z_mad0 z_mad1 z_mad;

%% merge N
if ~any(strcmp(varlist0,'Nbuff'))
        m0.Nbuff = cell(4,1);
end

if ~any(strcmp(varlist1,'Nbuff'))
        m1.Nbuff = cell(4,1);
end

N0=m0.Nbuff(n0,1);
N0=N0{1};

if isempty(N0)
        N0= m0.N(r0(1):r0(2),c0(1):c0(2));
        m0.Nbuff(n0,1) = {N0};
end

N1=m1.Nbuff(n1,1);
N1=N1{1};

if isempty(N1)
        N1= m1.N(r1(1):r1(2),c1(1):c1(2));
        m1.Nbuff(n1,1) = {N1};
end

N0=single(N0);
N1=single(N1);

N=(N0.*W0)+(N1.*(1-W0));

n= N0 == 0 & N1 ~= 0;
N(n)=N1(n);
n= N0 ~= 0 & N1 == 0;
N(n)=N0(n);

N=uint8(N);

m0.N(r0(1):r0(2),c0(1):c0(2))=N;
m1.N(r1(1):r1(2),c1(1):c1(2))=N;
clear N0 N1 N;

%% merge Nmt
if ~any(strcmp(varlist0,'Nmtbuff'))
        m0.Nmtbuff = cell(4,1);
end

if ~any(strcmp(varlist1,'Nmtbuff'))
        m1.Nmtbuff = cell(4,1);
end

Nmt0=m0.Nmtbuff(n0,1);
Nmt0=Nmt0{1};

if isempty(Nmt0)
        Nmt0= m0.Nmt(r0(1):r0(2),c0(1):c0(2));
        m0.Nmtbuff(n0,1) = {Nmt0};
end

Nmt1=m1.Nmtbuff(n1,1);
Nmt1=Nmt1{1};

if isempty(Nmt1)
        Nmt1= m1.Nmt(r1(1):r1(2),c1(1):c1(2));
        m1.Nmtbuff(n1,1) = {Nmt1};
end

Nmt0=single(Nmt0);
Nmt1=single(Nmt1);

Nmt=(Nmt0.*W0)+(Nmt1.*(1-W0));

n= Nmt0 == 0 & Nmt1 ~= 0;
Nmt(n)=Nmt1(n);
n= Nmt0 ~= 0 & Nmt1 == 0;
Nmt(n)=Nmt0(n);

Nmt=uint8(Nmt);

m0.Nmt(r0(1):r0(2),c0(1):c0(2))=Nmt;
m1.Nmt(r1(1):r1(2),c1(1):c1(2))=Nmt;
clear Nmt0 Nmt1 Nmt;

%% merge tmin
if ~any(strcmp(varlist0,'tminbuff'))
        m0.tminbuff = cell(4,1);
end

if ~any(strcmp(varlist1,'tminbuff'))
        m1.tminbuff = cell(4,1);
end

tmin0=m0.tminbuff(n0,1);
tmin0=tmin0{1};

if isempty(tmin0)
        tmin0= m0.tmin(r0(1):r0(2),c0(1):c0(2));
        m0.tminbuff(n0,1) = {tmin0};
end

tmin1=m1.tminbuff(n1,1);
tmin1=tmin1{1};

if isempty(tmin1)
        tmin1= m1.tmin(r1(1):r1(2),c1(1):c1(2));
        m1.tminbuff(n1,1) = {tmin1};
end

tmin0=single(tmin0);
tmin1=single(tmin1);

tmin=(tmin0.*W0)+(tmin1.*(1-W0));

n= tmin0 == 0 & tmin1 ~= 0;
tmin(n)=tmin1(n);
n= tmin0 ~= 0 & tmin1 == 0;
tmin(n)=tmin0(n);

tmin=uint16(tmin);

m0.tmin(r0(1):r0(2),c0(1):c0(2))=tmin;
m1.tmin(r1(1):r1(2),c1(1):c1(2))=tmin;
clear tmin0 tmin1 tmin;

%% merge tmax
if ~any(strcmp(varlist0,'tmaxbuff'))
        m0.tmaxbuff = cell(4,1);
end

if ~any(strcmp(varlist1,'tmaxbuff'))
        m1.tmaxbuff = cell(4,1);
end

tmax0=m0.tmaxbuff(n0,1);
tmax0=tmax0{1};

if isempty(tmax0)
        tmax0= m0.tmax(r0(1):r0(2),c0(1):c0(2));
        m0.tmaxbuff(n0,1) = {tmax0};
end

tmax1=m1.tmaxbuff(n1,1);
tmax1=tmax1{1};

if isempty(tmax1)
        tmax1= m1.tmax(r1(1):r1(2),c1(1):c1(2));
        m1.tmaxbuff(n1,1) = {tmax1};
end

tmax0=single(tmax0);
tmax1=single(tmax1);

tmax=(tmax0.*W0)+(tmax1.*(1-W0));

n= tmax0 == 0 & tmax1 ~= 0;
tmax(n)=tmax1(n);
n= tmax0 ~= 0 & tmax1 == 0;
tmax(n)=tmax0(n);

tmax=uint16(tmax);

m0.tmax(r0(1):r0(2),c0(1):c0(2))=tmax;
m1.tmax(r1(1):r1(2),c1(1):c1(2))=tmax;
clear tmax0 tmax1 tmax;

fprintf('tile merge complete\n')

clear cleanup_locks



function [r0_out,c0_out,r1_out,c1_out,failure] = adjustFeatherZone(m0,m1,r0_in,c0_in,r1_in,c1_in,edge0,edge1,n0,n1,maxdist,mindist)

r0_out = r0_in;
c0_out = c0_in;
r1_out = r1_in;
c1_out = c1_in;
failure = false;

varlist0 = who(m0);
varlist1 = who(m1);


% If any *buff variable already exists, check if overlap area dimensions have changed.
% If overlap area dimensions have changed, reset all arrays to unfeathered states
% and clear the *buff variable arrays for the pair of merge edges.

data_array_names = {
    'z',
    'z_mad',
    'N',
    'Nmt',
    'tmin',
    'tmax',
};

if ~all(ismember(data_array_names, varlist0))
    disp(data_array_names);
    error("One or more expected data arrays do not exist in 1st tile struct")
end
if ~all(ismember(data_array_names, varlist1))
    disp(data_array_names);
    error("One or more expected data arrays do not exist in 2nd tile struct")
end

overlap_sz_0 = [diff(r0_in)+1, diff(c0_in)+1];
overlap_sz_1 = [diff(r1_in)+1, diff(c1_in)+1];

check_varlist = {varlist0, varlist1};
check_m = {m0, m1};
check_n = {n0, n1};
check_sz = {overlap_sz_0, overlap_sz_1};

% Add data arrays that might be missing from the above list
for check_idx = 1:length(check_m)
    m_struct = check_m{check_idx};
    m_varlist = check_varlist{check_idx};
    z_size = size(m_struct, 'z');
    data_array_names_same_sz = m_varlist(cellfun(@(name) isequal(size(m_struct, name), z_size), m_varlist));
    data_array_names = union(data_array_names, data_array_names_same_sz);
end

buff_array_names = cellfun(@(x) [x,'buff'], data_array_names, 'UniformOutput',false);

need_to_reset_buffers = false;
for array_idx = 1:length(buff_array_names)
    buff_array_name = buff_array_names{array_idx};

    for check_idx = 1:length(check_m)
        m_varlist = check_varlist{check_idx};
        m = check_m{check_idx};
        n = check_n{check_idx};
        overlap_sz = check_sz{check_idx};

        if ~ismember(buff_array_name, m_varlist)
            continue;
        end
        buff_arrays = getfield(m, buff_array_name);

        if ~isempty(buff_arrays)
            if ~isempty(buff_arrays{n})
                buff_array_sz = size(buff_arrays{n});

                if ~isequal(buff_array_sz, overlap_sz)
                    disp("Tile overlap area size has changed since last merge/unmerge");
%                    disp("Will reset data arrays and clear buffer arrays for this pair of edges");
                    disp("Will check edge merged status before clearing buffer arrays for this pair of edges");
                    need_to_reset_buffers = true;
                    break;
                end
            end
        end
        clear buff_arrays
    end

    if need_to_reset_buffers
        break;
    end
end

if need_to_reset_buffers
    check_varlist = {varlist0, varlist1};
    check_m = {m0, m1};
    check_edge = {edge0, edge1};
    for check_idx = 1:length(check_m)
        m_varlist = check_varlist{check_idx};
        m = check_m{check_idx};
        edge = check_edge{check_idx};
        mergedVar = ['merged', upper(edge(1)), edge(2:end)];
        if any(strcmp(m_varlist, mergedVar)) && getfield(m, mergedVar) == true
            fprintf("Tile %g %s=true, so buffer arrays will not be automatically cleared\n", check_idx, mergedVar);
            fprintf("If you are confident that this tile has not been cropped, you should separately run undoMergeTileBuffersSingle on this tile to reset its data arrays\n")
            failure = true;
        end
    end
    if failure
        return;
    end

%    fprintf("Resetting data arrays for 1st tile %s edge\n", edge0);
%    undoMergeTileBuffersSingle(m0, edge0, 'ignoreMergedVars',true)
    fprintf("Clearing buffer arrays for 1st tile %s edge\n", edge0);
    for array_idx = 1:length(buff_array_names)
        buff_array_name = buff_array_names{array_idx};
        buff_arrays = getfield(m0, buff_array_name);
        buff_arrays{n0} = [];
        setfield(m0, buff_array_name, buff_arrays);
    end

%    fprintf("Resetting data arrays for 2nd tile %s edge\n", edge1);
%    undoMergeTileBuffersSingle(m1, edge1, 'ignoreMergedVars',true)
    fprintf("Clearing buffer arrays for 2nd tile %s edge\n", edge1);
    for array_idx = 1:length(buff_array_names)
        buff_array_name = buff_array_names{array_idx};
        buff_arrays = getfield(m1, buff_array_name);
        buff_arrays{n1} = [];
        setfield(m1, buff_array_name, buff_arrays);
    end
end


% Get zbuff array for both tiles

z0 = [];
if ismember('zbuff', varlist0)
    z0 = m0.zbuff(n0,1);
    z0 = z0{1};
end
if isempty(z0)
    z0 = m0.z(r0_in(1):r0_in(2),c0_in(1):c0_in(2));
	z0 = real(z0);
end

z1 = [];
if ismember('zbuff', varlist1)
    z1 = m1.zbuff(n1,1);
    z1 = z1{1};
end
if isempty(z1)
    z1 = m1.z(r1_in(1):r1_in(2),c1_in(1):c1_in(2));
	z1 = real(z1);
end


% Shrink the feather zone to maximum feather distance

shrink_px = 0;

if ismember(edge0, {'right', 'left'})
    coord_res = m0.x(1,2) - m0.x(1,1);
    overlap_coord_min = m0.x(1,c0_in(1));
    overlap_coord_max = m0.x(1,c0_in(2));
elseif ismember(edge0, {'top', 'bottom'})
    coord_res = m0.y(1,1) - m0.y(1,2);
    overlap_coord_min = m0.y(1,r0_in(2));
    overlap_coord_max = m0.y(1,r0_in(1));
end

fprintf("Overlap coordinate range: (%g, %g), total width = %g meters\n", overlap_coord_min, overlap_coord_max, overlap_coord_max - overlap_coord_min);

overlap_halfwidth = (overlap_coord_max - overlap_coord_min) / 2;
if overlap_halfwidth < mindist
    fprintf("Overlap area halfwidth (%gm) is less than minimum feather halfwidth (%gm)\n", overlap_halfwidth, mindist);
    failure = true;
    return;
end

if overlap_halfwidth > maxdist
    shrink_dist = overlap_halfwidth - maxdist;
    shrink_px = ceil(shrink_dist / coord_res);
    if ismember(edge0, {'right', 'left'})
        c0_out = [c0_in(1)+shrink_px, c0_in(2)-shrink_px];
        c1_out = [c1_in(1)+shrink_px, c1_in(2)-shrink_px];
    elseif ismember(edge0, {'top', 'bottom'})
        r0_out = [r0_in(1)+shrink_px, r0_in(2)-shrink_px];
        r1_out = [r1_in(1)+shrink_px, r1_in(2)-shrink_px];
    end
    fprintf("Shrunk feather zone by %g pixels on each side of overlap area to be within maximum feather halfwidth of %gm\n", shrink_px, maxdist);
end

if ismember(edge0, {'right', 'left'})
    feather_coord_min = m0.x(1,c0_out(1));
    feather_coord_max = m0.x(1,c0_out(2));
elseif ismember(edge0, {'top', 'bottom'})
    feather_coord_min = m0.y(1,r0_out(2));
    feather_coord_max = m0.y(1,r0_out(1));
end

fprintf("Base feather coordinate range: (%g, %g), total width = %g meters\n", feather_coord_min, feather_coord_max, feather_coord_max - feather_coord_min);


% Shrink the feather zone as necessary to be within NoData gap on edges
max_shrink_dist = overlap_halfwidth - mindist;
max_shrink_px = floor(max_shrink_dist / coord_res);
[buff_nrows, buff_ncols] = size(z0);
if ismember(edge0, {'right', 'left'})
    v0 = [1+shrink_px, buff_ncols-shrink_px];
elseif ismember(edge0, {'top', 'bottom'})
    v0 = [1+shrink_px, buff_nrows-shrink_px];
end
shrink_it = 0;
while true
    if ismember(edge0, {'right', 'left'})
        if     all(isnan(z0(:, v0(1)))) || all(isnan(z0(:, v0(2))))
            ;
        elseif all(isnan(z1(:, v0(1)))) || all(isnan(z1(:, v0(2))))
            ;
        else
            break
        end
    elseif ismember(edge0, {'top', 'bottom'})
        if     all(isnan(z0(v0(1), :))) || all(isnan(z0(v0(2), :)))
            ;
        elseif all(isnan(z1(v0(1), :))) || all(isnan(z1(v0(2), :)))
            ;
        else
            break
        end
    end
    shrink_it = shrink_it + 1;
    shrink_px = shrink_px + 1;
    v0 = [v0(1)+1, v0(2)-1];
    if shrink_it == 1
        disp("Found all-NaN zone on at least one tile's edge");
        disp("Will shrink feather zone one pixel at a time until out of all-NaN zone");
    end
    if shrink_px > max_shrink_px
        fprintf("Cannot continue to shrink feather zone more than %g pixels, reached minimum feather halfwidth of %gm\n", max_shrink_px, mindist);
        failure = true;
        return
    elseif v0(1) >= v0(2)
        fprintf("Cannot continue to shrink feather zone more than %g pixels, reached center of overlap area\n", max_shrink_px);
        failure = true;
        return
    end
end

if shrink_it > 0 && ~failure
    fprintf("Shrunk feather zone a total of %g pixels on each side of overlap area to be within NaN zones on tile edges\n", shrink_px);

    if ismember(edge0, {'right', 'left'})
        c0_out = [c0_in(1)+shrink_px, c0_in(2)-shrink_px];
        c1_out = [c1_in(1)+shrink_px, c1_in(2)-shrink_px];
    elseif ismember(edge0, {'top', 'bottom'})
        r0_out = [r0_in(1)+shrink_px, r0_in(2)-shrink_px];
        r1_out = [r1_in(1)+shrink_px, r1_in(2)-shrink_px];
    end

    if ismember(edge0, {'right', 'left'})
        feather_coord_min = m0.x(1,c0_out(1));
        feather_coord_max = m0.x(1,c0_out(2));
    elseif ismember(edge0, {'top', 'bottom'})
        feather_coord_min = m0.y(1,r0_out(2));
        feather_coord_max = m0.y(1,r0_out(1));
    end

    fprintf("Reduced feather coordinate range: (%g, %g), total width = %g meters\n", feather_coord_min, feather_coord_max, feather_coord_max - feather_coord_min);
end



function [locked] = lock_tiles(f0, f1, action)
    
locked = false;

f0_lockfile = [f0, '.lock'];
f1_lockfile = [f1, '.lock'];

lockfile_list = {f0_lockfile, f1_lockfile};

for lockfile_idx=1:length(lockfile_list)
    lockfile = lockfile_list{lockfile_idx};

    if isfile(lockfile)
        if strcmp(action, 'check')
            locked = true;
            return;
        elseif strcmp(action, 'unlock')
            delete(lockfile);
        end
    else
        if strcmp(action, 'lock')
            f = fopen(lockfile, 'w');
            fclose(f);
            locked = true;
        end 
    end
end
