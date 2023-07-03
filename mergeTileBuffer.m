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
        error('Tile 1 file does not exist');
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
        error('Tile 2 file does not exist');
    end
elseif isvalid(f1) % not a string, is it a valid file handle?
    m1 = f1; % yes, so set m to f and get the filename for f
    f1 = m1.Properties.Source;
else
    error('input arg must be a filename or valid matfile handle')
end

overwriteBuffer=false;
preciseCorners=true;
cornersSource='adj_side';
if ~isempty(varargin)
    if any(strcmpi(varargin,'overwrite'))
        overwriteBuffer=true;
    end
    if any(strcmpi(varargin,'preciseCorners'))
        preciseCorners=true;
    end
    if any(strcmpi(varargin,'cornersSourceSameSide'))
        cornersSource='same_side';
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


[~,fname,fext] = fileparts(f0);
f0_fname = [fname,fext];
[~,fname,fext] = fileparts(f1);
f1_fname = [fname,fext];

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

n_top = 1;
n_bottom = 2;
n_left = 3;
n_right = 4;

n_top_left = 1;
n_top_right = 2;
n_bottom_right = 3;
n_bottom_left = 4;

%if c0(1) > 1 % merging on f0's right boundary [old method can't account
%for deviations in tile sizes]
%if y1m > min(m0.y) && y1m < max(m0.y) && x1m > max(m0.x) [rev1 fails in
%cases where m0 is < 1/2 size of m1
if diff(c0) < diff(r0) && c0(2) == sz0(2)
    disp('Detected merge attempt on 1st tile RIGHT edge & 2nd tile LEFT edge')
    if  any(strcmp(varlist0,'mergedRight')) && any(strcmp(varlist1,'mergedLeft')) && ~overwriteBuffer
        rightNeedsRemerge = any(strcmp(varlist0,'rightNeedsRemerge')) && m0.rightNeedsRemerge;
        leftNeedsRemerge = any(strcmp(varlist1,'leftNeedsRemerge')) && m1.leftNeedsRemerge;
        if m0.mergedRight && m1.mergedLeft && ~(rightNeedsRemerge || leftNeedsRemerge)
            disp('already merged');
            return;
        end
    end

    n0 = n_right;
    n1 = n_left;

    [w_r0,w_c0,w_r1,w_c1,failure,failure_skip_merge,failure_merge_anyway] = adjustFeatherZone(m0,m1,r0,c0,r1,c1,'right','left',n0,n1,max_feather_halfwidth,min_feather_halfwidth);
    if failure
        if failure_skip_merge
            fprintf('Skipping merge of tiles due to acceptable error: %s, %s\n', f0_fname, f1_fname);
            return;
        elseif failure_merge_anyway
            fprintf('Passing over acceptable error and merging tiles: %s, %s\n', f0_fname, f1_fname);
        else
            error('adjustFeatherZone failure is not an acceptable error, raising process-killing error')
        end
    end

    W0 = [ones(1,w_c0(1)-c0(1)), linspace(1,0,diff(w_c0)+1), zeros(1,c0(2)-w_c0(2))];
    W0 = repmat(W0,diff(r0)+1,1);

%elseif c0(2) < sz0(2) % merging on f0's left boundary
%elseif y1m > min(m0.y) && y1m < max(m0.y) && x1m < min(m0.x)
elseif diff(c0) < diff(r0) && c0(1) == 1
    disp('Detected merge attempt on 1st tile LEFT edge & 2nd tile RIGHT edge')
    if  any(strcmp(varlist0,'mergedLeft')) && any(strcmp(varlist1,'mergedRight')) && ~overwriteBuffer
        leftNeedsRemerge = any(strcmp(varlist0,'leftNeedsRemerge')) && m0.leftNeedsRemerge;
        rightNeedsRemerge = any(strcmp(varlist1,'rightNeedsRemerge')) && m1.rightNeedsRemerge;
        if m0.mergedLeft && m1.mergedRight && ~(leftNeedsRemerge || rightNeedsRemerge)
            disp('already merged');
            return;
        end
    end

    n0 = n_left;
    n1 = n_right;

    [w_r0,w_c0,w_r1,w_c1,failure,failure_skip_merge,failure_merge_anyway] = adjustFeatherZone(m0,m1,r0,c0,r1,c1,'left','right',n0,n1,max_feather_halfwidth,min_feather_halfwidth);
    if failure
        if failure_skip_merge
            fprintf('Skipping merge of tiles due to acceptable error: %s, %s\n', f0_fname, f1_fname);
            return;
        elseif failure_merge_anyway
            fprintf('Passing over acceptable error and merging tiles: %s, %s\n', f0_fname, f1_fname);
        else
            error('adjustFeatherZone failure is not an acceptable error, raising process-killing error')
        end
    end

    W0 = [zeros(1,w_c0(1)-c0(1)), linspace(0,1,diff(w_c0)+1), ones(1,c0(2)-w_c0(2))];
    W0 = repmat(W0,diff(r0)+1,1);

%elseif r0(2) < sz0(1) % merging on f0's top boundary
%elseif x1m > min(m0.x) && x1m < max(m0.x) && y1m > max(m0.y)
elseif diff(c0) > diff(r0) && r0(1) == 1
    disp('Detected merge attempt on 1st tile TOP edge & 2nd tile BOTTOM edge')
    if  any(strcmp(varlist0,'mergedTop')) && any(strcmp(varlist1,'mergedBottom')) && ~overwriteBuffer
        topNeedsRemerge = any(strcmp(varlist0,'topNeedsRemerge')) && m0.topNeedsRemerge;
        bottomNeedsRemerge = any(strcmp(varlist1,'bottomNeedsRemerge')) && m1.bottomNeedsRemerge;
        if m0.mergedTop && m1.mergedBottom && ~(topNeedsRemerge || bottomNeedsRemerge)
            disp('already merged');
            return;
        end
    end

    n0 = n_top;
    n1 = n_bottom;

    [w_r0,w_c0,w_r1,w_c1,failure,failure_skip_merge,failure_merge_anyway] = adjustFeatherZone(m0,m1,r0,c0,r1,c1,'top','bottom',n0,n1,max_feather_halfwidth,min_feather_halfwidth);
    if failure
        if failure_skip_merge
            fprintf('Skipping merge of tiles due to acceptable error: %s, %s\n', f0_fname, f1_fname);
            return;
        elseif failure_merge_anyway
            fprintf('Passing over acceptable error and merging tiles: %s, %s\n', f0_fname, f1_fname);
        else
            error('adjustFeatherZone failure is not an acceptable error, raising process-killing error')
        end
    end

    W0 = [zeros(1,w_r0(1)-r0(1)), linspace(0,1,diff(w_r0)+1), ones(1,r0(2)-w_r0(2))];
    W0 = repmat(W0(:),1,diff(c0)+1);

%elseif r0(1) > 1 % merging on f0's bottom boundary
%elseif x1m > min(m0.x) && x1m < max(m0.x) && y1m < min(m0.y)
elseif diff(c0) > diff(r0) && r0(2) == sz0(1)
    disp('Detected merge attempt on 1st tile BOTTOM edge & 2nd tile TOP edge')
    if  any(strcmp(varlist0,'mergedBottom')) && any(strcmp(varlist1,'mergedTop')) && ~overwriteBuffer
        bottomNeedsRemerge = any(strcmp(varlist0,'bottomNeedsRemerge')) && m0.bottomNeedsRemerge;
        topNeedsRemerge = any(strcmp(varlist1,'topNeedsRemerge')) && m1.topNeedsRemerge;
        if m0.mergedBottom && m1.mergedTop && ~(bottomNeedsRemerge || topNeedsRemerge)
            disp('already merged');
            return;
        end
    end

    n0 = n_bottom;
    n1 = n_top;

    [w_r0,w_c0,w_r1,w_c1,failure,failure_skip_merge,failure_merge_anyway] = adjustFeatherZone(m0,m1,r0,c0,r1,c1,'bottom','top',n0,n1,max_feather_halfwidth,min_feather_halfwidth);
    if failure
        if failure_skip_merge
            fprintf('Skipping merge of tiles due to acceptable error: %s, %s\n', f0_fname, f1_fname);
            return;
        elseif failure_merge_anyway
            fprintf('Passing over acceptable error and merging tiles: %s, %s\n', f0_fname, f1_fname);
        else
            error('adjustFeatherZone failure is not an acceptable error, raising process-killing error')
        end
    end

    W0 = [ones(1,w_r0(1)-r0(1)), linspace(1,0,diff(w_r0)+1), zeros(1,r0(2)-w_r0(2))];
    W0 = repmat(W0(:),1,diff(c0)+1);

else
    error('cant determine side being merged')
end

fprintf('merging tiles...\n')

W0_feather_zone = (W0 > 0) & (W0 < 1);

%% merge z
z0 = loadBuffArray(m0,n0,r0,c0,'z',true);
z1 = loadBuffArray(m1,n1,r1,c1,'z',true);

if preciseCorners
    cacheOriginalCorners(m0, n0, 'z', cornersSource);
    cacheOriginalCorners(m1, n1, 'z', cornersSource);
    if n0 == n_right || n0 == n_left
        z0 = resetBuffCorners(m0, n0, 'z', z0);
        z1 = resetBuffCorners(m1, n1, 'z', z1);
    end
end

z=(z0.*W0)+(z1.*(1-W0));

n=isnan(z0) & ~isnan(z1);
z(n)=z1(n);
n=~isnan(z0) & isnan(z1);
z(n)=z0(n);

m0.z(r0(1):r0(2),c0(1):c0(2))=z;
m1.z(r1(1):r1(2),c1(1):c1(2))=z;

if preciseCorners
    if n0 == n_right || n0 == n_left
        burnCornersIntoAdjBuffs(m0, n0, 'z', z);
        burnCornersIntoAdjBuffs(m1, n1, 'z', z);
    end
end

clear z0 z1 z;

%%
if ~isempty(varargin)
    if any(strcmpi(varargin,'zonly'))
        fprintf('z only flag, stopping\n')
        return;
    end
end

%% merge z_mad
z_mad0 = loadBuffArray(m0,n0,r0,c0,'z_mad',true);
z_mad1 = loadBuffArray(m1,n1,r1,c1,'z_mad',true);

if preciseCorners
    cacheOriginalCorners(m0, n0, 'z_mad', cornersSource);
    cacheOriginalCorners(m1, n1, 'z_mad', cornersSource);
    if n0 == n_right || n0 == n_left
        z_mad0 = resetBuffCorners(m0, n0, 'z_mad', z_mad0);
        z_mad1 = resetBuffCorners(m1, n1, 'z_mad', z_mad1);
    end
end

z_mad=(z_mad0.*W0)+(z_mad1.*(1-W0));

n=isnan(z_mad0) & ~isnan(z_mad1);
z_mad(n)=z_mad1(n);
n=~isnan(z_mad0) & isnan(z_mad1);
z_mad(n)=z_mad0(n);

m0.z_mad(r0(1):r0(2),c0(1):c0(2))=z_mad;
m1.z_mad(r1(1):r1(2),c1(1):c1(2))=z_mad;

if preciseCorners
    if n0 == n_right || n0 == n_left
        burnCornersIntoAdjBuffs(m0, n0, 'z_mad', z_mad);
        burnCornersIntoAdjBuffs(m1, n1, 'z_mad', z_mad);
    end
end

clear z_mad0 z_mad1 z_mad;

%% merge N
N0 = loadBuffArray(m0,n0,r0,c0,'N',false);
N1 = loadBuffArray(m1,n1,r1,c1,'N',false);

N0=single(N0);
N1=single(N1);

if preciseCorners
    cacheOriginalCorners(m0, n0, 'N', cornersSource);
    cacheOriginalCorners(m1, n1, 'N', cornersSource);
    if n0 == n_right || n0 == n_left
        N0 = resetBuffCorners(m0, n0, 'N', N0);
        N1 = resetBuffCorners(m1, n1, 'N', N1);
    end
end

%N = (N0.*W0) + (N1.*(1-W0));
N = max(N0, N1);
N(W0 == 1) = N0(W0 == 1);
N(W0 == 0) = N1(W0 == 0);

n= N0 == 0 & N1 ~= 0;
N(n)=N1(n);
n= N0 ~= 0 & N1 == 0;
N(n)=N0(n);

N=uint8(N);

m0.N(r0(1):r0(2),c0(1):c0(2))=N;
m1.N(r1(1):r1(2),c1(1):c1(2))=N;

if preciseCorners
    if n0 == n_right || n0 == n_left
        burnCornersIntoAdjBuffs(m0, n0, 'N', N);
        burnCornersIntoAdjBuffs(m1, n1, 'N', N);
    end
end

clear N0 N1 N;

%% merge Nmt
Nmt0 = loadBuffArray(m0,n0,r0,c0,'Nmt',false);
Nmt1 = loadBuffArray(m1,n1,r1,c1,'Nmt',false);

Nmt0=single(Nmt0);
Nmt1=single(Nmt1);

if preciseCorners
    cacheOriginalCorners(m0, n0, 'Nmt', cornersSource);
    cacheOriginalCorners(m1, n1, 'Nmt', cornersSource);
    if n0 == n_right || n0 == n_left
        Nmt0 = resetBuffCorners(m0, n0, 'Nmt', Nmt0);
        Nmt1 = resetBuffCorners(m1, n1, 'Nmt', Nmt1);
    end
end

%Nmt = (Nmt0.*W0) + (Nmt1.*(1-W0));
Nmt = max(Nmt0, Nmt1);
Nmt(W0 == 1) = Nmt0(W0 == 1);
Nmt(W0 == 0) = Nmt1(W0 == 0);

n= Nmt0 == 0 & Nmt1 ~= 0;
Nmt(n)=Nmt1(n);
n= Nmt0 ~= 0 & Nmt1 == 0;
Nmt(n)=Nmt0(n);

Nmt=uint8(Nmt);

m0.Nmt(r0(1):r0(2),c0(1):c0(2))=Nmt;
m1.Nmt(r1(1):r1(2),c1(1):c1(2))=Nmt;

if preciseCorners
    if n0 == n_right || n0 == n_left
        burnCornersIntoAdjBuffs(m0, n0, 'Nmt', Nmt);
        burnCornersIntoAdjBuffs(m1, n1, 'Nmt', Nmt);
    end
end

clear Nmt0 Nmt1 Nmt;

%% merge tmin
tmin0 = loadBuffArray(m0,n0,r0,c0,'tmin',false);
tmin1 = loadBuffArray(m1,n1,r1,c1,'tmin',false);

tmin0=single(tmin0);
tmin1=single(tmin1);

if preciseCorners
    cacheOriginalCorners(m0, n0, 'tmin', cornersSource);
    cacheOriginalCorners(m1, n1, 'tmin', cornersSource);
    if n0 == n_right || n0 == n_left
        tmin0 = resetBuffCorners(m0, n0, 'tmin', tmin0);
        tmin1 = resetBuffCorners(m1, n1, 'tmin', tmin1);
    end
end

%tmin = (tmin0.*W0) + (tmin1.*(1-W0));
tmin = min(tmin0, tmin1);
tmin(W0 == 1) = tmin0(W0 == 1);
tmin(W0 == 0) = tmin1(W0 == 0);

n= tmin0 == 0 & tmin1 ~= 0;
tmin(n)=tmin1(n);
n= tmin0 ~= 0 & tmin1 == 0;
tmin(n)=tmin0(n);

tmin=uint16(tmin);

m0.tmin(r0(1):r0(2),c0(1):c0(2))=tmin;
m1.tmin(r1(1):r1(2),c1(1):c1(2))=tmin;

if preciseCorners
    if n0 == n_right || n0 == n_left
        burnCornersIntoAdjBuffs(m0, n0, 'tmin', tmin);
        burnCornersIntoAdjBuffs(m1, n1, 'tmin', tmin);
    end
end

clear tmin0 tmin1 tmin;

%% merge tmax
tmax0 = loadBuffArray(m0,n0,r0,c0,'tmax',false);
tmax1 = loadBuffArray(m1,n1,r1,c1,'tmax',false);

tmax0=single(tmax0);
tmax1=single(tmax1);

if preciseCorners
    cacheOriginalCorners(m0, n0, 'tmax', cornersSource);
    cacheOriginalCorners(m1, n1, 'tmax', cornersSource);
    if n0 == n_right || n0 == n_left
        tmax0 = resetBuffCorners(m0, n0, 'tmax', tmax0);
        tmax1 = resetBuffCorners(m1, n1, 'tmax', tmax1);
    end
end

%tmax = (tmax0.*W0) + (tmax1.*(1-W0));
tmax = max(tmax0, tmax1);
tmax(W0 == 1) = tmax0(W0 == 1);
tmax(W0 == 0) = tmax1(W0 == 0);

n= tmax0 == 0 & tmax1 ~= 0;
tmax(n)=tmax1(n);
n= tmax0 ~= 0 & tmax1 == 0;
tmax(n)=tmax0(n);

tmax=uint16(tmax);

m0.tmax(r0(1):r0(2),c0(1):c0(2))=tmax;
m1.tmax(r1(1):r1(2),c1(1):c1(2))=tmax;

if preciseCorners
    if n0 == n_right || n0 == n_left
        burnCornersIntoAdjBuffs(m0, n0, 'tmax', tmax);
        burnCornersIntoAdjBuffs(m1, n1, 'tmax', tmax);
    end
end

clear tmax0 tmax1 tmax;

%% merge waterFillMask
if any(strcmp(varlist0,'waterFillMask')) || any(strcmp(varlist1,'waterFillMask'))

    waterFillMask0 = loadBuffArray(m0,n0,r0,c0,'waterFillMask',false);
    waterFillMask1 = loadBuffArray(m1,n1,r1,c1,'waterFillMask',false);

    waterFillMask0=logical(waterFillMask0);
    waterFillMask1=logical(waterFillMask1);

    if preciseCorners
        cacheOriginalCorners(m0, n0, 'waterFillMask', cornersSource);
        cacheOriginalCorners(m1, n1, 'waterFillMask', cornersSource);
        if n0 == n_right || n0 == n_left
            waterFillMask0 = resetBuffCorners(m0, n0, 'waterFillMask', waterFillMask0);
            waterFillMask1 = resetBuffCorners(m1, n1, 'waterFillMask', waterFillMask1);
        end
    end

    waterFillMask = waterFillMask0 | waterFillMask1;

    m0.waterFillMask(r0(1):r0(2),c0(1):c0(2))=waterFillMask;
    m1.waterFillMask(r1(1):r1(2),c1(1):c1(2))=waterFillMask;

    if preciseCorners
        if n0 == n_right || n0 == n_left
            burnCornersIntoAdjBuffs(m0, n0, 'waterFillMask', waterFillMask);
            burnCornersIntoAdjBuffs(m1, n1, 'waterFillMask', waterFillMask);
        end
    end

    clear waterFillMask0 waterFillMask1 waterFillMask;
end

fprintf('tile merge complete\n')

if n0 == n_right
    m0.mergedRight = true;
elseif n0 == n_left
    m0.mergedLeft = true;
elseif n0 == n_top
    m0.mergedTop = true;
elseif n0 == n_bottom
    m0.mergedBottom = true;
end
if n1 == n_right
    m1.mergedRight = true;
elseif n1 == n_left
    m1.mergedLeft = true;
elseif n1 == n_top
    m1.mergedTop = true;
elseif n1 == n_bottom
    m1.mergedBottom = true;
end

if preciseCorners
    if n0 == n_right
        m0.rightNeedsRemerge = false;
    elseif n0 == n_left
        m0.leftNeedsRemerge = false;
    elseif n0 == n_top
        m0.topNeedsRemerge = false;
    elseif n0 == n_bottom
        m0.bottomNeedsRemerge = false;
    end
    if n1 == n_right
        m1.rightNeedsRemerge = false;
    elseif n1 == n_left
        m1.leftNeedsRemerge = false;
    elseif n1 == n_top
        m1.topNeedsRemerge = false;
    elseif n1 == n_bottom
        m1.bottomNeedsRemerge = false;
    end
end

clear cleanup_locks;



function z0 = loadBuffArray(m0, n0, r0, c0, data_array_name, make_real)

n_top = 1;
n_bottom = 2;
n_left = 3;
n_right = 4;

varlist0 = who(m0);
data_buff_name = [data_array_name,'buff'];

if ~ismember(data_buff_name, varlist0)
	setfield(m0, data_buff_name, cell(4,1));
end

z0 = getfield(m0, data_buff_name, {n0,1});
z0 = z0{1};

if isempty(z0)
    if n0 == n_top || n0 == n_bottom
	    eval(['z0 = m0.',data_array_name,'(r0(1):r0(2),:);']);
	elseif n0 == n_left || n0 == n_right
	    eval(['z0 = m0.',data_array_name,'(:,c0(1):c0(2));']);
    end
    if make_real
	    z0 = real(z0);
	end
	setfield(m0, data_buff_name, {n0,1}, {z0});
end

if n0 == n_top || n0 == n_bottom
    z0 = z0(:,c0(1):c0(2));
elseif n0 == n_left || n0 == n_right
    z0 = z0(r0(1):r0(2),:);
end



function cacheOriginalCorners(m0, n0, data_arr_name, source)

n_top = 1;
n_bottom = 2;
n_left = 3;
n_right = 4;

n_top_left = 1;
n_top_right = 2;
n_bottom_right = 3;
n_bottom_left = 4;

if n0 == n_top
    working_tile_position = 'bottom';
elseif n0 == n_bottom
    working_tile_position = 'top';
elseif n0 == n_left
    working_tile_position = 'right';
elseif n0 == n_right
    working_tile_position = 'left';
end

data_buff_name = [data_arr_name, 'buff'];
data_corners_name = [data_arr_name, 'buffcorners'];

varlist0 = who(m0);
if ~ismember(data_corners_name, varlist0)
    setfield(m0, data_corners_name, cell(4,1));
end

if n0 == n_right
    buff_r = [];
    
    corn_arr = getfield(m0, data_corners_name, {n_top_right,1});
    corn_arr = corn_arr{1};
    if isempty(corn_arr)
        buff_t = getfield(m0, data_buff_name, {n_top,1});
        buff_t = buff_t{1};
        if ~isempty(buff_t)
            if isempty(buff_r)
                buff_r = getfield(m0, data_buff_name, {n_right,1});
                buff_r = buff_r{1};
            end
            corn_nrows = size(buff_t, 1);
            corn_ncols = size(buff_r, 2);
            if strcmp(source, 'adj_side')
                corn_source_edge = 'top';
                corn_arr = buff_t(:, (end-corn_ncols+1):end);
            elseif strcmp(source, 'same_side')
                corn_source_edge = 'right';
                corn_arr = buff_r(1:corn_nrows, :);
            end
            fprintf("Caching %s tile's '%s' array top-right original corner from %s edge %s array\n", working_tile_position, data_arr_name, corn_source_edge, data_buff_name);
            setfield(m0, data_corners_name, {n_top_right,1}, {corn_arr});
        end
        clear buff_t;
    end
    
    corn_arr = getfield(m0, data_corners_name, {n_bottom_right,1});
    corn_arr = corn_arr{1};
    if isempty(corn_arr)
        buff_b = getfield(m0, data_buff_name, {n_bottom,1});
        buff_b = buff_b{1};
        if ~isempty(buff_b)
            if isempty(buff_r)
                buff_r = getfield(m0, data_buff_name, {n_right,1});
                buff_r = buff_r{1};
            end
            corn_nrows = size(buff_b, 1);
            corn_ncols = size(buff_r, 2);
            if strcmp(source, 'adj_side')
                corn_source_edge = 'bottom';
                corn_arr = buff_b(:, (end-corn_ncols+1):end);
            elseif strcmp(source, 'same_side')
                corn_source_edge = 'right';
                corn_arr = buff_r((end-corn_nrows+1):end, :);
            end
            fprintf("Caching %s tile's '%s' array bottom-right original corner from %s edge %s array\n", working_tile_position, data_arr_name, corn_source_edge, data_buff_name);
            setfield(m0, data_corners_name, {n_bottom_right,1}, {corn_arr});
        end
        clear buff_b;
    end

elseif n0 == n_left
    buff_l = [];
    
    corn_arr = getfield(m0, data_corners_name, {n_top_left,1});
    corn_arr = corn_arr{1};
    if isempty(corn_arr)
        buff_t = getfield(m0, data_buff_name, {n_top,1});
        buff_t = buff_t{1};
        if ~isempty(buff_t)
            if isempty(buff_l)
                buff_l = getfield(m0, data_buff_name, {n_left,1});
                buff_l = buff_l{1};
            end
            corn_nrows = size(buff_t, 1);
            corn_ncols = size(buff_l, 2);
            if strcmp(source, 'adj_side')
                corn_source_edge = 'top';
                corn_arr = buff_t(:, 1:corn_ncols);
            elseif strcmp(source, 'same_side')
                corn_source_edge = 'left';
                corn_arr = buff_l(1:corn_nrows, :);
            end
            fprintf("Caching %s tile's '%s' array top-left original corner from %s edge %s array\n", working_tile_position, data_arr_name, corn_source_edge, data_buff_name);
            setfield(m0, data_corners_name, {n_top_left,1}, {corn_arr});
        end
        clear buff_t;
    end
    
    corn_arr = getfield(m0, data_corners_name, {n_bottom_left,1});
    corn_arr = corn_arr{1};
    if isempty(corn_arr)
        buff_b = getfield(m0, data_buff_name, {n_bottom,1});
        buff_b = buff_b{1};
        if ~isempty(buff_b)
            if isempty(buff_l)
                buff_l = getfield(m0, data_buff_name, {n_left,1});
                buff_l = buff_l{1};
            end
            corn_nrows = size(buff_b, 1);
            corn_ncols = size(buff_l, 2);
            if strcmp(source, 'adj_side')
                corn_source_edge = 'bottom';
                corn_arr = buff_b(:, 1:corn_ncols);
            elseif strcmp(source, 'same_side')
                corn_source_edge = 'left';
                corn_arr = buff_l((end-corn_nrows+1):end, :);
            end
            fprintf("Caching %s tile's '%s' array bottom-left original corner from %s edge %s array\n", working_tile_position, data_arr_name, corn_source_edge, data_buff_name);
            setfield(m0, data_corners_name, {n_bottom_left,1}, {corn_arr});
        end
        clear buff_b;
    end
    
elseif n0 == n_top
    buff_t = [];
    
    corn_arr = getfield(m0, data_corners_name, {n_top_right,1});
    corn_arr = corn_arr{1};
    if isempty(corn_arr)
        buff_r = getfield(m0, data_buff_name, {n_right,1});
        buff_r = buff_r{1};
        if ~isempty(buff_r)
            if isempty(buff_t)
                buff_t = getfield(m0, data_buff_name, {n_top,1});
                buff_t = buff_t{1};
            end
            corn_nrows = size(buff_t, 1);
            corn_ncols = size(buff_r, 2);
            if strcmp(source, 'adj_side')
                corn_source_edge = 'right';
                corn_arr = buff_r(1:corn_nrows, :);
            elseif strcmp(source, 'same_side')
                corn_source_edge = 'top';
                corn_arr = buff_t(:, (end-corn_ncols+1):end);
            end
            fprintf("Caching %s tile's '%s' array top-right original corner from %s edge %s array\n", working_tile_position, data_arr_name, corn_source_edge, data_buff_name);
            setfield(m0, data_corners_name, {n_top_right,1}, {corn_arr});
        end
        clear buff_r;
    end
    
    corn_arr = getfield(m0, data_corners_name, {n_top_left,1});
    corn_arr = corn_arr{1};
    if isempty(corn_arr)
        buff_l = getfield(m0, data_buff_name, {n_left,1});
        buff_l = buff_l{1};
        if ~isempty(buff_l)
            if isempty(buff_t)
                buff_t = getfield(m0, data_buff_name, {n_top,1});
                buff_t = buff_t{1};
            end
            corn_nrows = size(buff_t, 1);
            corn_ncols = size(buff_l, 2);
            if strcmp(source, 'adj_side')
                corn_source_edge = 'left';
                corn_arr = buff_l(1:corn_nrows, :);
            elseif strcmp(source, 'same_side')
                corn_source_edge = 'top';
                corn_arr = buff_t(:, 1:corn_ncols);
            end
            fprintf("Caching %s tile's '%s' array top-left original corner from %s edge %s array\n", working_tile_position, data_arr_name, corn_source_edge, data_buff_name);
            setfield(m0, data_corners_name, {n_top_left,1}, {corn_arr});
        end
        clear buff_l;
    end
    
elseif n0 == n_bottom
    buff_b = [];
    
    corn_arr = getfield(m0, data_corners_name, {n_bottom_right,1});
    corn_arr = corn_arr{1};
    if isempty(corn_arr)
        buff_r = getfield(m0, data_buff_name, {n_right,1});
        buff_r = buff_r{1};
        if ~isempty(buff_r)
            if isempty(buff_b)
                buff_b = getfield(m0, data_buff_name, {n_bottom,1});
                buff_b = buff_b{1};
            end
            corn_nrows = size(buff_b, 1);
            corn_ncols = size(buff_r, 2);
            if strcmp(source, 'adj_side')
                corn_source_edge = 'right';
                corn_arr = buff_r((end-corn_nrows+1):end, :);
            elseif strcmp(source, 'same_side')
                corn_source_edge = 'bottom';
                corn_arr = buff_b(:, (end-corn_ncols+1):end);
            end
            fprintf("Caching %s tile's '%s' array bottom-right original corner from %s edge %s array\n", working_tile_position, data_arr_name, corn_source_edge, data_buff_name);
            setfield(m0, data_corners_name, {n_bottom_right,1}, {corn_arr});
        end
        clear buff_r;
    end
    
    corn_arr = getfield(m0, data_corners_name, {n_bottom_left,1});
    corn_arr = corn_arr{1};
    if isempty(corn_arr)
        buff_l = getfield(m0, data_buff_name, {n_left,1});
        buff_l = buff_l{1};
        if ~isempty(buff_l)
            if isempty(buff_b)
                buff_b = getfield(m0, data_buff_name, {n_bottom,1});
                buff_b = buff_b{1};
            end
            corn_nrows = size(buff_b, 1);
            corn_ncols = size(buff_l, 2);
            if strcmp(source, 'adj_side')
                corn_source_edge = 'left';
                corn_arr = buff_l((end-corn_nrows+1):end, :);
            elseif strcmp(source, 'same_side')
                corn_source_edge = 'bottom';
                corn_arr = buff_b(:, 1:corn_ncols);
            end
            fprintf("Caching %s tile's '%s' array bottom-left original corner from %s edge %s array\n", working_tile_position, data_arr_name, corn_source_edge, data_buff_name);
            setfield(m0, data_corners_name, {n_bottom_left,1}, {corn_arr});
        end
        clear buff_l;
    end

end



function [buff_arr] = resetBuffCorners(m0, n0, data_arr_name, buff_arr)

n_top = 1;
n_bottom = 2;
n_left = 3;
n_right = 4;

n_top_left = 1;
n_top_right = 2;
n_bottom_right = 3;
n_bottom_left = 4;

if n0 == n_top
    working_tile_position = 'bottom';
elseif n0 == n_bottom
    working_tile_position = 'top';
elseif n0 == n_left
    working_tile_position = 'right';
elseif n0 == n_right
    working_tile_position = 'left';
end

data_buff_name = [data_arr_name, 'buff'];
data_corners_name = [data_arr_name, 'buffcorners'];

varlist0 = who(m0);
if ~ismember(data_corners_name, varlist0)
    setfield(m0, data_corners_name, cell(4,1));
    error("Tile '%s' variable has not been set", data_corners_name);
end

if ~ismember('origCornersList', varlist0)
    m0.origCornersList = struct;
end
if ~ismember(data_buff_name, fields(m0.origCornersList))
    m0.origCornersList = setfield(m0.origCornersList, data_buff_name, struct);
end
origCornersList = getfield(m0.origCornersList, data_buff_name);
origCornersVarlist = fields(origCornersList);

modified_buff_arr = false;

if n0 == n_right

    rightTopCornerIsOrig = ismember('rightTopCornerIsOrig', origCornersVarlist) && origCornersList.rightTopCornerIsOrig;
    if ~rightTopCornerIsOrig
        corn_arr = getfield(m0, data_corners_name, {n_top_right,1});
        corn_arr = corn_arr{1};
        if ~isempty(corn_arr)
            fprintf("Resetting %s tile's right-side %s array top-right corner values to original corner values\n", working_tile_position, data_buff_name);
            [corn_nrows, corn_ncols] = size(corn_arr);
            buff_arr(1:corn_nrows, :) = corn_arr;
            origCornersList.rightTopCornerIsOrig = true;
            modified_buff_arr = true;
        end
    end

    rightBottomCornerIsOrig = ismember('rightBottomCornerIsOrig', origCornersVarlist) && origCornersList.rightBottomCornerIsOrig;
    if ~rightBottomCornerIsOrig
        corn_arr = getfield(m0, data_corners_name, {n_bottom_right,1});
        corn_arr = corn_arr{1};
        if ~isempty(corn_arr)
            fprintf("Resetting %s tile's right-side %s array bottom-right corner values to original corner values\n", working_tile_position, data_buff_name);
            [corn_nrows, corn_ncols] = size(corn_arr);
            buff_arr((end-corn_nrows+1):end, :) = corn_arr;
            origCornersList.rightBottomCornerIsOrig = true;
            modified_buff_arr = true;
        end
    end

elseif n0 == n_left

    leftTopCornerIsOrig = ismember('leftTopCornerIsOrig', origCornersVarlist) && origCornersList.leftTopCornerIsOrig;
    if ~leftTopCornerIsOrig
        corn_arr = getfield(m0, data_corners_name, {n_top_left,1});
        corn_arr = corn_arr{1};
        if ~isempty(corn_arr)
            fprintf("Resetting %s tile's left-side %s array top-left corner values to original corner values\n", working_tile_position, data_buff_name);
            [corn_nrows, corn_ncols] = size(corn_arr);
            buff_arr(1:corn_nrows, :) = corn_arr;
            origCornersList.leftTopCornerIsOrig = true;
            modified_buff_arr = true;
        end
    end

    leftBottomCornerIsOrig = ismember('leftBottomCornerIsOrig', origCornersVarlist) && origCornersList.leftBottomCornerIsOrig;
    if ~leftBottomCornerIsOrig
        corn_arr = getfield(m0, data_corners_name, {n_bottom_left,1});
        corn_arr = corn_arr{1};
        if ~isempty(corn_arr)
            fprintf("Resetting %s tile's left-side %s array bottom-left corner values to original corner values\n", working_tile_position, data_buff_name);
            [corn_nrows, corn_ncols] = size(corn_arr);
            buff_arr((end-corn_nrows+1):end, :) = corn_arr;
            origCornersList.leftBottomCornerIsOrig = true;
            modified_buff_arr = true;
        end
    end

elseif n0 == n_top

    topRightCornerIsOrig = ismember('topRightCornerIsOrig', origCornersVarlist) && origCornersList.topRightCornerIsOrig;
    if ~topRightCornerIsOrig
        corn_arr = getfield(m0, data_corners_name, {n_top_right,1});
        corn_arr = corn_arr{1};
        if ~isempty(corn_arr)
            fprintf("Resetting %s tile's top-side %s array top-right corner values to original corner values\n", working_tile_position, data_buff_name);
            [corn_nrows, corn_ncols] = size(corn_arr);
            buff_arr(:, (end-corn_ncols+1):end) = corn_arr;
            origCornersList.topRightCornerIsOrig = true;
            modified_buff_arr = true;
        end
    end

    topLeftCornerIsOrig = ismember('topLeftCornerIsOrig', origCornersVarlist) && origCornersList.topLeftCornerIsOrig;
    if ~topLeftCornerIsOrig
        corn_arr = getfield(m0, data_corners_name, {n_top_left,1});
        corn_arr = corn_arr{1};
        if ~isempty(corn_arr)
            fprintf("Resetting %s tile's top-side %s array top-left corner values to original corner values\n", working_tile_position, data_buff_name);
            [corn_nrows, corn_ncols] = size(corn_arr);
            buff_arr(:, 1:corn_ncols) = corn_arr;
            origCornersList.topLeftCornerIsOrig = true;
            modified_buff_arr = true;
        end
    end

elseif n0 == n_bottom

    bottomRightCornerIsOrig = ismember('bottomRightCornerIsOrig', origCornersVarlist) && origCornersList.bottomRightCornerIsOrig;
    if ~bottomRightCornerIsOrig
        corn_arr = getfield(m0, data_corners_name, {n_bottom_right,1});
        corn_arr = corn_arr{1};
        if ~isempty(corn_arr)
            fprintf("Resetting %s tile's bottom-side %s array bottom-right corner values to original corner values\n", working_tile_position, data_buff_name);
            [corn_nrows, corn_ncols] = size(corn_arr);
            buff_arr(:, (end-corn_ncols+1):end) = corn_arr;
            origCornersList.bottomRightCornerIsOrig = true;
            modified_buff_arr = true;
        end
    end

    bottomLeftCornerIsOrig = ismember('bottomLeftCornerIsOrig', origCornersVarlist) && origCornersList.bottomLeftCornerIsOrig;
    if ~bottomLeftCornerIsOrig
        corn_arr = getfield(m0, data_corners_name, {n_bottom_left,1});
        corn_arr = corn_arr{1};
        if ~isempty(corn_arr)
            fprintf("Resetting %s tile's bottom-side %s array bottom-left corner values to original corner values\n", working_tile_position, data_buff_name);
            [corn_nrows, corn_ncols] = size(corn_arr);
            buff_arr(:, 1:corn_ncols) = corn_arr;
            origCornersList.bottomLeftCornerIsOrig = true;
            modified_buff_arr = true;
        end
    end

end

if modified_buff_arr
    fprintf("Saving resetBuffCorners result to tile '%s' variable\n", data_buff_name);
    m0.origCornersList = setfield(m0.origCornersList, data_buff_name, origCornersList);
    setfield(m0, data_buff_name, {n0,1}, {buff_arr});
end



function burnCornersIntoAdjBuffs(m0, n0, data_arr_name, buff_arr)

n_top = 1;
n_bottom = 2;
n_left = 3;
n_right = 4;

n_top_left = 1;
n_top_right = 2;
n_bottom_right = 3;
n_bottom_left = 4;

if n0 == n_top
    working_tile_position = 'bottom';
elseif n0 == n_bottom
    working_tile_position = 'top';
elseif n0 == n_left
    working_tile_position = 'right';
elseif n0 == n_right
    working_tile_position = 'left';
end

data_buff_name = [data_arr_name, 'buff'];
data_corners_name = [data_arr_name, 'buffcorners'];

varlist0 = who(m0);
if ~ismember(data_corners_name, varlist0)
    setfield(m0, data_corners_name, cell(4,1));
    error("Tile '%s' variable has not been set", data_corners_name);
end

if n0 == n_right
    buff_r = buff_arr;
    clear buff_arr;

    corn_arr = getfield(m0, data_corners_name, {n_top_right,1});
    corn_arr = corn_arr{1};
    if ~isempty(corn_arr)
        [corn_nrows, corn_ncols] = size(corn_arr);
        buff_t = getfield(m0, data_buff_name, {n_top,1});
        buff_t = buff_t{1};
        if ~isempty(buff_t)
            fprintf("Burning %s tile's right-side %s array top-right corner values into top-side array\n", working_tile_position, data_buff_name);
            fprintf("WARNING! %s tile will need to be remerged on the top side!\n", working_tile_position);
            buff_t(:, (end-corn_ncols+1):end) = buff_r(1:corn_nrows, :);
            setfield(m0, data_buff_name, {n_top,1}, {buff_t});
            m0.topNeedsRemerge = true;
        end
    end

    corn_arr = getfield(m0, data_corners_name, {n_bottom_right,1});
    corn_arr = corn_arr{1};
    if ~isempty(corn_arr)
        [corn_nrows, corn_ncols] = size(corn_arr);
        buff_b = getfield(m0, data_buff_name, {n_bottom,1});
        buff_b = buff_b{1};
        if ~isempty(buff_b)
            fprintf("Burning %s tile's right-side %s array bottom-right corner values into bottom-side array\n", working_tile_position, data_buff_name);
            fprintf("WARNING! %s tile will need to be remerged on the bottom side!\n", working_tile_position);
            buff_b(:, (end-corn_ncols+1):end) = buff_r((end-corn_nrows+1):end, :);
            setfield(m0, data_buff_name, {n_bottom,1}, {buff_b});
            m0.bottomNeedsRemerge = true;
        end
    end

elseif n0 == n_left
    buff_l = buff_arr;
    clear buff_arr;

    corn_arr = getfield(m0, data_corners_name, {n_top_left,1});
    corn_arr = corn_arr{1};
    if ~isempty(corn_arr)
        [corn_nrows, corn_ncols] = size(corn_arr);
        buff_t = getfield(m0, data_buff_name, {n_top,1});
        buff_t = buff_t{1};
        if ~isempty(buff_t)
            fprintf("Burning %s tile's left-side %s array top-left corner values into top-side array\n", working_tile_position, data_buff_name);
            fprintf("WARNING! %s tile will need to be remerged on the top side!\n", working_tile_position);
            buff_t(:, 1:corn_ncols) = buff_l(1:corn_nrows, :);
            setfield(m0, data_buff_name, {n_top,1}, {buff_t});
            m0.topNeedsRemerge = true;
        end
    end

    corn_arr = getfield(m0, data_corners_name, {n_bottom_left,1});
    corn_arr = corn_arr{1};
    if ~isempty(corn_arr)
        [corn_nrows, corn_ncols] = size(corn_arr);
        buff_b = getfield(m0, data_buff_name, {n_bottom,1});
        buff_b = buff_b{1};
        if ~isempty(buff_b)
            fprintf("Burning %s tile's left-side %s array bottom-left corner values into bottom-side array\n", working_tile_position, data_buff_name);
            fprintf("WARNING! %s tile will need to be remerged on the bottom side!\n", working_tile_position);
            buff_b(:, 1:corn_ncols) = buff_l((end-corn_nrows+1):end, :);
            setfield(m0, data_buff_name, {n_bottom,1}, {buff_b});
            m0.bottomNeedsRemerge = true;
        end
    end

elseif n0 == n_top
    buff_t = buff_arr;
    clear buff_arr;

    corn_arr = getfield(m0, data_corners_name, {n_top_right,1});
    corn_arr = corn_arr{1};
    if ~isempty(corn_arr)
        [corn_nrows, corn_ncols] = size(corn_arr);
        buff_r = getfield(m0, data_buff_name, {n_right,1});
        buff_r = buff_r{1};
        if ~isempty(buff_r)
            fprintf("Burning %s tile's top-side %s array top-right corner values into right-side array\n", working_tile_position, data_buff_name);
            fprintf("WARNING! %s tile will need to be remerged on the right side!\n", working_tile_position);
            buff_r(1:corn_nrows, :) = buff_t(:, (end-corn_ncols+1):end);
            setfield(m0, data_buff_name, {n_right,1}, {buff_r});
            m0.rightNeedsRemerge = true;
        end
    end

    corn_arr = getfield(m0, data_corners_name, {n_top_left,1});
    corn_arr = corn_arr{1};
    if ~isempty(corn_arr)
        [corn_nrows, corn_ncols] = size(corn_arr);
        buff_l = getfield(m0, data_buff_name, {n_left,1});
        buff_l = buff_l{1};
        if ~isempty(buff_l)
            fprintf("Burning %s tile's top-side %s array top-left corner values into left-side array\n", working_tile_position, data_buff_name);
            fprintf("WARNING! %s tile will need to be remerged on the left side!\n", working_tile_position);
            buff_l(1:corn_nrows, :) = buff_t(:, 1:corn_ncols);
            setfield(m0, data_buff_name, {n_left,1}, {buff_l});
            m0.leftNeedsRemerge = true;
        end
    end

elseif n0 == n_bottom
    buff_b = buff_arr;
    clear buff_arr;

    corn_arr = getfield(m0, data_corners_name, {n_bottom_right,1});
    corn_arr = corn_arr{1};
    if ~isempty(corn_arr)
        [corn_nrows, corn_ncols] = size(corn_arr);
        buff_r = getfield(m0, data_buff_name, {n_right,1});
        buff_r = buff_r{1};
        if ~isempty(buff_r)
            fprintf("Burning %s tile's bottom-side %s array bottom-right corner values into right-side array\n", working_tile_position, data_buff_name);
            fprintf("WARNING! %s tile will need to be remerged on the right side!\n", working_tile_position);
            buff_r((end-corn_nrows+1):end, :) = buff_b(:, (end-corn_ncols+1):end);
            setfield(m0, data_buff_name, {n_right,1}, {buff_r});
            m0.rightNeedsRemerge = true;
        end
    end

    corn_arr = getfield(m0, data_corners_name, {n_bottom_left,1});
    corn_arr = corn_arr{1};
    if ~isempty(corn_arr)
        [corn_nrows, corn_ncols] = size(corn_arr);
        buff_l = getfield(m0, data_buff_name, {n_left,1});
        buff_l = buff_l{1};
        if ~isempty(buff_l)
            fprintf("Burning %s tile's bottom-side %s array bottom-left corner values into left-side array\n", working_tile_position, data_buff_name);
            fprintf("WARNING! %s tile will need to be remerged on the left side!\n", working_tile_position);
            buff_l((end-corn_nrows+1):end, :) = buff_b(:, 1:corn_ncols);
            setfield(m0, data_buff_name, {n_left,1}, {buff_l});
            m0.leftNeedsRemerge = true;
        end
    end

end


function [r0_out,c0_out,r1_out,c1_out,failure,failure_skip_merge,failure_merge_anyway] = adjustFeatherZone(m0,m1,r0_in,c0_in,r1_in,c1_in,edge0,edge1,n0,n1,maxdist,mindist)

r0_out = r0_in;
c0_out = c0_in;
r1_out = r1_in;
c1_out = c1_in;
failure = false;
failure_skip_merge = false;
failure_merge_anyway = false;
unequal_water_issue = false;

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
if ismember('waterFillMask', varlist0) || ismember('waterFillMask', varlist1)
    data_array_names{end+1} = 'waterFillMask';
end

if ~all(ismember(data_array_names, varlist0))
    disp("ERROR: One or more expected data arrays do not exist in 1st tile struct");
    missing_arrays = setdiff(data_array_names, varlist0);
    data_array_names
    missing_arrays
    failure = true;
%    if length(missing_arrays) == 1 && strcmp(missing_arrays{1}, 'waterFillMask')
%        disp("Will skip merging from this error assuming that tile is missing waterFillMask for good reason")
%        failure_skip_merge = true;
%    end
    return;
end
if ~all(ismember(data_array_names, varlist1))
    disp("ERROR: One or more expected data arrays do not exist in 2nd tile struct");
    missing_arrays = setdiff(data_array_names, varlist1);
    data_array_names
    missing_arrays
    failure = true;
%    if length(missing_arrays) == 1 && strcmp(missing_arrays{1}, 'waterFillMask')
%        disp("Will skip merging from this error assuming that tile is missing waterFillMask for good reason")
%        failure_skip_merge = true;
%    end
    return;
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


% Shrink the feather zone as necessary to be within NoData inequality areas on edges
% This can let us avoid feathering bugs if there happens to be a stripe of NaN pixels on edges.
max_shrink_dist = overlap_halfwidth - mindist;
max_shrink_px = floor(max_shrink_dist / coord_res);
[buff_nrows, buff_ncols] = size(z0);
nan_equality_arr = (isnan(z0) == isnan(z1));
if ismember(edge0, {'right', 'left'})
    % Remove rows from the top and bottom of the NaN equality array
    % where a later "column" merge should fix a potential NaN stripe on the top or bottom edges.
    nan_equality_arr = nan_equality_arr((1+max_shrink_px):(buff_nrows-max_shrink_px), :);
    nan_equality_vec = all(nan_equality_arr, 1);
elseif ismember(edge0, {'top', 'bottom'})
    % Remove cols from the left and right of the NaN equality array
    % where a later "row" merge should fix a potential NaN stripe on the top or bottom edges.
    nan_equality_arr = nan_equality_arr(:, (1+max_shrink_px):(buff_ncols-max_shrink_px));
    nan_equality_vec = all(nan_equality_arr, 2);
end
clear nan_equality_arr;

pre_iter_shrink = shrink_px;
nan_check_vec = nan_equality_vec;
checking_all_nan_stripes = false;
checking_unequal_nan_stripes = false;
check_for_unequal_water_issue = false;
while true
    min_shrink_px = pre_iter_shrink;
    test_shrink_px = pre_iter_shrink;
    displayed_warning = false;
    iter_failure = false;

    if ismember(edge0, {'right', 'left'})
        v0 = [1+pre_iter_shrink, buff_ncols-pre_iter_shrink];
    elseif ismember(edge0, {'top', 'bottom'})
        v0 = [1+pre_iter_shrink, buff_nrows-pre_iter_shrink];
    end

    while v0(1) <= v0(2)
        % if all(nan_check_vec(v0(1):v0(2)))
        %     break;
        if nan_check_vec(v0(1)) == 1 && nan_check_vec(v0(2)) == 1
            ;
        else
            if ~displayed_warning
                if checking_unequal_nan_stripes
                    disp("Found unequal NaN stripes between tiles within overlap area");
                    disp("Will shrink feather zone one pixel at a time until out of unequal NaN stripes zone");
                elseif checking_all_nan_stripes
                    disp("Found all-NaN zone on at least one tile's edge");
                    disp("Will shrink feather zone one pixel at a time until out of all-NaN zone");
                else
                    disp("Found unequal NaN locations between tiles within overlap area");
                    disp("Will shrink feather zone one pixel at a time until out of unequal NaN zone");
                end
                displayed_warning = true;
            end
            min_shrink_px = test_shrink_px + 1;
            if min_shrink_px > max_shrink_px
                fprintf("Cannot continue to shrink feather zone more than %g pixels, reached minimum feather halfwidth of %gm\n", test_shrink_px, mindist);
                iter_failure = true;
                break;
            elseif (v0(1) + 1) >= (v0(2) - 1)
                fprintf("Cannot continue to shrink feather zone more than %g pixels, reached center of overlap area\n", test_shrink_px);
                iter_failure = true;
                break;
            end
        end
        test_shrink_px = test_shrink_px + 1;
        v0 = [v0(1)+1, v0(2)-1];
    end

    if iter_failure
        if ~checking_all_nan_stripes
            checking_all_nan_stripes = true;
            disp("Trying again, checking for all-NaN zone (stripe) on tile edges instead of unequal NaN zone")
            if ismember(edge0, {'right', 'left'})
                nan_check_vec = ~(all(isnan(z0), 1) | all(isnan(z1), 1));
            elseif ismember(edge0, {'top', 'bottom'})
                nan_check_vec = ~(all(isnan(z0), 2) | all(isnan(z1), 2));
            end
            continue;
        elseif ~checking_unequal_nan_stripes
            checking_unequal_nan_stripes = true;
            disp("Trying again, checking unequal all-NaN stripes on tile edges")
            if ismember(edge0, {'right', 'left'})
                nan_check_vec = (all(isnan(z0), 1) == all(isnan(z1), 1));
            elseif ismember(edge0, {'top', 'bottom'})
                nan_check_vec = (all(isnan(z0), 2) == all(isnan(z1), 2));
            end
            continue;
        else
            failure = true;
            check_for_unequal_water_issue = true;
        end
    end

    break;
end

if failure
    if check_for_unequal_water_issue && ismember('land', varlist0) && ismember('land', varlist1)
        water0 = ~m0.land(r0_in(1):r0_in(2),c0_in(1):c0_in(2));
        water1 = ~m1.land(r1_in(1):r1_in(2),c1_in(1):c1_in(2));

        if ismember(edge0, {'right', 'left'})
            if ~isequal(all(water0, 1), all(water1, 1))
                unequal_water_issue = true;
            end
        elseif ismember(edge0, {'top', 'bottom'})
            if ~isequal(all(water0, 2), all(water1, 2))
                unequal_water_issue = true;
            end
        end

        % I didn't want to have to resort to the following blanket check,
        % but unfortunately we have to.
        if ~unequal_water_issue
            if all(water0, 'all') && all(water1, 'all')
                unequal_water_issue = true;
            end
        end
    end
    if unequal_water_issue
        fprintf("Detected unequal land-water presence in tile overlap area")
        fprintf("Deeming this a REMA v2 era acceptable error, assuming tiles straddle land-water boundary\n");
        failure_merge_anyway = true;
    else
        fprintf("Deeming this an ArcticDEM v4.1 era acceptable error, assuming tiles straddle land-water boundary\n");
        failure_merge_anyway = true;
    end
    return;
end

shrink_px = min_shrink_px;
post_iter_shrink = shrink_px;

if pre_iter_shrink ~= post_iter_shrink && ~failure
    fprintf("Shrunk feather zone a total of %g pixels on each side of overlap area to be within unequal NaN zones on tile edges\n", shrink_px);

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
