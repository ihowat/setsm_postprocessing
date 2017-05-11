function mergeTileBuffer(f0,f1)
% mergeTileBuffer: merge the edge buffers two tiles with linear feathering
%
% mergeTileBuffer(f0,f1) where f0 and f1 are the file names or mat file
% handles of neighboring tiles.

%% first test if input args are either valid filenames or mat file handles
fprintf('Merging tile %s with %s\n', f0, f1)

if isstr(f0) % it's a string, might be a filename
    if exist(f0,'file') % yup its a file
        m0 = matfile(f0); % load it
    else
        error('File does not exist');
    end
elseif isvalid(f0) % not a string, is it a valid file handle?
    m0 = f0; % yes, so set m to f and get the filename for f
    f0 = m0.Properties.Source;
else error('input arg must be a filename or valid matfile handle')
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
else error('input arg must be a filename or valid matfile handle')
end

% make sure writeable
m0.Properties.Writable = true;
m1.Properties.Writable = true;

% get variable names in m files to check if already merged
m0vars = whos(m0); m0vars = {m0vars.name};
m1vars = whos(m1); m1vars = {m1vars.name};

%% crop tiles to buffer
c0 = m0.x >= min(m1.x) & m0.x <= max(m1.x);
r0 = m0.y >= min(m1.y) & m0.y <= max(m1.y);
c0 = [find(c0,1,'first'),find(c0,1,'last')];
r0 = [find(r0,1,'first'),find(r0,1,'last')];

if isempty(c0) || isempty(r0); error('no overlap betwene tiles'); end

c1 = m1.x >= min(m0.x) & m1.x <= max(m0.x);
r1 = m1.y >= min(m0.y) & m1.y <= max(m0.y);
c1 = [find(c1,1,'first'),find(c1,1,'last')];
r1 = [find(r1,1,'first'),find(r1,1,'last')];

%% Make weight array based on boundary being merged
info0=whos(m0,'z');
sz0=info0.size;

if c0(1) > 1; % merging on f0's right boundary
    
    % check for pre-merged
    m0check = any(strcmp(m0vars,'mergedRight'));
    m1check = any(strcmp(m1vars,'mergedLeft'));
    
    if m0check; warning('tile %s already merged, skipping\n', f0); end;
    if m1check; warning('tile %s already merged, skipping\n', f1); end;
    if m0check | m1check; return; end

    W0 = linspace(1,0,diff(c0)+1);
    W0 = repmat(W0,diff(r0)+1,1);
    
    m0.mergedRight=true;
    m1.mergedLeft=true;
    
elseif c0(2) < sz0(2) % merging on f0's left boundary
    
    % check for pre-merged
    m0check = any(strcmp(m0vars,'mergedLeft'));
    m1check = any(strcmp(m1vars,'mergedRight'));
    
    if m0check; warning('tile %s already merged, skipping\n', f0); end;
    if m1check; warning('tile %s already merged, skipping\n', f1); end;
    if m0check | m1check; return; end
    
    W0 = linspace(0,1,diff(c0)+1);
    W0 = repmat(W0,diff(r0)+1,1);
    
    m0.mergedLeft=true;
    m1.mergedRight=true;
    
elseif r0(2) < sz0(2) % merging on f0's top boundary
    
    % check for pre-merged
    m0check = any(strcmp(m0vars,'mergedTop'));
    m1check = any(strcmp(m1vars,'mergedBottom'));
    
    if m0check; warning('tile %s already merged, skipping\n', f0); end;
    if m1check; warning('tile %s already merged, skipping\n', f1); end;
    if m0check | m1check; return; end
    
    W0 = linspace(0,1,diff(r0)+1);
    W0 = repmat(W0(:),1,diff(c0)+1);
    
    m0.mergedTop=true;
    m1.mergedBottom=true;

elseif r0(1) > 1 % merging on f0's bottom boundary
    
    % check for pre-merged
    m0check = any(strcmp(m0vars,'mergedBottom'));
    m1check = any(strcmp(m1vars,'mergedTop'));
    
    if m0check; warning('tile %s already merged, skipping\n', f0); end;
    if m1check; warning('tile %s already merged, skipping\n', f1); end;
    if m0check | m1check; return; end
    
    W0 = linspace(1,0,diff(r0)+1);
    W0 = repmat(W0(:),1,diff(c0)+1);
    
    m0.mergedBottom=true;
    m1.mergedTop=true;
    
else
    
    error('buffer region is same size as full grid')
    
end

%% merge z
z0= m0.z(r0(1):r0(2),c0(1):c0(2));
z1= m1.z(r1(1):r1(2),c1(1):c1(2));

z=(z0.*W0)+(z1.*(1-W0));

n=isnan(z0) & ~isnan(z1);
z(n)=z1(n);
n=~isnan(z0) & isnan(z1);
z(n)=z0(n);

m0.z(r0(1):r0(2),c0(1):c0(2))=z;
m1.z(r1(1):r1(2),c1(1):c1(2))=z;
clear z0 z1 z;

%% merge mt
mt0= m0.mt(r0(1):r0(2),c0(1):c0(2));
mt1= m1.mt(r1(1):r1(2),c1(1):c1(2));

mt = mt0 | mt1;

m0.mt(r0(1):r0(2),c0(1):c0(2))=mt;
m1.mt(r1(1):r1(2),c1(1):c1(2))=mt;
clear mt0 mt1 mt;

%% merge or
or0= m0.or(r0(1):r0(2),c0(1):c0(2));
or1= m1.or(r1(1):r1(2),c1(1):c1(2));

or0=single(or0);
or1=single(or1);

or=(or0.*W0)+(or1.*(1-W0));

n= or0 == 0 & or1 ~= 0;
or(n)=or1(n);
n= or0 ~= 0 & or1 == 0;
or(n)=or0(n);

or=int16(or);

m0.or(r0(1):r0(2),c0(1):c0(2))=or;
m1.or(r1(1):r1(2),c1(1):c1(2))=or;
clear or0 or1 or;

%% merge dy
dy0= m0.dy(r0(1):r0(2),c0(1):c0(2));
dy1= m1.dy(r1(1):r1(2),c1(1):c1(2));

dy0=single(dy0);
dy1=single(dy1);

dy=(dy0.*W0)+(dy1.*(1-W0));

n= dy0 == 0 & dy1 ~= 0;
dy(n)=dy1(n);
n= dy0 ~= 0 & dy1 == 0;
dy(n)=dy0(n);

dy=int16(dy);

m0.dy(r0(1):r0(2),c0(1):c0(2))=dy;
m1.dy(r1(1):r1(2),c1(1):c1(2))=dy;
clear dy0 dy1 dy;
