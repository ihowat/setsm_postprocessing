function mergeTileBuffer(f0,f1,varargin)
% mergeTileBuffer: merge the edge buffers two tiles with linear feathering
%
% mergeTileBuffer(f0,f1) where f0 and f1 are the file names or mat file
% handles of neighboring tiles.

% top = 1, Bottom = 2, left = 3, right = 4;

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
if nargin==3
    overwriteBuffer=true;
end

% make sure writeable
m0.Properties.Writable = true;
m1.Properties.Writable = true;

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
varlist0 = who(m0);
varlist1 = who(m1);

if c0(1) > 1 % merging on f0's right boundary
    
    if  any(strcmp(varlist0,'mergedRight')) && any(strcmp(varlist1,'mergedLeft')) && ~overwriteBuffer
        if m0.mergedRight && m1.mergedLeft
            disp('already merged');
            return;
        end
    end
      
    W0 = linspace(1,0,diff(c0)+1);
    W0 = repmat(W0,diff(r0)+1,1);
    
    m0.mergedRight=true;
    m1.mergedLeft=true;
    
    n0 = 4;
    n1 = 3;

elseif c0(2) < sz0(2) % merging on f0's left boundary
    
    if  any(strcmp(varlist0,'mergedLeft')) && any(strcmp(varlist1,'mergedRight')) && ~overwriteBuffer
        if m0.mergedLeft && m1.mergedRight
            disp('already merged');
            return;
        end
    end
    
    
    W0 = linspace(0,1,diff(c0)+1);
    W0 = repmat(W0,diff(r0)+1,1);
    
    m0.mergedLeft=true;
    m1.mergedRight=true;
    
    n0 = 3;
    n1 = 4;

elseif r0(2) < sz0(1) % merging on f0's top boundary
    
    if  any(strcmp(varlist0,'mergedTop')) && any(strcmp(varlist1,'mergedBottom')) && ~overwriteBuffer
        if m0.mergedTop && m1.mergedBottom
            disp('already merged');
            return;
        end
    end
    
    W0 = linspace(0,1,diff(r0)+1);
    W0 = repmat(W0(:),1,diff(c0)+1);
    
    m0.mergedTop=true;
    m1.mergedBottom=true;
    
    n0 =1;
    n1 =2;
    
elseif r0(1) > 1 % merging on f0's bottom boundary
    
    if  any(strcmp(varlist0,'mergedBottom')) && any(strcmp(varlist1,'mergedTop')) && ~overwriteBuffer
        if m0.mergedBottom && m1.mergedTop
            disp('already merged');
            return;
        end
    end
    
    W0 = linspace(1,0,diff(r0)+1);
    W0 = repmat(W0(:),1,diff(c0)+1);
    
    m0.mergedBottom=true;
    m1.mergedTop=true;
    
    n0 =2;
    n1 = 1;
    
else
    error('buffer region is same size as full grid')
end

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

N1=m1.Nbuff(n0,1);
N1=N1{1};

if isempty(N1)
        N1= m0.N(r0(1):r0(2),c0(1):c0(2));
        m1.Nbuff(n0,1) = {N1};
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

Nmt1=m1.Nmtbuff(n0,1);
Nmt1=Nmt1{1};

if isempty(Nmt1)
        Nmt1= m0.Nmt(r0(1):r0(2),c0(1):c0(2));
        m1.Nmtbuff(n0,1) = {Nmt1};
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

tmin1=m1.tminbuff(n0,1);
tmin1=tmin1{1};

if isempty(tmin1)
        tmin1= m0.tmin(r0(1):r0(2),c0(1):c0(2));
        m1.tminbuff(n0,1) = {tmin1};
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

tmax1=m1.tmaxbuff(n0,1);
tmax1=tmax1{1};

if isempty(tmax1)
        tmax1= m0.tmax(r0(1):r0(2),c0(1):c0(2));
        m1.tmaxbuff(n0,1) = {tmax1};
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
