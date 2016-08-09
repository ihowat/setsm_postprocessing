function [outFlag,outname]=alignTile(f)

% set out flag
outFlag=false; 
outname=[];

%% first test if input args are either valid filenames or mat file handles
nfiles=length(f);
for i=1:nfiles;
    
    if ischar(f{i}) % it's a string, might be a filename
        if exist(f{i},'file') % yup its a file
            m = matfile(f{i}); % load it
        else
            error('File does not exist');
        end
    elseif isvalid(f{i}) % not a string, is it a valid file handle?
        m = f{i}; % yes, so set m to f and get the filename for f
        f{i} = m.Properties.Source;
    else error('input arg must be a filename or valid matfile handle')
    end
    
    eval(['f',num2str(i-1),'=f{i};']);
    eval(['m',num2str(i-1),'=m;']);
    
    clear m
end

info0=whos(m0,'z');
sz0=info0.size;

%% Coregister to each neighbor
for i=2:nfiles;

    eval(['m=m',num2str(i-1),';']);

    % crop reference to buffer
    c0 = m0.x >= min(m.x) & m0.x <= max(m.x);
    r0 = m0.y >= min(m.y) & m0.y <= max(m.y);
    c0 = [find(c0,1,'first'),find(c0,1,'last')];
    r0 = [find(r0,1,'first'),find(r0,1,'last')];

    if isempty(c0) || isempty(r0); error('no overlap between tiles'); end

    % crop floating tile to buffer
    c1 = m.x >= min(m0.x) & m.x <= max(m0.x);
    r1 = m.y >= min(m0.y) & m.y <= max(m0.y);
    c1 = [find(c1,1,'first'),find(c1,1,'last')];
    r1 = [find(r1,1,'first'),find(r1,1,'last')];

    % Coregistration cluster array
    C0=m0.C(r0(1):r0(2),c0(1):c0(2));
    
    coregClusters=max(C0(:))-1;
    
    if coregClusters == -1; fprintf('no overlapping data in %s\n',f{i}); continue; end
    
    if coregClusters == 0; error('Grid appears to be already registered'); end
    
    %% pull buffer data and align
    x0 = m0.x(1,c0(1):c0(2));
    y0 = m0.y(r0(1):r0(2),1);
    z0 = m0.z(r0(1):r0(2),c0(1):c0(2));
    mt0= m0.mt(r0(1):r0(2),c0(1):c0(2));
    x1 =  m.x(1,c1(1):c1(2));
    y1 =  m.y(r1(1):r1(2),1);
    z1 =  m.z(r1(1):r1(2),c1(1):c1(2));
    mt1=  m.mt(r0(1):r0(2),c0(1):c0(2));
       
    % cluster coregistraton loop
    for j=1:coregClusters
    
        fprintf('Registering coregistered cluster %d of %d\n',j,coregClusters);
    
        % make a mask of this cluster
        mt0c= mt0 & (C0 == j+1);
     
        % co-register floating tile to reference tile
        [~,dtrans{j,i-1},rmse{j,i-1}] = coregisterdems(x0, y0, z0, x1, y1, z1, mt0c, mt1);

    end
    
    clear m r0 c0 r1 c1 C0 x0 y0 z0 mt0 x1 y1 z1 mt1 mt0c
 
end
   
% take average dtrans for each cluster
for i=1:size(dtrans,1)
    dtrans{i,1} = nanmean([dtrans{i,:}],2);
    rmse{i,1} = nanmean([rmse{i,:}],2);
end

dtrans=dtrans(:,1);
rmse=rmse(:,1);

% reorganize dtrans and rmse for consistency with gcp reg files
dtrans=dtrans';
dtrans = [dtrans{:}];

rmse=rmse';
rmse = [rmse{:}];

if ~any(~isnan(rmse)); 
    fprintf('coregistration failure\n'); 
    return
end

%% create output file
outname=strrep(m0.Properties.Source,'dem.mat','_reg_dem.mat');

clear m1
m1 = matfile(outname,'Writable',true);

m1.dtrans = dtrans;
m1.rmse   = rmse;

%% Regridding 

% get cluster array
C=m0.C;

% Apply registration to z grid variable
z= nan(sz0); % initialize output
for i=1:size(dtrans,2)
    
    if any(isnan(dtrans(:,i))); continue; end

    % make a mask of this cluster
    N = C == i+1;
 
    % apply registration
    ztemp = applyRegistration(dtrans(:,i),m0,N);
    clear N
    
    n=isnan(z) & ~isnan(ztemp);
    z(n) = ztemp(n);
    
    clear ztemp n
    
end

m1.z=z;

clear z

% Apply registration to mt grid variable

% cluster coregistraton loop
mt = false(sz0); % initialize output
for i=1:size(dtrans,2)
    
    if any(isnan(m.dtrans(:,i))); continue; end

    % make a mask of this cluster
    N = C == i+1;
 
    % Apply registration
    mttemp = applyRegistration(m.dtrans(:,i),m0,N,'mt');
    
    clear N
    
    % add in any matches
    mt = mt | mttemp;
    
    clear mttemp
    
end

m1.mt=mt;

clear mt

%% Apply registration to or grid variable

% cluster coregistraton loop
or = zeros(sz0,'int16'); % initialize output
for i=1:size(dtrans,2)
    
    if any(isnan(m.dtrans(:,i))); continue; end
    
    % make a mask of this cluster
    N = C == i+1;
 
    % Apply registration
    ortemp = applyRegistration(m.dtrans(:,i),m0,N,'or');
    
    clear N
    
    n= or==0 & ortemp ~= 0;
    or(n) = ortemp(n);
    
    clear ortemp n
    
end

m1.or=or;

clear or

%% Apply registration to dy grid variable

% cluster coregistraton loop
dy= zeros(sz0,'int16'); % initialize output
for i=1:size(dtrans,2)
    
    if any(isnan(m.dtrans(:,i))); continue; end
    
    % make a mask of this cluster
    N = C == i+1;
 
    % Apply registration
    dytemp = applyRegistration(m.dtrans(:,i),m0,N,'dy');
    
    clear N
    
    n= dy==0 & dytemp ~= 0;
    dy(n) = dytemp(n);
    
    clear dytemp n
    
end

m1.dy=dy;

clear dy

tileRegMeta(m1,{'Neighbor Align'})

outFlag=true;
    









