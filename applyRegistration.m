function applyRegistration(demMatFile,registrationFile,varargin)

% File to register

if ~exist(demMatFile,'file')
    warning('%s doesnt exist, skipping',demMatFile)
    return
end

% make out mat name to write to
outName=strrep(demMatFile,'.mat','_reg.mat');

if exist(outName,'file') && ~any(strcmpi(varargin,'overwrite'))
    fprintf('%s exists, skipping\n',outName)
    return
end

% get the region and tile names from the mat file name
S=strsplit(demMatFile,'/');

tileName=S{end}(1:5);

% load demMatFile into a matfile object
m = matfile(demMatFile);

%find this tile
if ~isempty(registrationFile)
    regData=load(registrationFile);
    n = find(contains(regData.tileNames,tileName));
    
    if isempty(n)
        warning('could not find registration data for %s, skipping\n',tileName)
        return
    end
    
    % extract data for this tile from registration data
    regData = rmfield(regData,{'tileFiles','tileNames'});
    regData.reg = structfun( @(x) x(n,:), regData.reg,'uniformoutput',0);
    regData.unreg = structfun( @(x) x(n,:), regData.unreg,'uniformoutput',0);
    
    m.Properties.Writable = true;
    m.reg=regData.reg;
    m.unreg=regData.unreg;
    
else
    if ~any(strcmp(fields(m),'reg'))
        warning('no reg field in %s\n',demMatFile)
        return
    end
    
    regData.reg=m.reg;
    regData.unreg=m.unreg;
    
end

if isnan(regData.reg.p(1))
    warning('NaN offset for %s, skipping',tileName)
    return
end

% get grid size
[nrows,ncols] = size(m,'z');

% make coregistration cluster mask - not currently used
C = true(nrows,ncols);

% if only vertical registration only z is changed.
if regData.reg.p(2) == 0 && regData.reg.p(3) == 0
    
    % copy the matfile to the outname
    eval(['!cp ',demMatFile,' ',outName]);
    
    % subtract the vertical offset from the dem
    m1=matfile(outName);
    m1.Properties.Writable = true;
    m1.z = m1.z - regData.reg.p(1);
    
    % else do interpolation on each var
else
    % load x and y coordintes
    x = m.x;
    y = m.y;
    
    reg=regData.reg;
    unreg=regData.unreg;
    
    % save to matfile
    save(outName,'reg','unreg','x','y','-v7.3');
    
    m1=matfile(outName);
    m1.Properties.Writable = true;
    
    % register z and add to mat file
    m1.z = applyRegistration2Var(regData.reg.p,m,C,'gridvar','z','subsetSize',5000);
    
    
    % register z_mad and add to mat file
    m1.z_mad =  applyRegistration2Var(regData.reg.p,m,C,'gridvar','z_mad','subsetSize',5000);
    
    % register N and add to mat file
    m1.N =  applyRegistration2Var(regData.reg.p,m,C,'gridvar','N','subsetSize',5000);
    
    % register Nmt and add to mat file
    m1.Nmt =  applyRegistration2Var(regData.reg.p,m,C,'gridvar','Nmt','subsetSize',5000);
    
    % register tmin and add to mat file
    m1.tmin =  applyRegistration2Var(regData.reg.p,m,C,'gridvar','tmin','subsetSize',5000);
    
    % register tmax and add to mat file
    m1.tmax = applyRegistration2Var(regData.reg.p,m,C,'gridvar','tmax','subsetSize',5000);
    
    % register land mask to mat file if exists
    if any(strcmp(fields(m),'land'))
        m1.land = applyRegistration2Var(regData.reg.p,m,C,'gridvar','land','subsetSize',5000);
    end
end

function z = applyRegistration2Var(dtrans,m,N,varargin)
% ApplyRegistration2Var applies 3-D transformation to grid. V4 compatible
%   z = applyRegistration2Var(dtrans,m,z,N) dtrans is [x,y,z] offsets to be
%   *subtracted* . m is matfile object. N is coregistration mask.
%   z = applyRegistration2Var(...,'gridVar','varname')
%   z = applyRegistration2Var(...,'subsetSize',length)

%check dtrans
if numel(dtrans) ~= 3
    error('dtrans must be a 3 element vector')
end
if any(isnan(dtrans))
    error('dtrans contains NaNs')
end

% parse/error check gridVar argument
gridVar='z'; %default
subsetSize=5000; % dimension of interpoaltion subset(pixels=subsetSize^2)


if ~isempty(varargin)
  
    n = find(strcmpi(varargin,'gridvar'));
    if ~isempty(n)
        gridVar=varargin{n+1};
    end
    
    n = find(strcmpi(varargin,'subsetsize'));
    if ~isempty(n)
        subsetSize=varargin{n+1};
    end
end

% extract grid vectors for speed rather in pulling on each iter
x=m.x;
y=m.y;
y = y(:);
res = x(2)-x(1);

% to speed up/converve memory, we'll perform gridded interpolation over
% 1000x1000 pixel subsets, overlapping by one pixel. 
% Make a list of min max indices.
sz = [length(y) length(x)];

if subsetSize < 0.75*sz(2)
    
    intsCol=1:subsetSize:sz(2);
    intsCol=[intsCol(1:end-1)',intsCol(2:end)'];
    intsCol(end)=sz(2);
    intsCol(1:end-1,2)=intsCol(1:end-1,2)-1;
    
else
    intsCol =[1 sz(2)];
end


if subsetSize < .75*sz(1)
    
    intsRow=1:subsetSize:sz(1);
    intsRow=[intsRow(1:end-1)',intsRow(2:end)'];
    intsRow(end)=sz(1);
    intsRow(1:end-1,2)=intsRow(1:end-1,2)-1;
    
else
    intsRow =[1 sz(1)];
end

% need to add pixel buffer to int to account for shifts, scaling with dtrans
buff=max(ceil(abs(dtrans(2:3)./res)));
if size(intsCol,1) == 1
    buffsCol = [0,0];
else
    buffsCol=[[0 buff];repmat([-buff buff],size(intsCol,1)-2,1);[-buff 0]];
end

if size(intsRow,1) == 1
    buffsRow=[0,0];
else
    buffsRow=[[0 buff];repmat([-buff buff],size(intsRow,1)-2,1);[-buff 0]];
end

%% nested interpolation loops with interpolation dependent on grid variable

% initialize output
switch gridVar
    case {'z','z_mad'}
        z = nan(sz,'single');
    case {'tmin','tmax'}
        z= zeros(sz,'uint16');
    case {'N','Nmt','land'}
        z = zeros(sz,'uint8');  
end

% subset interpolation loop
for i=1:size(intsRow,1)
    
    % make subset row vector with buffer
    intRowBuff=intsRow(i,1)+buffsRow(i,1):intsRow(i,2)+buffsRow(i,2);
    
    for j=1:size(intsCol,1)
        
        % make subset column vector with buffer
        intColBuff=intsCol(j,1)+buffsCol(j,1):intsCol(j,2)+buffsCol(j,2);
        
        % extract subset from unregistered tile with buffer and add offset
        eval(['z0 = m.',gridVar,'(intRowBuff,intColBuff);']);
        
        % convert to single for interpolation
        z0 = single(z0);
        
        % apply coregistration cluster mask
        z0(~N(intRowBuff,intColBuff)) = NaN;
        
        % check for empty subset to skip
        if ~any(~isnan(z0(:)) & z0(:) ~= 0); continue; end
        
        % make unbuffered row and column vectors to insert into new
        % array
        intRow=intsRow(i,1):intsRow(i,2);
        intCol=intsCol(j,1):intsCol(j,2);
        
        x1 = x(intColBuff)-dtrans(2);
        y1 = y(intRowBuff)-dtrans(3);
        
        % interpolate depending on variable
        switch gridVar
            
            case 'z'
                
                % apply xy,z offsets and interpolate
                z(intRow,intCol) =interp2(x1,y1, z0-dtrans(1),...
                    x(intCol),y(intRow),'*linear',NaN);

            case 'z_mad' % error grid
                
                % apply only x and y offsets and interpolate
                z(intRow,intCol) =interp2(x1,y1, z0,...
                    x(intCol),y(intRow),'*linear',NaN);
                
            case {'N','Nmt','land'} % stack, matchtag numbers and land mask
                
                z0(isnan(z0))=0;
                
                % apply horizontal offsets and interpolate
                z0=interp2(x1,y1, z0, x(intCol),y(intRow),'*nearest',0);
                
                % convert back to uint8 and insert into array
                z0(isnan(z0)) = 0; 
                z(intRow,intCol)=uint8(z0);
                
                clear z0
                
            case 'or'
                
                % apply horizontal offsets and interpolate
                z0=interp2(x1,y1, z0, x(intCol),y(intRow),'*cubic',0);
                
                % convert back to int16 and insert into array
                z0(isnan(z0)) = 0; 
                z(intRow,intCol)=int16(z0);
                
            case {'tmin','tmax'}
                
                % apply horizontal offsets and interpolate
                z0=interp2(x1,y1, z0, x(intCol),y(intRow),'*linear',0);
                
                z0(isnan(z0)) = 0; % convert back to uint8
                z(intRow,intCol)=uint16(z0);
                
            otherwise
                error('grid variable input argument not recognized, must be ''z'',''z_mad'',''tmin'',''tmax'',''N'',''land'' or ''Nmt''')
                  
        end
        
        clear z0
        
    end
end




