function z = applyRegistration(dtrans,m,N,varargin)
% ApplyRegistration applies 3-D transformation to grid. V4 compatible
%   z = applyRegistration(dtrans,m,z,N) dtrans is [x,y,z] offsets to be
%   *subtracted* . m is matfile object. N is coregistration mask.
%   z = applyRegistration(...,'gridVar','varname')
%   z = applyRegistration(...,'subsetSize',length)

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

intsCol=1:subsetSize:sz(2); 
intsCol=[intsCol(1:end-1)',intsCol(2:end)'];
intsCol(end)=sz(2);
intsCol(1:end-1,2)=intsCol(1:end-1,2)-1;

intsRow=1:subsetSize:sz(1); 
intsRow=[intsRow(1:end-1)',intsRow(2:end)'];
intsRow(end)=sz(1);
intsRow(1:end-1,2)=intsRow(1:end-1,2)-1;

% need to add pixel buffer to int to account for shifts, scaling with dtrans
buff=max(ceil(abs(dtrans(2:3)./res)));
buffsCol=[[0 buff];repmat([-buff buff],size(intsCol,1)-2,1);[-buff 0]];
buffsRow=[[0 buff];repmat([-buff buff],size(intsRow,1)-2,1);[-buff 0]];

%% nested interpolation loops with interpolation dependent on grid variable

% initialize output
switch gridVar
    case {'z','z_mad'}
        z = nan(sz,'single');
    case {'tmin','tmax'}
        z= zeros(sz,'uint16');
    case {'N','Nmt'}
        z = zeros(sz,'uint8');  
end

% subset interpolation loop
for i=1:length(intsRow)
    
    % make subset row vector with buffer
    intRowBuff=intsRow(i,1)+buffsRow(i,1):intsRow(i,2)+buffsRow(i,2);
    
    for j=1:length(intsCol)
        
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
                    x(intCol),y(intRow),'*linear');

            case 'z_mad' % error grid
                
                % apply only x and y offsets and interpolate
                z(intRow,intCol) =interp2(x1,y1, z0,...
                    x(intCol),y(intRow),'*linear');
                
            case {'N','Nmt'} % stack and matchtag numbers
                
                z0(isnan(z0))=0;
                
                % apply horizontal offsets and interpolate
                z0=interp2(x1,y1, z0, x(intCol),y(intRow),'*nearest');
                
                % convert back to uint8 and insert into array
                z0(isnan(z0)) = 0; 
                z(intRow,intCol)=uint8(z0);
                
                clear z0
                
            case 'or'
                
                % apply horizontal offsets and interpolate
                z0=interp2(x1,y1, z0, x(intCol),y(intRow),'*cubic');
                
                % convert back to int16 and insert into array
                z0(isnan(z0)) = 0; 
                z(intRow,intCol)=int16(z0);
                
            case {'tmin','tmax'}
                
                % apply horizontal offsets and interpolate
                z0=interp2(x1,y1, z0, x(intCol),y(intRow),'*linear');
                
                z0(isnan(z0)) = 0; % convert back to uint8
                z(intRow,intCol)=uint16(z0);
                
            otherwise
                error('grid variable input argument not recognized, must be ''z'',''z_mad'',''tmin'',''tmax'',''N'' or ''Nmt''')
                  
        end
        
        clear z0
        
    end
end




