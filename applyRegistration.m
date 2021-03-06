function z = applyRegistration(dtrans,m,N,varargin)
%   z = applyRegistration(dtrans,m,z,N)
%   z = applyRegistration(dtrans,m,z,gridVar)

%check dtrans
if numel(dtrans) ~= 3
    error('dtrans must be a 3 element vector')
end
if any(isnan(dtrans))
    error('dtrans contains NaNs')
end

% parse/error check gridVar argument
gridVar='z'; %default
if ~isempty(varargin)
    gridVar=varargin{1};
    if ~ischar(gridVar)
        error('grid variable argument must be a string');
    end
end

% extract grid vectors for speed rather in pulling on each iter
x=m.x;
y=m.y;
res = x(2)-x(1);

% to speed up/converve memory, we'll perform gridded interpolation over
% 2509x2509pixel subsets, giving 20x20 subsets = 400 total. Each subset
% overlaps by one pixel for merging.
varinfo=whos(m,gridVar);
int=1:2510:varinfo.size(2);
int=[int(1:end-1)',int(2:end)'];
int(1:end-1,2)=int(1:end-1,2)-1;

% need to add pixel buffer to int to account for shifts, scaling with dtrans
buff=max(ceil(abs(dtrans(2:3)./res)));
buff=[[0 buff];repmat([-buff buff],size(int,1)-2,1);[-buff 0]];

%% nested interpolation loops with interpolation dependent on grid variable

% initialize output
switch gridVar
    case 'z'
        z = nan(varinfo.size,'single');
    case 'mt'
        z = false(varinfo.size);
    case {'or','dy'}
        z= zeros(varinfo.size,'int16');
end

% subset interpolation loop
for i=1:length(int)
    
    % make subset row vector with buffer
    intRowBuff=int(i,1)+buff(i,1):int(i,2)+buff(i,2);
    
    for j=1:length(int)
        
        % make subset column vector with buffer
        intColBuff=int(j,1)+buff(j,1):int(j,2)+buff(j,2);
        
        % extract subset from unregistered tile with buffer and add offset
        z0 = single(m.z(intRowBuff,intColBuff));
        
        % apply coregistration cluster mask
        z0(~N(intRowBuff,intColBuff)) = NaN;
        
        % check for empty subset to skip
        if ~any(~isnan(z0(:)) & z0(:) ~= 0); continue; end
        
        % make unbuffered row and column vectors to insert into new
        % array
        intRow=int(i,1):int(i,2);
        intCol=int(j,1):int(j,2);
        
        % interpolate depending on variable
        switch gridVar
            
            case 'z'
                
                % apply xy,z offsets and interpolate
                z(intRow,intCol)=interp2(x(intColBuff)+dtrans(2),...
                    y(intRowBuff)+dtrans(3), z0 + dtrans(1),...
                    x(intCol),y(intRow),'*linear');
                
            case 'mt'
                
                z0(isnan(z0))=0;
                
                % apply horizontal offsets and interpolate
                z0 = interp2(x(intColBuff)+dtrans(2),...
                    y(intRowBuff)+dtrans(3), z0,...
                    x(intCol),y(intRow),'*nearest');
                
                z0(isnan(z0)) = 0; % convert back to uint8
                
                z(intRow,intCol)=logical(z0);
                
                clear z0
                
            case 'or'
                
                % apply horizontal offsets and interpolate
                z0 = interp2(x(intColBuff)+dtrans(2),...
                    y(intRowBuff)+dtrans(3), z0,...
                    x(intCol),y(intRow),'*cubic');
                
                z0(isnan(z0)) = 0; % convert back to uint8
                z(intRow,intCol)=int16(z0);
                
            case 'dy'
                
                % apply horizontal offsets and interpolate
                z0 = interp2(x(intColBuff)+dtrans(2),...
                    y(intRowBuff)+dtrans(3), z0,...
                    x(intCol),y(intRow),'*linear');
                
                z0(isnan(z0)) = 0; % convert back to uint8
                z(intRow,intCol)=int16(z0);
                
            otherwise
                error('grid variable input argument not recognized, must be ''z'',''mt'',''or'' or ''dy''')
                  
        end
        
        clear z0
        
    end
end




