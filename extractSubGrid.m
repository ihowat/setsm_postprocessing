function [x,y,z,missingFlag,m] =extractSubGrid(fileNames,subx0,subx1,suby0,suby1,res,varargin)
% extractSubGrid extract dem strips data to 3D "stack" of subsets
%
% [x,y,z] =extractSubGrid(fileNames,subx0,subx1,suby0,suby1,res) extracts
% the z values from the strips listed in cellstr fileNames within the
% rectangular cartesian grid defined by ranges subx,suby, at the grid
% spacing res. Each fileName(k) will be added as  z(:,:,k);
%
% [...,missingFlag,m] =extractSubGrid(...) also resturns a logical vector
% where missingFlag(k) = true if z(:,:,k) has no non-nan values, and the
% stack of matchtag rasters.

%% parameters/defaults
verbose = false;
applyBitmask=true;

% parse optional arguments
n = find(strcmpi('verbose',varargin));
if ~isempty(n)
    verbose = varargin{n+1};
    if ~islogical(verbose)
        error('Verbose option must be logical')
    end
end

n = find(strcmpi('applyBitmask',varargin));
if ~isempty(n)
    applyBitmask = varargin{n+1};
    if ~islogical(applyBitmask)
        error('applyBitmask option must be logical')
    end
end

%% make subset grid
x=subx0:res:subx1;
y=suby1:-res:suby0;
y=y(:);
z=nan(length(y),length(x),length(fileNames),'single');

m=[];
% create matchtag output if specified
if nargout == 5
    m = false(size(z));
end

% subset extraction loop
missingFlag=false(size(fileNames));
for i=1:length(fileNames)
    
    if verbose
        fprintf('extracting subset %d of %d\n',i,length(fileNames))
    end
    
    % first check map grid to ensure extents cover destination grid
    zsub=readGeotiff(fileNames{i},...
        'map_subset',[x(1) x(end) y(end) y(1)],'mapinfoonly');
    
    % add one pixel to reange if not
    buff = [-res res -res res].* ...
        [zsub.x(1) > x(1) zsub.x(end) < x(end) ...
        zsub.y(end)> y(end) zsub.y(1) < y(1) ];
    
    % load dem with with buffered range
    zsub=readGeotiff(fileNames{i},...
        'map_subset',[x(1) x(end) y(end) y(1)] +buff);
    
    % skip if too small a subset or no data
    if ~any(zsub.z(:) ~= -9999) || any(size(zsub.z) < 4)
        missingFlag(i)=true;
        continue
    end
    
    % convert nodata to NaN for interpolation
    zsub.z(zsub.z == -9999) = NaN;
    
    %% Apply bitmask if reqested
    if  applyBitmask
        % check for bitmask and apply
        bitMaskFile=strrep(fileNames{i},'_dem','_bitmask');
        if exist(bitMaskFile,'file')
            
            if verbose
                fprintf('applying bitmask file: %s\n',bitMaskFile)
            end
            
            masksub=readGeotiff(bitMaskFile,...
                'map_subset',[x(1) x(end) y(end) y(1)] +buff);
            
            % skip layer if all data missing (in bitmask 0 = good data)
            if ~any(masksub.z(:) == 0)
                missingFlag(i)=true;
                continue
            end
            
            %         % convert bitmask to edge mask
            %         masksub.z = masksub.z ~= 1 & masksub.z ~= 3 & masksub.z ~= 5 ...
            %             & masksub.z ~= 7;
            % set to good data masks
            
            % convert all good data (0) to true
            masksub.z = masksub.z==0;
            
            % if all good, don't bother applying
            if any(~masksub.z(:))
                
                % interpolate if not same size/greid
                if length(zsub.x) ~= length(masksub.x) || ...
                        length(zsub.y) ~= length(masksub.y) || ...
                        zsub.x(1) ~= masksub.x(1) || ...
                        zsub.y(1) ~= masksub.y(1)
                    
                    a = interp2(masksub.x,masksub.y(:),single(masksub.z),...
                        zsub.x,zsub.y(:),'*nearest');
                    a(isnan(a))= 0;
                    masksub.z = logical(a);
                    clear a
                end
                
                %apply mask
                zsub.z(~masksub.z) = NaN;
            end
        else
            error('bitmask file: %s not found',bitMaskFile)
        end
    end
    
    %% interpolate to grid if needed
    if length(zsub.x) ~= length(x) || length(zsub.y) ~= length(y) || ...
            zsub.x(1) ~= x(1) || ...
            zsub.y(1) ~= y(1)
        z(:,:,i) = interp2(zsub.x,zsub.y(:),zsub.z,x,y(:),'*linear');
    else
        z(:,:,i) = zsub.z;
    end

    %% return mt's if requested
    if ~isempty(m)
        
        if res == 10
            matchtagFileName = strrep(fileNames{i},'dem_10m.tif','matchtag.tif');
        elseif res == 2
            matchtagFileName = strrep(fileNames{i},'dem.tif','matchtag.tif');
        end
        
        % check if filename exists, set to true if not 
        if ~exist(matchtagFileName,'file')
            fprintf('matchtag file doesnt exist, setting this layer to true\n')
            m(:,:,i) = true;
            continue
        end
        
        % check coverage
        msub=readGeotiff(matchtagFileName,'mapinfoonly',...
            'map_subset',[x(1) x(end) y(end) y(1)]);
        
        % add one pixel to reange if not covering full range
        buff = [-res res -res res].* ...
            [msub.x(1) > x(1) msub.x(end) < x(end) ...
            msub.y(end)> y(end) msub.y(1) < y(1) ];
        
        %read buffered data
        msub=readGeotiff(matchtagFileName,...
            'map_subset',[x(1) x(end) y(end) y(1)] +buff);
        
        
        if any(~msub.z(:))
            % if all true, set to true
            m(:,:,i) = true;
        elseif any(msub.z(:)) && any(~msub.z(:))
            % interpolate if not on same grid
            if length(msub.x) ~= length(x) || length(msub.y) ~= length(y)|| ...
                    zsub.x(1) ~= msub.x(1) || ...
                    zsub.y(1) ~= mub.y(1)
                m(:,:,i) = interp2(msub.x,msub.y(:),msub.z,x,y(:),'*nearest');
            else
                m(:,:,i) = msub.z;
            end
        end
    end
end
