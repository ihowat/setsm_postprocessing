function [x,y,z,missingFlag,m] =extractSubGrid(fileNames,subx0,subx1,suby0,suby1,res)
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

verbose = false;
% make subset grid
x=subx0:res:subx1;
y=suby1:-res:suby0;
y=y(:);
z=nan(length(y),length(x),length(fileNames),'single');
m=[];
if nargout == 5
    m = false(size(z));
end

% subset extraction loop

missingFlag=false(size(fileNames));
for i=1:length(fileNames)
    
    if verbose
        fprintf('extracting subset %d of %d\n',i,length(fileNames))
    end
    % get a square subset of pixels surrounding this point
    zsub=readGeotiff(fileNames{i},'map_subset',[min(x),max(x),min(y),max(y)]);

    % skip if too small a subset or no data
    if ~any(zsub.z(:) ~= -9999) || any(size(zsub.z) < 4)
        missingFlag(i)=true;
        continue
    end
    
    % convert nodata to NaN for interpolation
    zsub.z(zsub.z == -9999) = NaN;
    
    if length(zsub.x) ~= length(x) || length(zsub.y) ~= length(y)
        z(:,:,i) = interp2(zsub.x,zsub.y(:),zsub.z,x,y(:),'*linear');
    else
        z(:,:,i) = zsub.z;
    end
    
    if exist('m','var')
        if res == 10
            matchtagFileName = strrep(fileNames{i},'dem_10m.tif','matchtag.tif');
        elseif res == 2
            matchtagFileName = strrep(fileNames{i},'dem.tif','matchtag.tif');
        end
        
        msub=readGeotiff(matchtagFileName,'map_subset',[min(x),max(x),min(y),max(y)]);
        
        
        if length(msub.x) ~= length(x) || length(msub.y) ~= length(y)
            m(:,:,i) = interp2(msub.x,msub.y(:),msub.z,x,y(:),'*nearest');
        else
            m(:,:,i) = msub.z;
        end
    end
    
end
