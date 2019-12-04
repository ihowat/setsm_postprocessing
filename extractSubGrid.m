function [x,y,z,missingFlag] =extractSubGrid(fileNames,subx0,subx1,suby0,suby1,res)

verbose = false;
% make subset grid
x=subx0:res:subx1;
y=suby1:-res:suby0;
y=y(:);

% subset extraction loop

z=nan(length(y),length(x),length(fileNames),'single');
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
    
end