function update10mSubTileOutput(subTileFiles)
% update10mSubTileOutput add matchtag, stats and time arrays
%
% update10mSubTileOutput(subTileFiles) adds the arrays 'mt','mta','za_mad',
% 'Nmt','tmax' and 'tmin' to the 10m subtile files to make them consistent
% with the 2m files.

% convert to input to cell if string
if ~iscell(subTileFiles)
    subTileFiles = {subTileFiles};
end

for n=1:length(subTileFiles)
    
    clear fileNames fileNames0 x y za dZ dX dY fa mt mta Nmt za_mad tmin tmax 
    
    fprintf('subtile %d of %d\n',n,length(subTileFiles))
        
    outName = subTileFiles{n};
    
    % check to see if already updated & skip if yes
    mvars  = who('-file',outName);
    if any(ismember(mvars,'za_mad'))
      fprintf('za_mad found, already updated, returning\n')
       return
     end
     
     if any(~ismember({'za','za','dZ','dX','dY','fa'},mvars))  
       fprintf('%s missing variables,skipping\n',outName)
       continue
     end
     
    % load needed fields
    load(outName,'fileNames','fileNames0','x','y','za','dZ','dX','dY','fa')
    
    % path change for laptop
    if ismac
        fileNames0 = strrep(fileNames0,'/fs/byo/howat-','/Users/ihowat/');
    end
    
    % extract the matchtag stack
    [~,~,mt] =extractMatchtagSubGrid(fileNames0,min(x),max(x),...
        min(y),max(y),10);

    % merge segments from same strips
    % dont get offsets between segs in same strip:make a vector of z's
    % that ar belonging to the same strip
    [~,stripid] =  cellfun(@fileparts,fileNames0,'uniformoutput',0);
    stripid =  cellfun(@(x) x(1:47),stripid,'uniformoutput',0);
    unique_stripid = unique(stripid);
    
    r=[];
    it=1;
    for it=1:length(unique_stripid)
        
        segs = find(strcmp(stripid,unique_stripid(it)));
        
        if length(segs) > 1
               
            mt(:,:,segs(1)) = any(mt(:,:,segs),3);
            
            segs=segs(:);
            r = [r(:);segs(2:end)];
        end
    end
    
    mt(:,:,r) = [];
    clear it segs r
    
    % make adjusted mt arrays
    mta = false(size(mt));
    
    for k=1:size(mt,3)
        
        if ~isnan(dZ(k))
            mtak = interp2(x + dX(k),y + dY(k), single(mt(:,:,k)),...
                x,y,'*nearest');
        
            mtak(isnan(mtak)) = 0;
            mtak = logical(mtak);
        
            mta(:,:,k) = mtak;
        end
    end
    
    za(~fa) = NaN;
    mta(~fa) = false;
    za_mad = mad(za,1,3);
    Nmt = uint8(sum(mta,3));

    
    % make date vector
    [~,name] =  cellfun(@fileparts,fileNames,'uniformoutput',0);
    t=cellfun(@(x) datenum(x(6:13),'yyyymmdd'),name)';
    
    
    t=t-datenum('1/1/2000 00:00:00');
    t=reshape(t,1,1,[]);
    t=repmat(t,size(za_mad));
    t(isnan(za))=NaN;
    tmax = max(t,[],3);
    tmin = min(t,[],3);
    %tmean = nanmean(t,3);
    
    tmax = uint16(tmax);
    tmin = uint16(tmin);
    %tmean = uint16(tmean);
    
    
    fprintf('saving za_mad mt mta Nmt tmax tmin to %s\n',outName)
    save(outName,'mt','mta','za_mad','Nmt','tmax','tmin','-append');
    
    
end


function [x,y,m] = extractMatchtagSubGrid(fileNames,subx0,subx1,suby0,suby1)
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

% make subset grid
res=10;
x=subx0:res:subx1;
y=suby1:-res:suby0;
y=y(:);
m=false(length(y),length(x),length(fileNames));


% subset extraction loop
for i=1:length(fileNames)
    
    if contains(fileNames{i},'_smooth')
            matchtagFileName = strrep(fileNames{i},'dem_smooth_10m.tif','matchtag_10m.tif');
    else
            matchtagFileName = strrep(fileNames{i},'dem_10m.tif','matchtag.tif');        
    end
    
    msub=readGeotiff(matchtagFileName,'map_subset',[min(x),max(x),min(y),max(y)]);

    if length(msub.x) ~= length(x) || length(msub.y) ~= length(y)
        m(:,:,i) = interp2(msub.x,msub.y(:),msub.z,x,y(:),'*nearest');
    else
        m(:,:,i) = msub.z;
    end

end
