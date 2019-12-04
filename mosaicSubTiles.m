function mosaicSubTiles(subTileDir)
% '15-Nov-2019 12:44:26': changed from mean to median of heightd
% differences for vertical alignment

%regionDir = '26_29';

% change upper level path if coming in from mac
if ismac
    uppath='/Users/ihowat';
else
    uppath='/home/howat.4';
end

if strcmp(subTileDir(end),'/')
	subTileDir(end)=[];
end

outName=[subTileDir,'.mat'];
subTileFiles=dir([subTileDir,'/*.mat']);
subTileFiles=cellfun( @(x) [subTileDir,'/',x],{subTileFiles.name},'uniformoutput',0);

if exist(outName,'file')
    % need to add subtile count to output for restart
    fprintf('%s exists, skipping\n',outName)
end

[~,subtileNum] = cellfun(@fileparts,subTileFiles,'uniformoutput',0);
subtileNum = cellfun(@(x) str2num(x(7:end)),subtileNum);

[~,n] = sort(subtileNum);
subTileFiles = subTileFiles(n);

% loop through tiles and get coordinate range
x0 = nan(size(subTileFiles)); x1=x0; y0 = x0; y1 = x0; dx = x0; dy =x0;
for n = 1:length(subTileFiles)
    
    clear x y
    load(subTileFiles{n},'x','y')
    x0(n) = min(x);
    x1(n) = max(x);
    y0(n) = min(y);
    y1(n) = max(y);
    dx(n) = x(2)-x(1);
    dy(n) = y(2)-y(1);
    
end

if length(unique(dx)) > 1 || length(unique(dy)) > 1
    error('multiple dx, dy sizes for subgrids')
end

dx = dx(1);
dy = dy(1);

% build final grid
x = min(x0):dx:max(x1);
y = max(y1):dy:min(y0);
z = nan(length(y),length(x));
N = z;

Nn = sum(~isnan(z(:)));
for filen=1:length(subTileFiles)

    fprintf('adding subtile %d of %d\n',filen,length(subTileFiles))

    % load next tile
    
    m  = who('-file',subTileFiles{filen});  

    if ~any(ismember(m,'za'))
        fprintf('no za found, skipping\n')
        continue
    end	
 
    if any(ismember(m,'za_med'))
        z2=load(subTileFiles{filen},'x','y','za_med','N','land');
    else
        z2=load(subTileFiles{filen},'x','y','za','land','fa');
        
        %apply filter
        z2.za(~z2.fa) = NaN;
        
        %fa = pairwiseDifferenceFilter(z2.za,'mask',z2.land,'minmad',1);
        %save(subTileFiles{filen},'fa','-append');
        %z2.za(~fa) = NaN;
        
        % get medians and number of repeats
        z2.za_med = nanmedian(z2.za,3);
        
        z2.N=sum(~isnan(z2.za),3);
        
        % remove unneeded data to free memory
        z2 = rmfield(z2,'za');
        z2 = rmfield(z2,'fa');
    end
  
    %filter out water
    if ~isfield(z2,'land')
	z2.land=true(size(z2.za_med));
    end

    z2.za_med(~z2.land) = NaN;
  
    if  ~any(z2.N(:)) || ~any(z2.land(:))
        fprintf('no data in tile, skipping\n')
        continue
    end
    
    % get corner indexes for this grid in the master
    col0  = find(z2.x(1) == x);
    col1  = find(z2.x(end) == x);
    row0  = find(z2.y(1) == y);
    row1  = find(z2.y(end) == y);
    
    % subet the current master grid to this one
    z1 = z(row0:row1,col0:col1);
    N1 = N(row0:row1,col0:col1);
    
    % if overlapping existing data
    if any(~isnan(z1(:)))
        
        dz = nanmedian(z2.za_med(~isnan(z1(:))) -z1(~isnan(z1(:))));
        
        if ~isnan(dz)
            
            f = 1./filen;
            
            z2.za_med = z2.za_med - f.*dz;
            z = z + (1-f).*dz;
            z1 = z1 + (1-f).*dz;
            
            % fill missing data in this tile with current data
            z2.za_med(isnan(z2.za_med) & ~isnan(z1)) = z1(isnan(z2.za_med) & ~isnan(z1));
            
            %determine number of rows/cols of overlap
            buffA = single(~(~isnan(z2.za_med) & ~isnan(z1)));
            buffA(~buffA) = NaN;
            
            buffA(1,isnan(buffA(1,:))) = 0;
            buffA(end,isnan(buffA(end,:))) = 0;
            
            buffA(isnan(buffA(:,1)),1) = 0;
            buffA(isnan(buffA(:,end)),end) = 0;
            
            buffA=inpaint_nans(double(buffA),2);
            
            notMissing = ~isnan(z1);
            
            z2.za_med(notMissing) = z2.za_med(notMissing).*buffA(notMissing) +...
                z1(notMissing).*(1- buffA(notMissing));
            
        else
	
		z2.za_med(~isnan(z1(:))) = z1(~isnan(z1(:)));

	end
        
    end
    
    z(row0:row1,col0:col1) = z2.za_med;
    N(row0:row1,col0:col1) = z2.N;
    
    Nn1 = sum(~isnan(z(:)));
    
    if Nn1 < Nn
        error('more nans in mosaic with this iteration\n');
    end
    Nn = Nn1;
    
    % imagesc(z2.za_med - z(rows2,cols2))
    %     clf
    %     imagesc(x,y,z); axis xy equal
    %     drawnow
    %    % input('')
end

save(outName,'x','y','z','N')
z(isnan(z)) = -9999;
outNameTif = strrep(outName,'.mat','.tif');
writeGeotiff(outNameTif,x,y,z,4,-9999,'polar stereo north')
