function mosaicSubTiles(subTileDir,dx,x0,x1,y0,y1,outName)
% mosaicSubTiles mosaic subtiles and write mat and geotiff output
%
% mosaicSubTiles(subTileDir,dx,x0,x1,y0,y1,outName) mosaics all of the
% subtile .mat files in the directory subTileDir into a mosaic with grid
% resolution dx (2 or 10) and rectangular coordinate range x0 to x1 and y0
% to y1. Writes matfile output to outName, with tiff output as
% strrep(outName,'.mat','.tif')
%
% version 1:
% Subtile alignment is from median difference between that subtile and
% the mosaic, weighted by the the number of subtiles arleady added - so the
% subtile is shifted by dz*(1-f) and the mosaic is shifted by dz*f, where f
% = 1/N and N is the Nth subtile begin added. I tried added in descending
% order of number strips in subtile, but problems arose by matches at
% corner overlaps. A bundle adjustment procedure is in development.

fprintf('Indexing subtiles\n')

% make a cellstr of resolved subtile filenames
subTileFiles=dir([subTileDir,'/*_',num2str(dx),'m.mat']);
subTileFiles=cellfun( @(x) [subTileDir,'/',x],{subTileFiles.name},'uniformoutput',0);

% Get column-wise number of subtile from file names - assumes the subtile
% names is {tilex}_{tily}_{subtilenum}_....
[~,subTileName] = cellfun(@fileparts,subTileFiles,'uniformoutput',0);
subTileName=cellfun(@(x) strsplit(x,'_'),subTileName,'uniformoutput',0);
subTileNum = cellfun(@(x) str2num(x{3}),subTileName);

% sort subtilefiles by ascending subtile number order
[subTileNum,n] = sort(subTileNum);
subTileFiles = subTileFiles(n);
NsubTileFiles = length(subTileFiles);

% find buffer size from first 2 tiles
n = diff(subTileNum);
n(mod(subTileNum(1:end-1),100) == 0) = 0;
n = find(n == 1,1,'first');
if ~isempty(n)
    buffcheck1=load(subTileFiles{n},'y');
    buffcheck2=load(subTileFiles{n+1},'y');
    buff = (length(buffcheck2.y)-find(buffcheck1.y(1) == buffcheck2.y))/2;
    buff = round(buff);

    fprintf('Using subtile buffer of %d pixels\n',buff)
else
    fprintf('Too few subtiles exist to determine buffer, quitting')
    return
end
    
% make a polyshape out of boundary for checking subtile overlap
tilePoly = polyshape([x0 x0 x1 x1]',[y0 y1 y1 y0]');

% build tile coordinate vectors
x = x0:dx:x1;
y = y1:-dx:y0;

% build tile output arrays
z = nan(length(y),length(x));
N = zeros(length(y),length(x),'uint8');

if dx == 2
    z_mad = z;
    tmax = zeros(length(y),length(x),'uint16');
    tmin =zeros(length(y),length(x),'uint16');
end

% initialize subtile count for use in n-weighted alignment
subtile_n=1;

% initialize count of pixels with data in mosaic for error checking
Nn=0;
for filen=1:NsubTileFiles
    
    fprintf('adding subtile %d of %d: %s\n',filen,NsubTileFiles,subTileFiles{filen})
    
    % get list of variables within this subtile mat file
    mvars  = who('-file',subTileFiles{filen});
    
    if ~any(ismember(mvars,'za_med'))
        fprintf('no za_med found, skipping\n')
        continue
    end
    
    % Check if subtile overlaps mosaic boundary:
    % open matfile to load coordinate vectors withoutloading arrays
    m=matfile(subTileFiles{filen});
    
    % make a polyshape out of this subtile boundary
    filenPoly = polyshape([min(m.x) min(m.x) max(m.x) max(m.x)]',...
        [min(m.y) max(m.y) max(m.y) min(m.y)]');
    
    % test overlap
    if ~overlaps(tilePoly,filenPoly)
        fprintf('subtile out of bounds, skipping\n')
        continue
    end
    
    % load subtile into structure
    if dx == 2
        if ismember(mvars,'land')
            zsub=load(subTileFiles{filen},'x','y','za_med','land','N','za_mad','tmax','tmin');
        else
            % if 2m subtile .mat doesnt include a land array, load it from
            % the 10m subtile and resize it
            zsub=load(subTileFiles{filen},'x','y','za_med','N','za_mad','tmax','tmin');
            subTileFile10m = strrep(subTileFiles{filen},'_2m','_10m');
            land = load(subTileFile10m,'land');
            if any(~land.land(:))
                zsub.land = imresize(land.land,size(zsub.za_med),'nearest');
                clear land subTileFile10m
            else
                zsub.land = true(size(zsub.za_med));
            end
        end
    elseif dx == 10
        zsub=load(subTileFiles{filen},'x','y','za_med','N','land');
    end

    if ~any(~isnan(zsub.za_med(:)))
        fprintf('all nans, skipping\n')
        continue
    end
    
    % get corner indexes for this grid in the master
    col0  = find(zsub.x(1) == x);
    col1  = find(zsub.x(end) == x);
    row0  = find(zsub.y(1) == y);
    row1  = find(zsub.y(end) == y);
    
    if isempty(col0) || isempty(col1) || isempty(row0) || isempty(row1)
        fprintf('subtile and tile coordinates dont match, subtile may span boundary,skipping\n')
        continue
    end
    
    % subset the current mosaic by range of subtile
    z1 = z(row0:row1,col0:col1);
    N1 = N(row0:row1,col0:col1);
    
    if dx == 2
        z_mad1 = z_mad(row0:row1,col0:col1);
        tmax1 = tmax(row0:row1,col0:col1);
        tmin1 = tmin(row0:row1,col0:col1);
    end
    
    % find overlapping non-nan and non-water pixels
    n_overlap = ~isnan(z1(:)) & ~isnan(zsub.za_med(:)) & zsub.land(:);
    
    % check if overlapping pixels exist to determine if blending is needed
    if any(n_overlap)
        
        %get median difference between this subtile and mosaic subset
        dz_med = median(zsub.za_med(n_overlap) - z1(n_overlap));
        
        % subtile number-weighted alignement
        f = 1./subtile_n; % weight of this subtile relative to the mosaic 
        zsub.za_med = zsub.za_med - (1-f).*dz_med; % shift the subtile
        z = z + f.*dz_med; % shift the tile
        z1 = z1 + f.*dz_med; % shift the tile subset
        
        % fill missing data in this subtile with current data in mosaic
        % subset
        zsub.za_med(isnan(zsub.za_med) & ~isnan(z1)) =...
            z1(isnan(zsub.za_med) & ~isnan(z1));
        
        if dx == 2
            % also blend za_mad
            zsub.za_mad(isnan(zsub.za_mad) & ~isnan(z_mad1)) =...
                z1(isnan(zsub.za_mad) & ~isnan(z_mad1));
        end
        
        % Create the blending array by setting zeros at the far edge of
        % subtile/tile overlap and ones at the other edge, linearly
        % interpolating between the two.
        % find where data is missing in both subtile and tile
        buffA = single(~(~isnan(zsub.za_med) & ~isnan(z1)));
        
        % set pixels with data in both subtile and tiles to NaN as an
        % interpolation flag for inpaint_nans
        buffA(~buffA) = NaN;
        
        % set boundaries of blend array to zero
        buffA(1,isnan(buffA(1,:))) = 0;
        buffA(end,isnan(buffA(end,:))) = 0;
        
        buffA(isnan(buffA(:,1)),1) = 0;
        buffA(isnan(buffA(:,end)),end) = 0;
        
        % interpolate linearly across NaNs 
        buffA=inpaint_nans(double(buffA),2);
        
        % find where there is data in the tile subset
        notMissing = ~isnan(z1);
        
        % blend the subtile and mosaic subset, applying the edge-distance
        % weighting
        zsub.za_med(notMissing) = zsub.za_med(notMissing).*buffA(notMissing) +...
            z1(notMissing).*(1- buffA(notMissing));
        
        if dx == 2
            zsub.za_mad(notMissing) = zsub.za_mad(notMissing).*buffA(notMissing) +...
                z_mad1(notMissing).*(1- buffA(notMissing));
        end
        
    else
        
        % if no pixels overlap, just add the subset data into the subtile,
        % replacing the NaNs
        zsub.za_med(~isnan(z1(:))) = z1(~isnan(z1(:)));
 
        if dx == 2
            zsub.za_mad(~isnan(z1(:))) = z_mad1(~isnan(z1(:)));
        end
        
    end
    
    % place the belended substile into the tile
    z(row0:row1,col0:col1) = zsub.za_med; 
    
    % place N grid into the tile just within the tile borders
    N(row0+buff:row1-buff,col0+buff:col1-buff) =...
        zsub.N(buff+1:end-buff,buff+1:end-buff);
    
    if dx == 2
        z_mad(row0:row1,col0:col1) = zsub.za_mad;
        tmax(row0+buff:row1-buff,col0+buff:col1-buff) =...
            zsub.tmax(buff+1:end-buff,buff+1:end-buff);
        tmin(row0+buff:row1-buff,col0+buff:col1-buff) =...
            zsub.tmin(buff+1:end-buff,buff+1:end-buff);
    end
    
    % count the number of pixels with data after this merge
    Nn1 = sum(~isnan(z(:)));
    
    % if the new number of pixels with data is now less than before,
    % then data was overwritten with NaNs, which is an error.
    if Nn1 < Nn
        error('more nans in mosaic with this iteration\n');
    end
    
    % reset count of pixels with data for next iteration
    Nn = Nn1;
    
    % update count of subtiles added to mosaic
    subtile_n= subtile_n+1;

end

if ~(isempty(nonzeros(N)))
    % save matfile outputs 
    if dx == 2
        save(outName,'x','y','z','N','z_mad','tmax','tmin','-v7.3')
    elseif dx == 10
        save(outName,'x','y','z','N','-v7.3')
    end

    % write tiff files
    z(isnan(z)) = -9999;
    outNameTif = strrep(outName,'.mat','_dem.tif');
    writeGeotiff(outNameTif,x,y,z,4,-9999,'polar stereo north')

    outNameTif = strrep(outName,'.mat','_N.tif');
    writeGeotiff(outNameTif,x,y,N,1,0,'polar stereo north')

    if dx == 2
        z_mad(isnan(z_mad)) = -9999;
        outNameTif = strrep(outName,'.mat','_mad.tif');
        writeGeotiff(outNameTif,x,y,z_mad,4,-9999,'polar stereo north')

        outNameTif = strrep(outName,'.mat','_tmax.tif');
        writeGeotiff(outNameTif,x,y,tmax,2,0,'polar stereo north')

        outNameTif = strrep(outName,'.mat','_tmin.tif');
        writeGeotiff(outNameTif,x,y,tmin,2,0,'polar stereo north')
    end
else
    fprintf('N array is empty, outputs not written\n');
    % Write semophone file so that rerunning does not try this tile again
    outNameSem = strrep(outName,'.mat','_empty.txt');
    fileID = fopen(outNameSem,'w');
    fclose(fileID);
end
