
if ismac
    addpath('/Users/ihowat/unity-home/demtools');
    tileDir=['/Users/ihowat/project/howat.4/earthdem/earthdem_mosaic_testing_1km/'];
    waterTileDir='~/data/pgc_projects/ak_water_rasters_v2';
    gcpFile='~/GLA14_rel634.csv';
else
    %addpath('/home/howat.4/demtools');
    tileDir=['/mnt/pgc/data/scratch/claire/pgc/arcticdem/mosaic/2m_v4.2/test'];
    waterTileDir='/mnt/pgc/data/scratch/claire/pgc/arcticdem/coastline/global_surface_water/tiles_v2';
    gcpFile='/mnt/pgc/data/scratch/claire/pgc/arcticdem/gcp/icesat/mat/GLA14_rel634.mat';
end

% set to true to overwite existing offsets, false will skip
overwriteFlag=false;

% cellstr of tile file path/names
tileFiles = dir([tileDir,'/*2m.mat']);
tileFiles = cellfun( @(x) [tileDir,'/',x], {tileFiles.name},'uniformoutput',0);

% load gcp file
[~,~,ext] = fileparts(gcpFile);
if strcmpi(ext,'.csv')
    gcp=loadGCPFile_is(gcpFile);
elseif  strcmpi(ext,'.mat')
    gcp=load(gcpFile);
else
    error('file extension for %s not recognized',gcpFile)
end

% tile file loop
i=1;
for i=1:length(tileFiles)
    
    fprintf('processing %d of %d: %s ',i,length(tileFiles),tileFiles{i})
    
    % load tile
    m=matfile(tileFiles{i},'Writable',true);
    
    if ~overwriteFlag
        % skip if outputs already exist in this tile file
        if ismember('icesat_dz_med',fields(m))
            fprintf('icesat_dz_med already exists in file, skipping\n')
            continue
        end
    end
    
    % load coordinate vectors
    x=m.x;
    y=m.y;
    
    % get coordinate ranges
    x0 = min(x);
    x1 = max(x);
    y0 = min(y);
    y1 = max(y);
    
    
    [~,tileName] = fileparts(tileFiles{i});
    tileName = tileName(1:5);
    
    % get water mask
    waterMask = getTileWaterMask(waterTileDir,tileName,x0-10,x1+10,y0-10,y1+10,10);
    
    % interpolate 10m water mask to 2m matching the dem tile
    waterMask = interp2(waterMask.x,waterMask.y(:),waterMask.z,x,y(:),'*nearest');
    
    % get median difference and stats
    [dz_med,dz_mad,N] = diffDEMFromGCPs(m,gcp,waterMask);
    
    fprintf(' dz_med=%.2f, dz_mad=%.2f, N=%d\n',dz_med,dz_mad,N)
    
    m.icesat_dz_med = dz_med;
    m.icesat_dz_mad = dz_mad;
    m.icesat_dz_N = N;
    
    if ~isnan(dz_med)
        z=  m.z;
        z = z - dz_med;
        m.z = z;
        clear z
    end
    
end





