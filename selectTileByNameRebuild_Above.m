function selectTileByName_Above(dstdir,tile_name,res,varargin)
% selectTilesByRegion mosaics strips covering tiles within the specified
% NGA ArcticDEM region

%tile_name='45_17';
%res=2;          % strip & tile region
projstr='canada albers equal area conic';

% file names
tilefile  = 'PGC_Imagery_Mosaic_Tiles_ABoVE_nocoast.mat'; %PGC/NGA Tile definition file
%dbasefile = ['REMAdatabase_',num2str(res),'m.mat']; % strip datbase file
dbasefile = ['aboveDEMdatabase_2m.mat']; % strip datbase file
%includeListFile = 'includeListRema.txt'; % list of strips to include

% make output directory
outdir=[dstdir,'/',tile_name];

if ~exist(outdir,'dir'); mkdir(outdir); end

%Get Arctic Tile Defs
tiles=load(tilefile);

% get target tile
[~,n]=intersect(tiles.I,tile_name);

% crop tile structure to overlapping tiles
tiles = structfun(@(x) ( x(n) ), tiles, 'UniformOutput', false);

% load database file into mat file object
meta=loadQcData(dbasefile);

% import and apply includeList if a file name exists
if exist('includeListFile','var');
    includeList=textread(includeListFile,'%s');
    includeList=strrep(includeList,'_8m_dem',''); % remove suffix to generalize
    includeList=strrep(includeList,'_2m_dem',''); % remove suffix to generalize
    
    % stripID from database structure to compare with includeList
    stripIDs = cellfun(@(x) x(find(x=='/',1,'last')+1:...
        end), meta.f, 'UniformOutput',false);
    
    stripIDs=strrep(stripIDs,['_',num2str(res),'m_meta.txt'],'');
    
    % compare includeList with stripIDs
    [~,n]=intersect(stripIDs,includeList);
    
    % exit if empty
    if isempty(n); error('all strips excluded');  end
    
    % remove excluded data from database
    meta = structfun(@(x) ( x(n) ), meta, 'UniformOutput', false);

end

% get gcps into gcp.x, y, and z
if ~isempty(varargin);
    gcpfile=varargin{1};
    gcp = loadGCPFile_is(gcpfile);
    % send to mosaicker
    mosaicStripsRebuild(meta,tiles,res,outdir,projstr,gcp);

else
    % send to mosaicker
    mosaicStripsRebuild(meta,tiles,res,outdir,projstr);
    
end

%% Crop Buffers and Write Tiles To Geotiffs 
tilef=dir([outdir,'/',tile_name,'_',int2str(res),'m_dem.mat']);
tilef = cellfun(@(x) [outdir,'/',x], {tilef.name}, 'UniformOutput',false);
tilef = tilef{1};
fprintf('writing tiffs from %s\n',tilef);

if exist(strrep(tilef,'dem.mat','reg_dem.mat'),'file')
    fi=strrep(tilef,'dem.mat','reg_dem.mat');
else
    fi=tilef;
end

fprintf('source: %s\n',fi);

% calc buffer to remove
buffer = floor(200 / res);

load(fi,'x','y');
% crop buffer tile
x=x(buffer+1:end-buffer);
y=y(buffer+1:end-buffer);

OutDemName = strrep(fi,'.mat','.tif');
if ~exist(OutDemName,'file');
    load(fi,'z');
    z=z(buffer+1:end-buffer,buffer+1:end-buffer);
    z(isnan(z)) = -9999;
    writeGeotiff(OutDemName,x,y,z,4,-9999,projstr)
    clear z
end

OutMatchtagName = strrep(fi,'dem.mat','matchtag.tif');
if ~exist(OutMatchtagName,'file');
    load(fi,'mt');
    mt =mt(buffer+1:end-buffer,buffer+1:end-buffer);
    writeGeotiff(OutMatchtagName,x,y,mt,1,0,projstr)
    clear mt
end

hillshade=strrep(OutDemName,'dem.tif','dem_shade.tif');
if ~exist(hillshade,'file');
    system(['gdaldem hillshade -compute_edges -b 1 -q -of GTiff -co tiled=yes -co compress=lzw -co bigtiff=if_safer ',...
        OutDemName,' ',hillshade]);
end

