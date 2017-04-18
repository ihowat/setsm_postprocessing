function selectTileByName_Antarctic(dstdir,tile_name,res,varargin)
% selectTilesByRegion mosaics strips covering tiles within the specified
% NGA ArcticDEM region

%tile_name='45_17';
%res=2;          % strip & tile region
projstr='polar stereo south';

% file names
tilefile  = 'PGC_Imagery_Mosaic_Tiles_Antarctic.mat'; %PGC/NGA Tile definition file
%dbasefile = ['REMAdatabase_',num2str(res),'m.mat']; % strip datbase file
dbasefile = ['REMAdatabase_2m.mat']; % strip datbase file
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
    mosaicStrips(meta,tiles,res,outdir,projstr,gcp);

else
    % send to mosaicker
    mosaicStrips(meta,tiles,res,outdir,projstr);
    
end
