function selectTileByName(dstdir,tile_name,res,varargin)
% selectTilesByRegion mosaics strips covering tiles within the specified
% NGA ArcticDEM region

%tile_name='45_17';
%res=2;
projstr='polar stereo north';

% file names
tilefile  = 'PGC_Imagery_Mosaic_Tiles_Arctic.mat'; %PGC/NGA Tile definition file
dbasefile = ['arcticDEMdatabase_2m.mat']; % strip datbase file
%dbasefile = ['GrITdatabase_2m.mat']; % strip datbase file
%includeListFile = 'includeList_baffin.txt'; % list of strips to include

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
    %gcp = loadGCPFile_is(gcpfile);
    gcp = loadGCPFile_is(gcpfile);
    % send to mosaicker
    mosaicStrips(meta,tiles,res,outdir,projstr,gcp);

else
    % send to mosaicker
    mosaicStrips(meta,tiles,res,outdir,projstr);
    
end