function batchRegisterTiles(tileDir,is2Dir,varargin)
% set paths

%     %   tileDir='/Users/ihowat/project/howat.4/rema_mosaic/rema_mosaic_v13e';
%     is2Dir='/Users/ihowat/data4/REMA/altimetryByTile';

%     %    tileDir='/fs/project/howat.4/howat.4/rema_mosaic_v13e/rema_19_victoria_land';
%     is2Dir='/fs/byo/howat-data4/REMA/altimetryByTile';

strt=1;
inc=1;
resolution='10m';

if length(varargin) == 2 && ~isnan(str2double(varargin{1})) && ~isnan(str2double(varargin{2}))
    strt=varargin{1};
    inc= varargin{2};

elseif ~isempty(varargin)
    n=find(strcmpi(varargin,'start'));
    if ~isempty(n)
        strt=varargin{n+1};
    end

    n=find(strcmpi(varargin,'inc'));
    if ~isempty(n)
        inc=varargin{n+1};
    end

    resolution_choices = {'10m', '2m'};
    n = find(strcmpi('resolution', varargin));
    if ~isempty(n)
        resolution = varargin{n+1};
    end
    if ~ismember(resolution, resolution_choices)
        error("'resolution' must be one of the following, but was '%s': {'%s'}", resolution, strjoin(resolution_choices, "', '"));
    end
end

fprintf('processing every %d tile, starting at %d\n',inc,strt)

tileFiles = dir(sprintf('%s/*_%s.mat', tileDir, resolution));
%tileFiles = dir([tileDir,'/*.mat']);
tileNames = {tileFiles.name};
tileFiles = fullfile({tileFiles.folder}, tileNames);
tileNames10m = cellfun(@(x) get10mTileName(x), tileNames, 'UniformOutput',false);

tileFiles=tileFiles(:);
tileNames10m=tileNames10m(:);

i=strt;
for i=strt:inc:length(tileFiles)

    % get this tile name and set is2 file name from it
    tileFile=tileFiles{i};
    if endsWith(tileFile, '_reg.mat')
        continue
    end
    tileName10m = tileNames10m{i};
    is2TileFile = [is2Dir,'/',tileName10m,'_is2.mat'];
    
    fprintf('%d of %d, %s\n',i,length(tileNames10m),tileFile)

    registerTileToIS2(tileFile,is2TileFile)
     
    applyRegistration(tileFile,[])
    
    regTileFile=strrep(tileFile,'.mat','_reg.mat');
    
    if exist(regTileFile,'file')
        fit2is2(regTileFile,is2TileFile)
    end
end
