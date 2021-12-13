function batchRegisterTiles(tileDir,is2Dir,varargin)
% set paths

%     %   tileDir='/Users/ihowat/project/howat.4/rema_mosaic/rema_mosaic_v13e';
%     is2Dir='/Users/ihowat/data4/REMA/altimetryByTile';

%     %    tileDir='/fs/project/howat.4/howat.4/rema_mosaic_v13e/rema_19_victoria_land';
%     is2Dir='/fs/byo/howat-data4/REMA/altimetryByTile';

strt=1;
inc=1;
if length(varargin) == 2
    strt=varargin{1};
    inc= varargin{2};
    fprintf('processing every %d tile, starting at %d\n',inc,strt)
end

tileFiles = dir([tileDir,'/*_10m.mat']);
tileNames = {tileFiles.name};
tileFiles = cellfun( @(x) [tileDir,'/',x], tileNames,'uniformoutput',0);
tileNames = strrep(tileNames,'_10m.mat','');

tileFiles=tileFiles(:);
tileNames=tileNames(:);

i=strt;
for i=strt:inc:length(tileFiles)
    
    % get this tile name and set is2 file name from it
    tileName = tileNames{i};
    tileFile=tileFiles{i};
    is2TileFile = [is2Dir,'/',tileName,'_is2.mat'];
    
    fprintf('%d of %d, %s\n',i,length(tileNames),tileFile)

%     registerTileToIS2(tileFile,is2TileFile)
%     
%     applyRegistration(tileFile,[])
    
    regTileFile=strrep(tileFile,'.mat','_reg.mat');
    
    if exist(regTileFile,'file')
        fit2is2(regTileFile,is2TileFile)
    end
end
    
    
    