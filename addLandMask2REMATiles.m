function addLandMask2REMATiles(tileDir,tileDefFile,varargin)

% fileNames=dir('*10m_reg.mat');
% fileNames={fileNames.name};

% if ismac
%     tileDefFile = '/Users/ihowat/data4/REMA/rema_tile_definitions.mat';
% else
%     tileDefFile = '/fs/byo/howat-data4/REMA/rema_tile_definitions.mat';
% end


%fileNames=dir([tileDir,'/*_10m.mat']);
fileNames=dir([tileDir,'/*.mat']);
fileNames = fullfile({fileNames.folder}, {fileNames.name});

tileDefs=load(tileDefFile);
if isfield(tileDefs,'I') && ~isfield(tileDefs,'tileName')
    tileDefs.tileName=tileDefs.I;
end


strt=1;
inc=1;
if length(varargin) == 2
    strt=varargin{1};
    inc= varargin{2};
    fprintf('processing every %d tile, starting at %d\n',inc,strt)
end


j=strt;
for j=strt:inc:length(fileNames)
    fprintf('%d: adding land mask to %s...',j,fileNames{j})
    addLandMask2REMATile(fileNames{j},tileDefs)
    fprintf('done\n')
end

fprintf('completed\n')




