function checkLandPercent(tileDir,varargin)

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


strt=1;
inc=1;
if length(varargin) == 2
    strt=varargin{1};
    inc= varargin{2};
    fprintf('processing every %d tile, starting at %d\n',inc,strt)
end


fprintf('percent land coverage of tiles is as follows...\n')

j=strt;
for j=strt:inc:length(fileNames)
    m=matfile(fileNames{j});
    if ~any(strcmp(fields(m),'land'))
        landFraction=NaN;
    else
        landFraction=nnz(m.land)/numel(m.land);
    end
    fprintf('%s: %g\n',fileNames{j},landFraction)
end

fprintf('completed\n')




