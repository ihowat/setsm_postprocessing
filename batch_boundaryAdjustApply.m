function batch_boundaryAdjustApply(tileNeighborIndexFile,varargin)

strt=1;
inc=1;
if ~isempty(varargin)
    n=find(strcmpi(varargin,'start'));
    if ~isempty(n)
        strt=varargin{n+1};
    end
    
    n=find(strcmpi(varargin,'inc'));
    if ~isempty(n)
        inc=varargin{n+1};
    end
end

fprintf('processing every %d tile, starting at %d\n',inc,strt)

load(tileNeighborIndexFile,'fileNames')

if ismac
    fileNames=strrep(fileNames,'/fs/project/howat.4/','/Users/ihowat/project/');
end

%% dz0 calc and write loop - this could be run simultaneously/parallel
i=strt;
for i=strt:inc:length(fileNames)
    fprintf('%d: applying dz0 to %s\n',i,fileNames{i})
    % get this tile name and set is2 file name from it
    boundaryAdjustApply(fileNames{i})
end

fprintf('complete\n')

    
