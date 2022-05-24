function batch_boundaryAdjustApply(tileNeighborIndexFile,varargin)

strt=1;
inc=1;
tiles=[];
quad_row=[];
quad_col=[];

if ~isempty(varargin)
    n=find(strcmpi(varargin,'start'));
    if ~isempty(n)
        strt=varargin{n+1};
    end
    n=find(strcmpi(varargin,'inc'));
    if ~isempty(n)
        inc=varargin{n+1};
    end

    n = find(strcmpi('tiles', varargin));
    if ~isempty(n)
        tiles = varargin{n+1};
    end

    n = find(strcmpi('quad_row', varargin));
    if ~isempty(n)
        quad_row = varargin{n+1};
    end
    n = find(strcmpi('quad_col', varargin));
    if ~isempty(n)
        quad_col = varargin{n+1};
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
    if ~isempty(tiles)
        if ~any(strcmp(get10mTileName(fileNames{i}), tiles))
            continue
        end
    end
    if ~isempty(quad_row) || ~isempty(quad_col)
        [qrow, qcol] = getQuadTileRowCol(fileNames{i});
        if ~isempty(quad_row) && qrow ~= quad_row
            continue
        end
        if ~isempty(quad_col) && qcol ~= quad_col
            continue
        end
    end

    fprintf('%d: applying dz0 to %s\n',i,fileNames{i})
    % get this tile name and set is2 file name from it
    boundaryAdjustApply(fileNames{i})
end

fprintf('complete\n')
