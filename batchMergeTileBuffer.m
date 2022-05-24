function batchMergeTileBuffer(tileNeighborIndexFile, varargin)
% batchMergeTileBuffer: mergeTileBuffer to all neighboring files
%
% batchMergeTileBuffer((tileNeighborIndexFile) where tileNeighborIndexFile 
% contains the list of fileNames to be merged (fileNames) and the tile
% neighbor index array nN=[indTop, indBottom, indLeft, indRight], created
% by tileNeighborIndex

tiles=[];
quad_row=[];
quad_col=[];
rows=false;
cols=false;

if ~isempty(varargin)
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

    n = find(strcmpi('rows', varargin));
    if ~isempty(n)
        rows = true;
    end
    n = find(strcmpi('cols', varargin));
    if ~isempty(n)
        cols = true;
    end
end

load(tileNeighborIndexFile,'fileNames','nN');

if ~isempty(tiles)
    tileNames10m = cellfun(@(x) get10mTileName(x), fileNames, 'UniformOutput',false);
    selected_tiles = zeros(size(fileNames));
    for i=1:length(tiles)
        selected_tiles = selected_tiles | strcmp(tileNames10m, tiles{i});
    end

    if ~isempty(quad_row) || ~isempty(quad_col)
        [file_quad_rows, file_quad_cols] = cellfun(@(x) getQuadTileRowCol(x), fileNames, 'UniformOutput',false);
        file_quad_rows = cell2mat(file_quad_rows);
        file_quad_cols = cell2mat(file_quad_cols);
        if ~isempty(quad_row)
            selected_tiles = selected_tiles & (file_quad_rows == quad_row);
        end
        if ~isempty(quad_col)
            selected_tiles = selected_tiles & (file_quad_cols == quad_col);
        end
    end
else
    selected_tiles = ones(size(fileNames));
end
if ~isequal(size(selected_tiles), size(nN))
    selected_tiles = selected_tiles.';
end

% find files where indices of top and right files exist
indTopExists = find(selected_tiles & ~isnan(nN(:,1)));
indRightExists = find(selected_tiles & ~isnan(nN(:,4)));

if cols
    n0 = indTopExists;
    n1 = nN(indTopExists,1);
elseif rows
    n0 = indRightExists;
    n1 = nN(indRightExists,4);
else
    % concatenate vectors of indices of bottom/top and left/right pairs
    n0 = [indTopExists; indRightExists]; %[bottom file0; left file0]
    n1 = [nN(indTopExists,1); nN(indRightExists,4)]; %[top file1; right file1]
end

% run mergeTileBuffer on each pair in list
for i=1:length(n0)
    mergeTileBuffer(fileNames{n0(i)},fileNames{n1(i)});
end

% If running a subset of tiles, make sure all border tiles
% have their edges merged to the outside world by checking
% and running the bottom and left merges if necessary.
if ~isempty(tiles)
    % find files where indices of bottom and left files exist
    indBottomExists = find(selected_tiles & ~isnan(nN(:,2)));
    indLeftExists = find(selected_tiles & ~isnan(nN(:,3)));

    if cols
        n0 = indBottomExists;
        n1 = nN(indBottomExists,2);
    elseif rows
        n0 = indLeftExists;
        n1 = nN(indLeftExists,3);
    else
        % concatenate vectors of indices of bottom/top and left/right pairs
        n0 = [indBottomExists; indLeftExists]; %[top file0; right file0]
        n1 = [nN(indBottomExists,2); nN(indLeftExists,3)]; %[bottom file1; left file1]
    end

    % run mergeTileBuffer on each pair in list
    for i=1:length(n0)
        mergeTileBuffer(fileNames{n0(i)},fileNames{n1(i)});
    end
end
