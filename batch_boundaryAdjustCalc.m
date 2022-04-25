function batch_boundaryAdjustCalc(tileNeighborIndexFile,varargin)
% boundaryAdjust adjust tile based on offsets with neighbors over buffers
%
% boundaryAdjust(fileNames) where fileNames is a cell of full filenames of
% tile file. The eight surrounding tiles are
% loading and the offset over the buffer is caculated. The values are then
% upscaled by the resizeFraction, which should be ~ the fractional width
% of the buffer, and added along the edges of  a correction field (dz0)
% that is then filled  through interpolation. dz0 is generated for all
% tiles in the list and then, afterward, is applied to all tiles in the
% list. A flag, adusted, is added with value false if dz0 has not been
% applied and true if it has.
%
% boundaryAdjust(fileNames,'noApply')only calculates/saves the dz0 field
% but does not apply it to z.

%resizeFraction=0.1;
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

load(tileNeighborIndexFile,'fileNames','nN')

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

    fprintf('%d: Calculating offset grid for %s .... ',i,fileNames{i})
    
    if ~any(~isnan(nN(i,:)))
        fprintf('no neighboring tiles, skipping\n')
        continue
    end
    
    neighborFiles = cell(8,1);
    neighborFiles(~isnan(nN(i,:))) = fileNames(nN(i,~isnan(nN(i,:))));
    
    boundaryAdjustCalc(fileNames{i},neighborFiles)
    fprintf('done\n')
    clear  neighborFiles
end

fprintf('complete\n')
