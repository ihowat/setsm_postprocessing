function undoMergeTileBuffersSingle(tile_mat, edgeName, varargin)
% undoMergeTileBuffers reverts the tile edge buffers to the original values
%
% undoMergeTileBuffers(fileNames) replaces the merged tile buffers, created
% by mergeTileBuffer, with the original values, setting buffer flags to
% false.


%% Parse/validate positional args

tile_matfile = '';
m_struct = [];
if isa(tile_mat, 'matlab.io.MatFile')
    m_struct = tile_mat;
elseif isa(tile_mat, 'char') || isa(tile_mat, 'string')
    tile_matfile = tile_mat;
    if ~isfile(tile_matfile)
        error("Source tile .mat file does not exist: %s\n", tile_mat);
    end
    m_struct = matfile(tile_matfile);
else
    error("Unexpected class for argument 'tile_mat': '%s'", class(tile_mat));
end

buff_n_top = 1;
buff_n_bottom = 2;
buff_n_left = 3;
buff_n_right = 4;
edge_n = [];

edgeName_choices = {'all', 'top', 'bottom', 'left', 'right'};
if ~ismember(edgeName, edgeName_choices)
    error("'edgeName' must be one of the following, but was '%s': {'%s'}", edgeName, strjoin(edgeName_choices, "', '"));
end
if strcmp(edgeName, 'all')
    edge_n = [];
elseif strcmp(edgeName, 'top')
    edge_n = buff_n_top;
elseif strcmp(edgeName, 'bottom')
    edge_n = buff_n_bottom;
elseif strcmp(edgeName, 'left')
    edge_n = buff_n_left;
elseif strcmp(edgeName, 'right')
    edge_n = buff_n_right;
end


%% Parse/validate optional args

whichArrays_choices = {'all', 'zonly', 'all-but-z'};
n = find(strcmpi('whichArrays', varargin));
if ~isempty(n)
    whichArrays = varargin{n+1};
else
    whichArrays = 'all';
end
if ~ismember(whichArrays, whichArrays_choices)
    error("'whichArrays' must be one of the following, but was '%s': {'%s'}", whichArrays, strjoin(whichArrays_choices, "', '"));
end

n = find(strcmpi('ignoreMergedVars', varargin));
if ~isempty(n)
    ignoreMergedVars = varargin{n+1};
else
    ignoreMergedVars = false;
end


%% Verify state of tile struct

m_varlist = who(m_struct);

mergedTop = ismember('mergedTop', m_varlist) && m_struct.mergedTop;
mergedBottom = ismember('mergedBottom', m_varlist) && m_struct.mergedBottom;
mergedLeft = ismember('mergedLeft', m_varlist) && m_struct.mergedLeft;
mergedRight = ismember('mergedRight', m_varlist) && m_struct.mergedRight;

fprintf("\n");
fprintf("State of tile before undo:\n");
fprintf("mergedTop = %s\n", mat2str(mergedTop));
fprintf("mergedBottom = %s\n", mat2str(mergedBottom));
fprintf("mergedLeft = %s\n", mat2str(mergedLeft));
fprintf("mergedRight = %s\n", mat2str(mergedRight));
fprintf("\n");

data_array_names = {
    'z',
    'z_mad',
    'N',
    'Nmt',
    'tmin',
    'tmax',
};
exclude_array_names = {
    'land',
};
if ~all(ismember(data_array_names, m_varlist))
    data_array_names
    error("One or more expected data arrays do not exist in tile struct")
end

% Add data arrays that might be missing from the above list
z_size = size(m_struct, 'z');
data_array_names_same_sz = m_varlist(cellfun(@(name) isequal(size(m_struct, name), z_size), m_varlist));
data_array_names = setdiff(union(data_array_names, data_array_names_same_sz), exclude_array_names);

buff_array_names = cellfun(@(name) [name,'buff'], data_array_names, 'UniformOutput',false);

% Ensure 'merged*' variable states are consistent with presence of '*buff' arrays
if mergedTop || mergedBottom || mergedLeft || mergedRight
    if ~all(ismember(buff_array_names, m_varlist))
        buff_array_names
        error("Tile 'merged{Right|Left|Top|Bottom}=true', but some '*buff' array(s) does not exist")
    end
elseif strcmp(edgeName, 'all')
    fprintf("Tile 'merged{Right&Left&Top&Bottom}=false', indicating tile has not been merged: %s\n", tile_matfile);
    if ignoreMergedVars
        fprintf(strjoin(["Will attempt to undo merge with possibly existing '*buff' arrays",...
                         "because 'ignoreMergedVars=true'\n"]));
    else
        fprintf(strjoin(["Provide 'ignoreMergedVars=true' if you want to attempt to undo merge with",...
                         "possibly existing '*buff' arrays\n"]));
        return;
    end
end

% Check merged state of indicated edge and report for user FYI
if ~strcmp(edgeName, 'all')
    indicated_edge_not_merged = false;
    if strcmp(edgeName, 'top') && ~mergedTop
        indicated_edge_not_merged = true;
    elseif strcmp(edgeName, 'bottom') && ~mergedBottom
        indicated_edge_not_merged = true;
    elseif strcmp(edgeName, 'left') && ~mergedLeft
        indicated_edge_not_merged = true;
    elseif strcmp(edgeName, 'right') && ~mergedRight
        indicated_edge_not_merged = true;
    end
    if indicated_edge_not_merged
        fprintf("Tile 'merged%s=false', indicating tile has not been merged on this edge\n", edgeName);
        if ignoreMergedVars
            fprintf(strjoin(["Will attempt to undo merge with possibly existing '*buff' arrays",...
                             "because 'ignoreMergedVars=true'\n"]));
        else
            fprintf(strjoin(["Provide 'ignoreMergedVars=true' if you want to attempt to undo merge with",...
                             "possibly existing '*buff' arrays\n"]));
            return;
        end
    end
end


%% Undo buffer merges

m_struct.Properties.Writable = true;

unmerged_array_names = {};
unmerge_failure_general = false;
unmerge_failure_top = false;
unmerge_failure_bottom = false;
unmerge_failure_left = false;
unmerge_failure_right = false;

if strcmp(whichArrays, 'all')
    reset_array_names = data_array_names;
elseif strcmp(whichArrays, 'zonly')
    reset_array_names = data_array_names(strcmp('z', data_array_names));
elseif strcmp(whichArrays, 'all-but-z')
    reset_array_names = data_array_names(~strcmp('z', data_array_names));
end

for array_idx = 1:length(reset_array_names)
    data_array_name = reset_array_names{array_idx};     % ex. 'z'
    buff_array_name = [data_array_name,'buff'];         % ex. 'zbuff'
    if ~ismember(buff_array_name, m_varlist)
        continue;
    end
    buff_arrays = getfield(m_struct, buff_array_name);  % four buffer arrays, one for each edge

    if isempty(buff_arrays) || all(cellfun(@(x) isempty(x), buff_arrays))
        m_struct
        fprintf("Tile '%s' is empty; if '%s' is merged, cannot undo\n", buff_array_name, data_array_name);
        unmerge_failure_general = true;
        continue;
    end

    if ~strcmp(edgeName, 'all') && isempty(buff_arrays{edge_n})
        m_struct
        fprintf("Tile '%s(%d)' (%s side buffer) is empty; if '%s' is merged, cannot undo\n",...
            buff_array_name, edge_n, edgeName, data_array_name);
        unmerge_failure_general = true;
        continue;
    end

    data_array = getfield(m_struct, data_array_name);
    unmerged_an_edge = false;

    if ismember(edgeName, {'all', 'top'})
        if isempty(buff_arrays{buff_n_top})
            if mergedTop
                fprintf("Tile '%s(%d)' (%s side buffer) is empty; if '%s' is merged, cannot undo\n",...
                    buff_array_name, buff_n_top, 'top', data_array_name);
                unmerge_failure_top = true;
            end
        elseif mergedTop || ignoreMergedVars
            fprintf("Undoing '%s' merge on top edge\n", data_array_name);
            [nrows,~] = size(buff_arrays{buff_n_top});
            data_array(1:nrows,:) = buff_arrays{buff_n_top};
            unmerged_an_edge = true;
        end
    end

    if ismember(edgeName, {'all', 'bottom'})
        if isempty(buff_arrays{buff_n_bottom})
            if mergedBottom
                fprintf("Tile '%s(%d)' (%s side buffer) is empty; if '%s' is merged, cannot undo\n",...
                    buff_array_name, buff_n_bottom, 'bottom', data_array_name);
                unmerge_failure_bottom = true;
            end
        elseif mergedBottom || ignoreMergedVars
            fprintf("Undoing '%s' merge on bottom edge\n", data_array_name);
            [nrows,~] = size(buff_arrays{buff_n_bottom});
            data_array((end-nrows+1):end,:) = buff_arrays{buff_n_bottom};
            unmerged_an_edge = true;
        end
    end

    if ismember(edgeName, {'all', 'left'})
        if isempty(buff_arrays{buff_n_left})
            if mergedLeft
                fprintf("Tile '%s(%d)' (%s side buffer) is empty; if '%s' is merged, cannot undo\n",...
                    buff_array_name, buff_n_left, 'left', data_array_name);
                unmerge_failure_left = true;
            end
        elseif mergedLeft || ignoreMergedVars
            fprintf("Undoing '%s' merge on left edge\n", data_array_name);
            [~,ncols] = size(buff_arrays{buff_n_left});
            data_array(:,1:ncols) = buff_arrays{buff_n_left};
            unmerged_an_edge = true;
        end
    end

    if ismember(edgeName, {'all', 'right'})
        if isempty(buff_arrays{buff_n_right})
            if mergedRight
                fprintf("Tile '%s(%d)' (%s side buffer) is empty; if '%s' is merged, cannot undo\n",...
                    buff_array_name, buff_n_right, 'right', data_array_name);
                unmerge_failure_right = true;
            end
        elseif mergedRight || ignoreMergedVars
            fprintf("Undoing '%s' merge on right edge\n", data_array_name);
            [~,ncols] = size(buff_arrays{buff_n_right});
            data_array(:,(end-ncols+1):end) = buff_arrays{buff_n_right};
            unmerged_an_edge = true;
        end
    end

    if unmerged_an_edge
        setfield(m_struct, data_array_name, data_array);
        unmerged_array_names{end+1, 1} = data_array_name;
    end
end

if ismember(edgeName, {'all', 'top'})
    if ismember('mergedTop', m_varlist) && ~unmerge_failure_top
        m_struct.mergedTop = false;
    end
end
if ismember(edgeName, {'all', 'bottom'}) && ~unmerge_failure_bottom
    if ismember('mergedBottom', m_varlist)
        m_struct.mergedBottom = false;
    end
end
if ismember(edgeName, {'all', 'left'}) && ~unmerge_failure_left
    if ismember('mergedLeft', m_varlist)
        m_struct.mergedLeft = false;
    end
end
if ismember(edgeName, {'all', 'right'}) && ~unmerge_failure_right
    if ismember('mergedRight', m_varlist)
        m_struct.mergedRight = false;
    end
end

if ~isempty(unmerged_array_names)
    disp("Undid merge on the following data arrays:")
    disp(unmerged_array_names)
end


%% Present state of tile struct

m_varlist = who(m_struct);

mergedTop = ismember('mergedTop', m_varlist) && m_struct.mergedTop;
mergedBottom = ismember('mergedBottom', m_varlist) && m_struct.mergedBottom;
mergedLeft = ismember('mergedLeft', m_varlist) && m_struct.mergedLeft;
mergedRight = ismember('mergedRight', m_varlist) && m_struct.mergedRight;

fprintf("\n");
fprintf("State of tile after undo:\n");
fprintf("mergedTop = %s\n", mat2str(mergedTop));
fprintf("mergedBottom = %s\n", mat2str(mergedBottom));
fprintf("mergedLeft = %s\n", mat2str(mergedLeft));
fprintf("mergedRight = %s\n", mat2str(mergedRight));
fprintf("\n");
