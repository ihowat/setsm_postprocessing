function resetMergeTileBuffers(tile_mat)

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

m0 = m_struct;
m0.Properties.Writable = true;

varlist0 = who(m0);

data_array_names = {
    'z',
    'z_mad',
    'N',
    'Nmt',
    'tmin',
    'tmax',
};
if ismember('waterFillMask', varlist0)
    data_array_names{end+1} = 'waterFillMask';
end

if ~all(ismember(data_array_names, varlist0))
    disp(data_array_names);
    error("One or more expected data arrays do not exist in tile struct")
end

buff_edge_names = cellfun(@(x) [x,'buff'], data_array_names, 'UniformOutput',false);
buff_corn_names = cellfun(@(x) [x,'buffcorners'], data_array_names, 'UniformOutput',false);

for array_idx = 1:length(buff_edge_names)
    buff_array_name = buff_edge_names{array_idx};
    if ismember(buff_array_name, varlist0)
        buff_arrays = getfield(m0, buff_array_name);
        for i = 1:length(buff_arrays)
            buff_arrays{i} = [];
        end
        setfield(m0, buff_array_name, buff_arrays);
    end
end

for array_idx = 1:length(buff_corn_names)
    buff_array_name = buff_corn_names{array_idx};
    if ismember(buff_array_name, varlist0)
        buff_arrays = getfield(m0, buff_array_name);
        for i = 1:length(buff_arrays)
            buff_arrays{i} = [];
        end
        setfield(m0, buff_array_name, buff_arrays);
    end
end

if ismember('origCornersList', varlist0)
    for i = 1:length(buff_edge_names)
        buff_edge_name = buff_edge_names{i};
        if ismember(buff_edge_name, fields(m0.origCornersList))
            setfield(m0.origCornersList, buff_edge_name, struct);
        end
    end
end

edge_list = {
    'left',
    'right',
    'top',
    'bottom'
};
for edge_idx = 1:length(edge_list)
    edge = edge_list{edge_idx};
    mergedVar = ['merged', upper(edge(1)), edge(2:end)];
    if ismember(mergedVar, varlist0)
        setfield(m0, mergedVar, false);
    end
    mergedVar = [edge,'NeedsRemerge'];
    if ismember(mergedVar, varlist0)
        setfield(m0, mergedVar, false);
    end
end
