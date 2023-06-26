function dumpInsertSubsetInfo(tile_mat)

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
varlist0 = who(m0);

baseFile = m0.baseFile;
insertFile = m0.insertFile;
if isa(insertFile, 'char')
    wrap_in_cell = true;
else
    wrap_in_cell = false;
end

xr = m0.xr;
yr = m0.yr;
if ismember('p', varlist0)
    p = m0.p;
else
    p = [];
end

if wrap_in_cell
    insertFile = {insertFile};
    xr = {xr};
    yr = {yr};
end

[~,baseFile,~] = fileparts(baseFile);

%fileID = fopen('/mnt/pgc/data/elev/dem/setsm/ArcticDEM/mosaic/v4.1/from_unity/greenland_10m_grimp_results_insert-subset-dump.txt', 'w');

for i = 1:length(insertFile)
    try
        i_insertFile = insertFile{i};
        [~,i_insertFile,~] = fileparts(i_insertFile);
        i_xr = xr{i};
        i_yr = yr{i};
        if ~isempty(p)
            i_p = p{i};
            i_p = mat2str(i_p);
        else
            i_p = '';
        end
        fprintf("%s,%s,%d,%d,%d,%d,%s\n", baseFile, i_insertFile, i_xr(1), i_xr(2), i_yr(1), i_yr(2), i_p);
%        fprintf(fileID, "%s,%s,%d,%d,%d,%d,%s\n", baseFile, i_insertFile, i_xr(1), i_xr(2), i_yr(1), i_yr(2), i_p);
    catch ME
        tile_mat
        i_xr
        i_yr
        i_p
        keyboard
        rethrow(ME);
    end
end

%fclose(fileID);
