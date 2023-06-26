function call_insertSubset(baseTileFile, annualRootDir, boxShp, pShp, varargin)

if ~isa(baseTileFile, 'char') && ~isa(baseTileFile, 'string')
    error("baseTileFile should be a char/string path to a tile .mat file");
end
if ~isfile(baseTileFile)
    error("baseTileFile tile matfile does not exist: %s\n", baseTileFile);
end

[~,baseTileFname,baseTileFnameExt] = fileparts(baseTileFile);
baseTileFname = [baseTileFname,baseTileFnameExt];
baseSupertileName = get10mTileName(baseTileFname);
if isempty(baseSupertileName)
    error("Cannot parse supertile name from baseTileFile: %s", baseTileFile);
end
if ~startsWith(baseTileFname, baseSupertileName)
    error("baseTileFname tile filename should start with baseSupertileName ('%s'), but does not: %s", baseSupertileName, baseTileFile);
end

if ~isfolder(annualRootDir)
    error("annualRootDir folder does not exist: %s\n", annualRootDir);
end

boxCsv = strrep(boxShp, '.shp', '.csv');
if ~isfile(boxShp)
    error("boxShp file does not exist: %s\n", boxShp);
end
if ~isfile(boxCsv)
    error("boxShp partner CSV file does not exist: %s\n", boxCsv);
end

if ~isfile(pShp)
    error("pShp file does not exist: %s\n", pShp);
end

n = find(strcmpi(varargin,'overwrite'));
if ~isempty(n)
    overwrite = true;
else
    overwrite = false;
end


box_mapstruct = shaperead(boxShp);
box_polyshape_arr = arrayfun(@(feat) polyshape(feat.X, feat.Y), box_mapstruct);
box_csvarr = readtable(boxCsv);

p_mapstruct = shaperead(pShp);
p_polyshape_arr = arrayfun(@(feat) polyshape(feat.X, feat.Y), p_mapstruct);

box_mapstruct_supertile = {box_mapstruct.supertile};
box_mapstruct_patchid = cell2mat({box_mapstruct.patchid});
p_mapstruct_patchid = cell2mat({p_mapstruct.patchid});


all_valid = true;
for i = 1:length(box_mapstruct)
    box_ms = box_mapstruct(i);
    box_patchid = box_ms.patchid;
    valid = true;

    % Check insertyear only for patches that are applicable to this tile, which happens later.
    if false && isnan(box_ms.insertyear)
        fprintf(2, "ERROR: Box shapefile feature [%d] patchid [%d] missing 'insertyear'\n", i, box_patchid);
        valid = false;
    end

    valid_geom = true;
    if ~isnan(box_ms.xmin) && (box_ms.xmin ~= nanmin(box_ms.X))
        valid_geom = false;
    end
    if ~isnan(box_ms.xmax) && (box_ms.xmax ~= nanmax(box_ms.X))
        valid_geom = false;
    end
    if ~isnan(box_ms.ymin) && (box_ms.ymin ~= nanmin(box_ms.Y))
        valid_geom = false;
    end
    if ~isnan(box_ms.ymax) && (box_ms.ymax ~= nanmax(box_ms.Y))
        valid_geom = false;
    end
    if ~valid_geom
        fprintf(2, "ERROR: Box shapefile feature [%d] patchid [%d] geometry does not align with xmin/xmax/ymin/ymax field(s)\n", i, box_patchid);
        valid = false;
    end

    box_ms_p = box_ms.p;
    if ~strcmp(box_ms_p, '')
        box_csv_idx = find(box_csvarr.patchid == box_patchid);
        if isempty(box_csv_idx)
            fprintf(2, "ERROR: No Box CSV feature matches Box shapefile feature [%d] patchid [%d]\n", i, box_patchid);
            valid = false;
        elseif length(box_csv_idx) > 1
            fprintf(2, "ERROR: More than one Box CSV feature match Box shapefile feature [%d] patchid [%d]\n", i, box_patchid);
            valid = false;
        else
            box_csv_p = box_csvarr.p{box_csv_idx};
            if strcmp(box_csv_p, '')
                fprintf(2, "ERROR: Box CSV feature 'p' is empty while not empty in Box shapefile feature [%d] patchid [%d]\n", i, box_patchid);
                valid = false;
            else
                box_p = box_csv_p;

                p_match_idx = find(p_mapstruct_patchid == box_patchid);
                p_match_idx = p_match_idx(:);
                if isempty(p_match_idx)
                    fprintf(2, "ERROR: No P shapefile feature matches Box feature [%d] patchid [%d]\n", i, box_patchid);
                    valid = false;
                elseif length(p_match_idx) > 1
                    fprintf(2, "ERROR: More than one P shapefile feature match Box feature [%d] patchid [%d]\n", i, box_patchid);
                    valid = false;
                else
                    box_p_num = str2num(box_p);
                    p_ps = p_polyshape_arr(p_match_idx);
                    p_ps_vert = p_ps.Vertices;
                    if ~isempty(setxor(box_p_num, p_ps_vert))
                        fprintf(2, "ERROR: Box shapefile feature [%d] 'p' vertex coords string mismatch P feature [%d] patchid [%d] vertex coords\n", i, p_match_idx, box_patchid);
                        valid = false;
                    end
                end
            end
        end
    end

    if ~valid
        all_valid = false;
    end
end
if ~all_valid
    error("Box shapefile has issues")
end


boxes_match_supername_idx = find(strcmpi(baseSupertileName, box_mapstruct_supertile));
boxes_match_supername_idx = boxes_match_supername_idx(:);

boxes_no_supername_idx = find(strcmp('', box_mapstruct_supertile));
boxes_no_supername_idx = boxes_no_supername_idx(:);

m = matfile(baseTileFile);
xmin = m.x(1,1);
xmax = m.x(1,end);
ymin = m.y(1,end);
ymax = m.y(1,1);
clear m;

tile_vert = [[xmin, ymin];...
             [xmin, ymax];...
             [xmax, ymax];...
             [xmax, ymin]];
tile_ps = polyshape(tile_vert(:,1), tile_vert(:,2));

box_overlaps = overlaps(tile_ps, box_polyshape_arr);
boxes_intersect_idx = find(box_overlaps);
boxes_intersect_idx = boxes_intersect_idx(:);

boxes_for_tile_idx = union(boxes_match_supername_idx, boxes_no_supername_idx);
boxes_for_tile_idx = intersect(boxes_for_tile_idx, boxes_intersect_idx);
boxes_for_tile_dropped_idx = setdiff(boxes_match_supername_idx, boxes_for_tile_idx);

boxes_for_tile_pid = box_mapstruct_patchid(boxes_for_tile_idx);
[~,~,p_for_tile_idx] = intersect(boxes_for_tile_pid, p_mapstruct_patchid);

if isempty(boxes_for_tile_idx)
    fprintf("No patches to apply\n");
    return;
end

% Order patches by patchid
boxes_for_tile_patchid = cell2mat({box_mapstruct(boxes_for_tile_idx).patchid});
[~,sort_I] = sort(boxes_for_tile_patchid);
boxes_for_tile_idx = boxes_for_tile_idx(sort_I);


patchTileFile = strrep(baseTileFile,'.mat','_patched.mat');
patchTempFile = strrep(baseTileFile,'.mat','_patched.mat.temp');

all_valid = true;
for i = 1:length(boxes_for_tile_idx)
    box_idx = boxes_for_tile_idx(i);
    box_ms = box_mapstruct(box_idx);
    box_ps = box_polyshape_arr(box_idx);
    box_patchid = box_ms.patchid;

    box_insertyear = box_ms.insertyear;
    if isnan(box_insertyear)
        fprintf(2, "ERROR: Box shapefile feature [%d] patchid [%d] missing 'insertyear'\n", box_idx, box_patchid);
        all_valid = false;
    end

    insertFile = fullfile(annualRootDir, num2str(box_insertyear), baseSupertileName, baseTileFname);
    if ~isfile(insertFile)
        fprintf(2, "ERROR: insertFile for Box feature [%d] patchid [%d] does not exist: %s\n", box_idx, box_patchid, insertFile);
        all_valid = false;
    end
end
if ~all_valid
    error("Box shapefile has issues")
end

%if isfile(patchTempFile)
%    fprintf("Removing existing patched.mat.temp file: %s\n", patchTempFile);
%    delete(patchTempFile);
%end
if isfile(patchTileFile)
    if overwrite
        fprintf("Removing existing patched.mat file (overwrite option provided): %s\n", patchTileFile);
        delete(patchTileFile);
    else
        fprintf("Skipping existing patched.mat file (overwrite option not provided): %s\n", patchTileFile);
        return;
    end
end

if ~isfile(patchTempFile)
    fprintf("Copying base tile .mat file to new patched.mat.tmp file: %s\n", patchTempFile);
    eval(['!cp ',baseTileFile,' ',patchTempFile]);
end

for i = 1:length(boxes_for_tile_idx)
    box_idx = boxes_for_tile_idx(i);
    box_ms = box_mapstruct(box_idx);
    box_ps = box_polyshape_arr(box_idx);
    box_patchid = box_ms.patchid;

    box_insertyear = box_ms.insertyear;
    if isnan(box_insertyear)
        error("Box shapefile feature [%d] patchid [%d] missing 'insertyear'", box_idx, box_patchid);
        all_valid = false;
    end

    insertFile = fullfile(annualRootDir, num2str(box_insertyear), baseSupertileName, baseTileFname);
    if ~isfile(insertFile)
        error("insertFile for Box feature [%d] patchid [%d] does not exist: %s", box_idx, box_patchid, insertFile);
    end

    box_vert = box_ps.Vertices;
    xr = [nanmin(box_vert(:,1)), nanmax(box_vert(:,1))];
    yr = [nanmin(box_vert(:,2)), nanmax(box_vert(:,2))];

    p_idx = find(p_mapstruct_patchid == box_patchid);
    if isempty(p_idx)
        p = box_ps.Vertices;
    else
        p_ps = p_polyshape_arr(p_idx);
        p = p_ps.Vertices;
    end

    insertSubset(patchTempFile, insertFile, xr, yr, p, patchTempFile);
%    patchTempFile
%    insertFile
%    xr
%    yr
%    p
end

fprintf("Renaming patched.mat.tmp file to patched.mat: %s\n", patchTileFile);
eval(['!mv ',patchTempFile,' ',patchTileFile]);
