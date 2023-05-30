
shp_dbase_path = '/mnt/pgc/data/projects/earthdem/strip_databases/unity_databases/arcticdem_v4.1_greenland/arcticdem_strips.shp';
mat_dbase_path = '/mnt/pgc/data/projects/earthdem/strip_databases/ArcticDEMdatabase4_2m_v4.1_20230425_all_plus_reproj.mat';
mat_dbase_cropped_path = '/mnt/pgc/data/projects/earthdem/strip_databases/ArcticDEMdatabase4_2m_v4.1_20230425_clip_iangreenland.mat';

fprintf('Loading shapefile dbase: %s\n', shp_dbase_path);
shp_dbase = shaperead(shp_dbase_path);

shp_strip = arrayfun(@(feat) feat.strip, shp_dbase, 'UniformOutput',false);
shp_strip_unique = unique(shp_strip);
clear shp_dbase shp_strip;

fprintf('Loading matfile dbase: %s\n', mat_dbase_path);
mat_dbase = load(mat_dbase_path);

[stripdir,~,~] = cellfun(@(x) fileparts(x), mat_dbase.fileName, 'UniformOutput',false);
[~,stripdname,~] = cellfun(@(x) fileparts(x), stripdir, 'UniformOutput',false);
mat_strip = cellfun(@(x) x(1:54), stripdname, 'UniformOutput',false);

C = intersect(mat_strip, shp_strip_unique);
ia = find(ismember(mat_strip, C));
mat_dbase_clip = structfun(@(x) x(ia), mat_dbase, 'UniformOutput',false);

fprintf('Writing cropped matfile dbase: %s\n', mat_dbase_cropped_path);
save(mat_dbase_cropped_path,'-struct','mat_dbase_clip','-v7.3');
