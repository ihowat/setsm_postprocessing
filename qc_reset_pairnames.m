
pairnames_txt = '/mnt/pgc/data/elev/dem/setsm/ArcticDEM/region/arcticdem_greenland_maxh_rerun_pairnames.txt';

sdir1 = '/mnt/pgc/data/elev/dem/setsm/ArcticDEM/region/arcticdem_05_greenland_northeast/strips_unf/2m';
sdir2 = '/mnt/pgc/data/elev/dem/setsm/ArcticDEM/region/arcticdem_06_greenland_northwest/strips_unf/2m';

sdir1_win = 'V:\pgc\data\elev\dem\setsm\ArcticDEM\region\arcticdem_05_greenland_northeast\strips_unf\2m';
sdir2_win = 'V:\pgc\data\elev\dem\setsm\ArcticDEM\region\arcticdem_06_greenland_northwest\strips_unf\2m';

stripdirs = [string(sdir1) string(sdir2)];
stripdirs_win = [string(sdir1_win) string(sdir2_win)];


fileID = fopen(pairnames_txt, 'r');
cols = textscan(fileID, '%s', inf, 'delimiter', ',');
[reset_pnames] = cols{:};
fclose(fileID);


for i_sdir = 1:length(stripdirs)
    sdir = char(stripdirs(i_sdir));
    sdir_win = char(stripdirs_win(i_sdir));
    
    qc_file = [sdir,'/qc.mat'];
    fprintf('Loading qc file %s\n', qc_file);
    qc = load(qc_file);
    qc_fnames = cellfun(@(s) s(max([strfind(s, '\') strfind(s, '/')])+1:end), qc.fileNames, 'UniformOutput',false);
    qc_pnames = cellfun(@(s) s(1:47), qc_fnames, 'UniformOutput',false);
    
    fewer_segments = 0; 
    
    for i_pname = 1:length(reset_pnames)
        reset_pname = reset_pnames{i_pname};
        fprintf('%s',reset_pname);
        
        sdir_pname_files = dir([sdir,'/',reset_pname,'*/*_dem.tif']);
        if ~isempty(sdir_pname_files)
            fprintf(', %d in strip folder', numel(sdir_pname_files));
        end
        sdir_pname_fnames = {sdir_pname_files.name}.';
        sdir_pname_fnames = cellfun(@(s) strrep(s,'dem.tif','dem_10m_shade_masked.tif'), sdir_pname_fnames, 'UniformOutput',false);
        
        qc_reset_i = find(strcmp(reset_pname, qc_pnames));
        if ~isempty(qc_reset_i)
            fprintf(', %d in qc.mat', numel(qc_reset_i));
        end
        
        if isempty(sdir_pname_files) && isempty(qc_reset_i)
            fprintf(', skipping\n');
            continue;
        elseif isempty(sdir_pname_files) || isempty(qc_reset_i)
            fprintf(', error? skipping\n');
            continue;
        end
        
        fprintf('\n');
        
        fewer_segments = fewer_segments + (numel(qc_reset_i) - numel(sdir_pname_files));
        
        qc_fields = fieldnames(qc);
        for i_field = 1:numel(qc_fields)
            field = qc.(qc_fields{i_field});
            field(qc_reset_i,:) = [];
            qc.(qc_fields{i_field}) = field;
        end
        
        for i_fname = 1:numel(sdir_pname_fnames)
            qc.fileNames = [qc.fileNames; [sdir_win,'/',sdir_pname_fnames{i_fname}]];
            qc.flag = [qc.flag; 0];
            qc.x = [qc.x; {cell(1)}];
            qc.y = [qc.y; {cell(1)}];
        end
        
    end
    
    fprintf('qc.mat file will now contain %d fewer strip segments\n', fewer_segments);
    
%    save(qc_file,'-struct','qc');
    
end

clear
