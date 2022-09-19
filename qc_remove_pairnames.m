
pairnames_txt = '/mnt/pgc/data/elev/dem/setsm/ArcticDEM/region/arcticdem_greenland_update_few_pairnames.txt';

qc_files_to_check = {
%    '/mnt/pgc/data/elev/dem/setsm/ArcticDEM/region/arcticdem_05_greenland_northeast/strips_unf/2m/qc.mat',
    '/mnt/pgc/data/elev/dem/setsm/ArcticDEM/region/arcticdem_06_greenland_northwest/strips_unf/2m/qc.mat',
};


fileID = fopen(pairnames_txt, 'r');
cols = textscan(fileID, '%s', inf, 'delimiter', ',');
[reset_pnames] = cols{:};
fclose(fileID);

for i = 1:length(qc_files_to_check)
    qcFile = qc_files_to_check{i};

    fprintf('Loading qc file %s\n', qcFile);
    qc = load(qcFile);
    qc_fnames = cellfun(@(s) s(max([strfind(s, '\') strfind(s, '/')])+1:end), qc.fileNames, 'UniformOutput',false);
    qc_pnames = cellfun(@(s) s(1:47), qc_fnames, 'UniformOutput',false);

    [C,IA] = intersect(qc_pnames, reset_pnames);
    Lia = ismember(qc_pnames, reset_pnames);

    if isempty(C)
        fprintf("Found no matching pairnames in qc file to remove\n");
        continue
    else
        fprintf("Found %d strip segments in qc file matching %d pairnames to remove\n", nnz(Lia), length(C));
    end

    qc_remove_i = Lia;

    qc_fields = fieldnames(qc);
    for i_field = 1:numel(qc_fields)
        field = qc.(qc_fields{i_field});
        field(qc_remove_i,:) = [];
        qc.(qc_fields{i_field}) = field;
    end

    qcFile_bak = strrep(qcFile, 'qc.mat', 'qc_bak.mat');
    copyfile(qcFile, qcFile_bak);

    save(qcFile,'-struct','qc');

    clear qc;

end
