
dbase_in    = '/mnt/pgc/data/projects/earthdem/strip_databases/ArcticDEMdatabase4_2m_v4.1_20230425_all_plus_reproj.mat';
matfile_out = '/mnt/pgc/data/projects/earthdem/strip_databases/ArcticDEMdatabase4_2m_v4.1_20230425_all_plus_reproj_remerge_info.mat';

fprintf('Loading strip database: %s\n', dbase_in);
meta = matfile(dbase_in);

[stripID,remerge_info_file] = cellfun(@(x) stripFileNameToParts(x), meta.fileName, 'UniformOutput',false);
[~,ia,~] = unique(stripID);
stripID = stripID(ia);
remerge_info_file = remerge_info_file(ia);

out.stripID = stripID;
out.remerge_info_file = remerge_info_file;

fprintf('Writing %d total records to %s\n', length(out.stripID), matfile_out);
save(matfile_out,'-struct','out','-v7.3');


function [stripID,remerge_info_file] = stripFileNameToParts(fileName)
    [stripdir,~,~] = fileparts(fileName);
    [~,stripdname,~] = fileparts(stripdir);
    stripID_parts = split(stripdname, '_');
    stripID = strjoin(stripID_parts(1:4), '_');
    remerge_info_file = [stripdir,'/',stripdname,'_remerge.info'];
end
