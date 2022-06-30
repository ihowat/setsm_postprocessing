
is2_dir  = '/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/from_unity/altimetryByTile';


is2_file_listing = dir([is2_dir,'/*_is2.mat']);
is2_files_cellarr = fullfile({is2_file_listing.folder}.', {is2_file_listing.name}.');

for file_idx = 1:length(is2_files_cellarr)
    is2_matfile = is2_files_cellarr{file_idx};
    is2_csvfile = strrep(is2_matfile, '.mat', '.csv');
    if strcmp(is2_matfile, is2_csvfile)
        error('Could not determine output CSV filename');
    end

    fprintf('Loading is2 matfile: %s\n', is2_matfile);
    is2 = load(is2_matfile);

    fprintf('Writing is2 csvfile: %s\n', is2_csvfile);
    is2_csvfile_fid = fopen(is2_csvfile, 'w');
    fprintf(is2_csvfile_fid, 'x,y,z,t\n');
    for i = 1:length(is2.x)
        fprintf(is2_csvfile_fid, '%f,%f,%f,%f\n', is2.x(i,1), is2.y(i,1), is2.z(i,1), is2.t(i,1));
    end
    fclose(is2_csvfile_fid);
end
