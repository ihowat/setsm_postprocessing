%batch_tileMeta

addpath('/mnt/pgc/data/scratch/erik/repos/setsm_postprocessing4');

tiledir='/mnt/pgc/data/elev/dem/setsm/EarthDEM/mosaic/v1/results/output_tiles';
project='EarthDEM';
tileVersion='1.0';

tiledlist=dir([tiledir,'/*']);
for i=1:length(tiledlist);
    tilefolder=tiledlist(i);
    tilefolder_path=[tilefolder.folder,'/',tilefolder.name];
    if isdir(tilefolder_path)
        tilematflist=dir([tilefolder_path,'/*2m.mat']);
        for j=1:length(tilematflist)
            tilematfile=tilematflist(j);
            tilematfile=[tilematfile.folder,'/',tilematfile.name];
            tilemetafile=strrep(tilematfile, '2m.mat', '2m_meta.txt');
            if isfile(tilemetafile)
                fprintf("removing metafile: %s\n", tilemetafile);
                delete(tilemetafile);
            end
            fprintf("tileMetav4('%s','project','%s','tileVersion','%s')\n", tilematfile, project, tileVersion);
            tileMetav4(tilematfile,'project',project,'tileVersion',tileVersion);
        end
    end
end
