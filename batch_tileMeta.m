%batch_tileMeta

udir='/mnt/pgc/data/elev/dem/setsm/ArcticDEM/mosaic/2m_v1.3';
f=dir([udir,'/*']);
for i=1:length(f);
    f1=f(i);
    if exist([udir,'/',f1.name],'dir') == 7
        g=dir([udir,'/',f1.name,'/*2m_dem.mat']);
        gp=[udir,'/',f1.name,'/',g.name];
        if ~isdir(gp)
            fprintf('tileMeta %s\n',gp);
            tileMeta(gp);
        end
    end
end
