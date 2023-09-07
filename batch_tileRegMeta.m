%batch_tileMeta

udir='/mnt/pgc/data/elev/dem/setsm/ArcticDEM/mosaic/2m_v1.0';
dataset='GLA14_rel634';
f=dir([udir,'/1*']);
for i=1:length(f);
    f1=f(i);
    if exist([udir,'/',f1.name],'dir') == 7
        g=dir([udir,'/',f1.name,'/*2m_reg_dem.mat']);
        gp=[udir,'/',f1.name,'/',g.name];
        if ~isdir(gp)
            fprintf('tileRegMeta %s\n',gp);
            tileRegMeta(gp,dataset);
        end
    end
end
