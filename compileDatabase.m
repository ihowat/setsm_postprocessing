% built a searchable index of all data

udir='/mnt/pgc/data/elev/dem/setsm/ArcticDEM/region';
rdir=dir([udir,'/arcticdem_3*']);
outname='/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing_pgc/arcticDEMdatabase_2m_20191125.mat';

k=0;
for i=1:length(rdir);
    fprintf('Building array of DEM paths in %s\n',rdir(i).name);
    pdir=dir([udir,'/',rdir(i).name,'/strips_unf/2m/*_2m_lsf']);
    for j=1:length(pdir);
        k=k+1;
        demDir{k}=[pdir(j).folder,'/',pdir(j).name];
    end
end

f               = cell(length(rdir),1);
creationDate    = cell(length(rdir),1);
stripDate       = cell(length(rdir),1);
res             = cell(length(rdir),1);
x               = cell(length(rdir),1);
y               = cell(length(rdir),1);
avg_rmse        = cell(length(rdir),1);
med_rmse        = cell(length(rdir),1);
max_rmse        = cell(length(rdir),1);
Nscenes         = cell(length(rdir),1);

datasets        = cell(length(rdir),1);
Ngcps           = cell(length(rdir),1);
NgcpsRock       = cell(length(rdir),1);
NgcpsIce        = cell(length(rdir),1);
quads           = cell(length(rdir),1);
avgIceGCPdt     = cell(length(rdir),1);
trans           = cell(length(rdir),1);
avgdz           = cell(length(rdir),1);
avgdz_rock      = cell(length(rdir),1);
avgdz_ice       = cell(length(rdir),1);
stddz           = cell(length(rdir),1);
stddz_rock      = cell(length(rdir),1);
stddz_ice       = cell(length(rdir),1);
dzpc            = cell(length(rdir),1);
dzpc_rock       = cell(length(rdir),1);
dzpc_ice        = cell(length(rdir),1);

i=1;
for i=1:length(demDir);

    fprintf('reading %s\n',demDir{i})

    filecheck=dir([demDir{i},'/*meta.txt']);

    if isempty(filecheck); continue; end

    [f{i},creationDate{i},stripDate{i},res{i},x{i},y{i},...
        avg_rmse{i},med_rmse{i},max_rmse{i},Nscenes{i},reg...
    ]=compileStripMeta(demDir{i});

    f{i} = cellstr([repmat([demDir{i},'/'],length(f{i}),1),char(f{i})]);


       datasets{i} = reg.datasets;
          Ngcps{i} = reg.Ngcps;
      NgcpsRock{i} = reg.NgcpsRock;
       NgcpsIce{i} = reg.NgcpsIce;
          quads{i} = reg.quads;
    avgIceGCPdt{i} = reg.avgIceGCPdt;
          trans{i} = reg.trans;
          avgdz{i} = reg.avgdz;
     avgdz_rock{i} = reg.avgdz_rock;
      avgdz_ice{i} = reg.avgdz_ice;
          stddz{i} = reg.stddz;
     stddz_rock{i} = reg.stddz_rock;
      stddz_ice{i} = reg.stddz_ice;
           dzpc{i} = reg.dzpc;
      dzpc_rock{i} = reg.dzpc_rock;
       dzpc_ice{i} = reg.dzpc_ice;

end

a.f=[]; i=1; for i=1:length(f); a.f = [a.f;f{i}]; end; clear f;

a.region=cell(size(a.f));
i=1;
for i=1:length(a.f); 
        a.region{i} = fileparts(a.f{i});
end

a.creationDate=cell2mat(creationDate(:)); clear creationDate
a.stripDate = cell2mat(stripDate(:)); clear stripDate
a.res= cell2mat(res(:)); clear res

a.x=[]; i=1; for i=1:length(x); a.x = [a.x;x{i}]; end; clear x;
a.y=[]; i=1; for i=1:length(y); a.y = [a.y;y{i}]; end; clear y;

a.avg_rmse= cell2mat(avg_rmse(:)); clear avg_rmse
a.med_rmse= cell2mat(med_rmse(:)); clear med_rmse
a.max_rmse= cell2mat(max_rmse(:)); clear max_rms
a.Nscenes= cell2mat(Nscenes(:)); clear Nscenes

a.datasets=[]; i=1; for i=1:length(datasets); 
    a.datasets = [a.datasets;datasets{i}]; end; clear datasets;

a.Ngcps= cell2mat(Ngcps(:)); clear Ngcps;
a.NgcpsRock= cell2mat(NgcpsRock(:)); clear NgcpsRock;
a.NgcpsIce= cell2mat(NgcpsIce(:)); clear NgcpsIce;
a.quads= cell2mat(quads(:)); clear quads;
a.avgIceGCPdt= cell2mat(avgIceGCPdt(:)); clear avgIceGCPdt;
a.trans= cell2mat(trans(:)); clear trans;
a.avgdz= cell2mat(avgdz(:)); clear avgdz;
a.avgdz_rock= cell2mat(avgdz_rock(:)); clear avgdz_rock;
a.avgdz_ice= cell2mat(avgdz_ice(:)); clear avgdz_ice;
a.stddz= cell2mat(stddz(:)); clear stddz;
a.stddz_rock= cell2mat(stddz_rock(:)); clear stddz_rock;
a.stddz_ice= cell2mat(stddz_ice(:)); clear stddz_ice;
a.dzpc= cell2mat(dzpc(:)); clear dzpc;
a.dzpc_rock= cell2mat(dzpc_rock(:)); clear dzpc_rock;
a.dzpc_ice= cell2mat(dzpc_ice(:)); clear dzpc_ice;

a.xmin = cellfun(@min, a.x);
a.xmax=cellfun(@max, a.x);
a.ymin = cellfun(@min, a.y);
a.ymax=cellfun(@max, a.y);

save(outname,'-struct','a','-v7.3');
