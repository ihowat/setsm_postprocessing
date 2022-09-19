% built a searchable index of all data

addpath ../setsm_postprocessing3/

udir='/mnt/pgc/data/elev/dem/setsm/ArcticDEM/region';
rdir=dir([udir,'/arcticdem_*']);
outname_existing='/mnt/pgc/data/projects/earthdem/strip_databases/arcticDEMdatabase3_2m_unf_20200812.mat';
%outname_appended=outname_existing;
outname_appended='/mnt/pgc/data/projects/earthdem/strip_databases/arcticDEMdatabase3_2m_unf_greenland.mat';

%%% CHECK THIS SETTING %%%
report_number_of_strips_to_append_but_dont_actually_append = false;
%%% CHECK THIS SETTING %%%

b = {};
demDname_exist = {};
if exist('outname_existing', 'var') && exist(outname_existing, 'file') == 2
    fprintf('Loading existing database file: %s\n', outname_existing)
    b = load(outname_existing);
    [demDir_exist,~,~] = cellfun(@fileparts, b.f, 'UniformOutput', false);
    [~,demDname_exist,~] = cellfun(@fileparts, demDir_exist, 'UniformOutput', false);
    demDname_exist = unique(demDname_exist);
end

demDir={};
k=0;
for i=1:length(rdir);
    fprintf('Building array of DEM paths in %s ... ', rdir(i).name);
    pdir=dir([udir,'/',rdir(i).name,'/strips_unf/2m/*_2m_lsf*']);
    new_strip_num=0;
    for j=1:length(pdir);

        strippair_dir=[pdir(j).folder,'/',pdir(j).name];

        finfilecheck=dir([strippair_dir,'/*.fin']);
        if isempty(finfilecheck); continue; end

        demfilecheck=dir([strippair_dir,'/*dem.tif']);
        if isempty(demfilecheck); continue; end

        %% don't consider this strip if it's already in existing database
        if ~isempty(demDname_exist) && any(strcmp(demDname_exist, pdir(j).name)); continue; end

        new_strip_num=new_strip_num+1;
        k=k+1;
        demDir{k}=strippair_dir;

    end
    fprintf('%d strips to add\n', new_strip_num)
end

if report_number_of_strips_to_append_but_dont_actually_append
    return
end

f               = cell(length(rdir),1);
region          = cell(length(rdir),1);
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
last_print_len=0;
for i=1:length(demDir);

    fprintf(repmat('\b', 1, last_print_len));
    last_print_len=fprintf('reading (%d/%d): %s',i,k,demDir{i});

    filecheck=dir([demDir{i},'/*meta.txt']);
    finfilecheck=dir([demDir{i},'/*.fin']);

    if isempty(filecheck) || isempty(finfilecheck); continue; end

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
fprintf('\n')

a.f=[]; i=1; for i=1:length(f); a.f = [a.f;f{i}]; end; clear f;

a.region=cell(size(a.f));
i=1;
for i=1:length(a.f); 
        a.region{i} = fileparts(fileparts(a.f{i}));
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

%% merge existing and new database sets
if ~isempty(b)
    f = fieldnames(b);
    for i = 1:length(f)
        b.(f{i}) = [b.(f{i}); a.(f{i})];
    end
    a = b;
end

save(outname_appended,'-struct','a','-v7.3');
