function convertQcDataToTxt(txtfile)
%loadQcData - Loads QC.mat file data from folders with DEMs in the meta db
%   into the meta db structure

udir='/mnt/pgc/data/elev/dem/setsm/ArcticDEM/region';
rdir=dir([udir,'/region_*']);
for i=1:length(rdir);
    regionDir{i}=[udir,'/',rdir(i).name,'/strips/2m'];
end

qc.fileNames=[];
qc.flag=[];
qc.x =[];
qc.y=[];

for i=1:length(regionDir)

    if exist([regionDir{i},'/qc.mat'],'file')

        qci=load([regionDir{i},'/qc.mat']);

        qc.fileNames=[qc.fileNames(:);qci.fileNames(:)];
        qc.flag=[qc.flag(:);qci.flag(:)];
        qc.x=[qc.x(:);qci.x(:)];
        qc.y=[qc.y(:);qci.y(:)];

        clear qci;
    end

end

fprintf('QC records: %d\n',length(qc.fileNames));

if ~isempty(qc.fileNames)
    
    f = regexprep(qc.fileNames, '[\w\d\\:]+/','');
    f=strrep(f,'_dem_browse.tif','_v3.0');
    f=strcat('SETSM_',f);
    
    flag = qc.flag;
    
    clear qc
    
    fid=fopen(txtfile,'w');
    for i=1:length(f)
        
        fprintf(fid,'%s,%d\n',f{i},flag(i));
        
    end
    
end


