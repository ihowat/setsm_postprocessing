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
    
    % get indices where flag = 3
    I = qc.flag == 3;
    % get basenames, x, y for those indicies
    f = regexprep(qc.fileNames(I), '[\w\d\\:]+/','');
    f=strrep(f,'_dem_browse.tif','_v2.0');
    f=strcat('SETSM_',f);
    
    x=qc.x(I);
    y=qc.y(I);
    
    clear qc
    
    fid=fopen(txtfile,'w');
    for i=1:length(x)
        coords={};
        for j=1:length(x{i})
            cx=x{i}{j};
            cy=y{i}{j};
            if ~(isempty(cx) || isempty(cy))
                c=mat2str([cx,cy]);
                c2=strrep(strrep(strrep(c,'[','(('),']','))'),';',',');
                coords{end+1}=c2;
            end
        end
        if ~isempty(coords)
            c3=strjoin(coords,',');
            fprintf(fid,'%s;MULTIPOLYGON(%s)\n',f{i},c3);
        end
    end
    
end


