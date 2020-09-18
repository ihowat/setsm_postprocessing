
%[~,homeDir]=system('echo ~');
%homeDir=homeDir(1:end-1);

res=2;
%outname =[homeDir,'/data4/REMA/arcticDEMdatabase_',num2str(res),'m.mat'];
outname='earthdem_database_unf.mat';

%upDir{1}=[homeDir,'/fs/project/howat.4/EarthDEM'];
%regionDirs{1}=dir([upDir{1},'/*']);
%regionDirs{1}=cellfun(@(x) [upDir{1},'/',x,'/strips_unf/2m'], {regionDirs{1}.name},...
%    'UniformOutput',false);
%upDir{2}=[homeDir,'/data5/REMA'];
%regionDirs{2}=dir([upDir{2},'/region_*']);
%regionDirs{2}=cellfun(@(x) [upDir{2},'/',x,'/strips_unf/2m'], {regionDirs{2}.name},...
 %   'UniformOutput',false);

%regionDirs = cat(2,regionDirs{:});

upDir=['/fs/project/howat.4/EarthDEM'];
regionDirs=dir([upDir,'/region*']);
regionDirs=cellfun(@(x) [upDir,'/',x,'/strips_unf/2m'], {regionDirs([regionDirs.isdir]).name},...
    'UniformOutput',false);

meta=[];

if exist(outname,'file')
    out0=load(outname);
end

i=1;
for i=1:length(regionDirs)
    
     regionDir=regionDirs{i};
    
     if exist(regionDir,'dir')
         
        metaFiles=dir([regionDir,'/*/*meta.txt']);
        metaFiles = strcat({metaFiles.folder}',repmat({'/'},length(metaFiles),1),{metaFiles.name}');
         
        if exist('out0','var')
            [~,IA] = intersect(metaFiles, out0.fileName);
            metaFiles(IA) = [];
            if isempty(metaFiles)
                continue
            end
        end
        
         j=1;
         for j=1:length(metaFiles)
             metaFile=metaFiles{j};
             fprintf('adding file %s\n',metaFile)
             if isempty(meta)
                 meta=readStripMeta(metaFile,'noSceneMeta');
             else
                 meta(length(meta)+1)=readStripMeta(metaFile,'noSceneMeta');
             end
             
         end 
     end
end

flds = fields(meta);
i=1;
for i =1:length(flds)
    fld = flds{i};
    eval(['out.',fld,'= {meta.',fld,'};']);
end
    
%out.fileName = strrep(out.fileName,homeDir,'');

[filePath,out.stripName] = cellfun(@fileparts,out.fileName,'uniformoutput',0);
out.stripName=strrep(out.stripName,'_meta','');

stripNameChar=char(out.stripName{:});

out.stripDate=datenum(stripNameChar(:,6:13),'yyyymmdd');
out.stripDate=out.stripDate(:)';

out.satID=cellstr(stripNameChar(:,1:4))';

out.creation_date = [out.creation_date{:}];
out.strip_creation_date = [out.strip_creation_date{:}];

fprintf('writing %d new records to database file\n',length(out.fileName));
if exist('out0','var')
    flds = fields(out);
    i=1;
    for i =1:length(flds)
        fld = flds{i};
        eval(['out.',fld,'= [out0.',fld,',out.',fld,'];']);
    end
end

fprintf('writing %d records to %s\n',length(out.fileName),outname);
save(outname,'-struct','out','-v7.3');
