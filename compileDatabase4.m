
addpath ../setsm_postprocessing4/
if exist('out0','var')
    clear out0
end

%[~,homeDir]=system('echo ~');
%homeDir=homeDir(1:end-1);

res=2;
%outname =[homeDir,'/data4/REMA/arcticDEMdatabase_',num2str(res),'m.mat'];
outname='/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing_pgc/arcticDEMdatabase4_2m_unf_20200519.mat';
%outname_appended=outname;
outname_appended='/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing_pgc/arcticDEMdatabase4_2m_unf_20200519.mat';

%upDir{1}=[homeDir,'/fs/project/howat.4/EarthDEM'];
%regionDirs{1}=dir([upDir{1},'/*']);
%regionDirs{1}=cellfun(@(x) [upDir{1},'/',x,'/strips_unf/2m'], {regionDirs{1}.name},...
%    'UniformOutput',false);
%upDir{2}=[homeDir,'/data5/REMA'];
%regionDirs{2}=dir([upDir{2},'/region_*']);
%regionDirs{2}=cellfun(@(x) [upDir{2},'/',x,'/strips_unf/2m'], {regionDirs{2}.name},...
 %   'UniformOutput',false);

%regionDirs = cat(2,regionDirs{:});

upDir=['/mnt/pgc/data/elev/dem/setsm/ArcticDEM/region'];
regionDirs=dir([upDir,'/arcticdem_*']);
regionDirs=cellfun(@(x) [upDir,'/',x,'/strips_unf/2m'], {regionDirs([regionDirs.isdir]).name},...
    'UniformOutput',false);

meta=[];

if ~exist(outname,'file')
    fprintf('Creating new database: %s\n', outname)
else
    fprintf('Loading existing database to be appended to: %s\n', outname)
    out0=load(outname);
    [stripDirs0,~,~] = cellfun(@fileparts, out0.fileName, 'UniformOutput', false);
    [~,stripDnames0,~] = cellfun(@fileparts, stripDirs0, 'UniformOutput', false);
    stripDnames0 = unique(stripDnames0);
end

i=1;
for i=1:length(regionDirs)
        
    regionDir=regionDirs{i};
    
    if exist(regionDir,'dir')

%        if exist('out0','var')
%            [~,IA] = intersect(metaFiles, out0.fileName);
%            metaFiles(IA) = [];
%            if isempty(metaFiles)
%                continue
%            end
%        end
%
%        metaFiles=dir([regionDir,'/*_2m_lsf*/*meta.txt']);
%        metaFiles = strcat({metaFiles.folder}',repmat({'/'},length(metaFiles),1),{metaFiles.name}');
%
%        j=1;
%        for j=1:length(metaFiles)
%            metaFile=metaFiles{j};
%            fprintf('adding file %s\n',metaFile)
%            if isempty(meta)
%                meta=readStripMeta(metaFile,'noStripMeta');
%            else
%                meta(length(meta)+1)=readStripMeta(metaFile,'noStripMeta');
%            end
%        end

        stripDir_pattern=[regionDir,'/*_2m_lsf*'];
        fprintf('Gathering strips with pattern: %s ... ', stripDir_pattern)

        stripDirs=dir(stripDir_pattern);
        stripDirs = strcat({stripDirs.folder}',repmat({'/'},length(stripDirs),1),{stripDirs.name}');
        if length(stripDirs) == 0
            fprintf('None found\n')
            continue
        end
        [~,stripDnames,~] = cellfun(@fileparts, stripDirs, 'UniformOutput', false);

        test_dup_stripids=stripDnames;
        k=1;
        for k=1:length(test_dup_stripids)
            stripid = test_dup_stripids{k};
            stripid_parts = strsplit(test_dup_stripids{k}, '_');
            strip_verkey = stripid_parts(end);
            strip_verkey = strip_verkey{1};
            if length(strip_verkey) == 7 && strcmp(strip_verkey(1), 'v') && ~isnan(str2double(strip_verkey(2:7)))
                stripid = strrep(stripid, ['_',strip_verkey], '');
                test_dup_stripids{k} = stripid;
            end
        end
        [~, uniqueIdx] = unique(test_dup_stripids);
        test_dup_stripids(uniqueIdx) = [];
        if length(test_dup_stripids) > 0
            fprintf('\nERROR: Multiple strips exist matching strip ID:\n')
            test_dup_stripids
            fprintf('Exiting early -- Nothing was written to database\n')
            return
        end

        if exist('out0','var')
            [~,IA] = intersect(stripDnames, stripDnames0);
            stripDirs(IA) = [];
            if isempty(stripDirs)
                fprintf('No new strips to add\n')
                continue
            end
        end

        stripDirs_filtered={};
        k=0;
        j=1;
        for j=1:length(stripDirs)
            strippair_dir=stripDirs{j};

            finfilecheck=dir([strippair_dir,'/*.fin']);
            if isempty(finfilecheck); continue; end

            demfilecheck=dir([strippair_dir,'/*dem.tif']);
            if isempty(demfilecheck); continue; end

            k=k+1;
            stripDirs_filtered{k}=strippair_dir;
        end
        stripDirs=stripDirs_filtered;

        num_strips_to_add=length(stripDirs);
        fprintf('%d to add\n', num_strips_to_add)
%        continue
        
        k=1;
        last_print_len=0;
        for k=1:length(stripDirs)
            stripDir=stripDirs{k};

            fprintf(repmat('\b', 1, last_print_len));
            last_print_len=fprintf('Reading strip (%d/%d): %s',k,num_strips_to_add,stripDir);

            finFile=dir([stripDir,'/*.fin']);
            if isempty(finFile); continue; end

            metaFiles=dir([stripDir,'/*meta.txt']);
            if isempty(metaFiles); continue; end
            metaFiles = strcat({metaFiles.folder}',repmat({'/'},length(metaFiles),1),{metaFiles.name}');

            j=1;
            for j=1:length(metaFiles)
                metaFile=metaFiles{j};
%                fprintf('adding file %s\n',metaFile)
                if isempty(meta)
                    meta=readStripMeta(metaFile,'noStripMeta');
                else
                    meta(length(meta)+1)=readStripMeta(metaFile,'noStripMeta');
                end

            end
        end
        fprintf('\n')
    end
end

if isempty(meta)
    fprintf('\nNo new records to add to database\n')
    return
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

fprintf('Writing %d new records to database file\n',length(out.fileName));
if exist('out0','var')
    flds = fields(out);
    i=1;
    for i =1:length(flds)
        fld = flds{i};
        eval(['out.',fld,'= [out0.',fld,',out.',fld,'];']);
    end
end

fprintf('Writing %d total records to %s\n',length(out.fileName),outname_appended);
save(outname_appended,'-struct','out','-v7.3');
