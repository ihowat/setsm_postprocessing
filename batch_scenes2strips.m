regionnum='21';
res='2m';
demdir=dir(['/data2/ArcticDEM/region_',regionnum,'*']);
demdir=['/data2/ArcticDEM/',demdir(1).name,'/tif_results/',res];
outdir=strrep(demdir,'tif_results','strips');

%demdir=['/data1/ArcticDEM/region_06_greenland_northwest_grit_2m/tif_results/2m']
demdir=['/data/nga/Egypt/new/results/'];
outdir=['/data/nga/Egypt/new/results/strips'];
fprintf('%s\n',demdir)

% demdir=['/data/nga/Arctic/Greenland/SETSM_DEMs/ellyn/v2'];
%demdir=['/data2/NGA_tests_v2/Thule/tif_results/8m'];
%demdir=['/data2/ArcticDEM/region_32_alaska_aleutians_set2/tif_results/2m'];

if ~exist(outdir,'dir'); mkdir(outdir); end

f = dir([demdir,'/*_dem.tif']);

stripids=char(f.name);
s=find(stripids(1,:)=='_');
stripids = stripids(:,s(2)+1:s(4)-1);
stripids =cellstr(stripids);
unique_stripids =unique(stripids);

f={f.name};

strt=1;
inc=1;
for i=strt:inc:length(unique_stripids)
    
    %% find scenes with this strip pair id
    n = find(strcmp(unique_stripids{i},stripids));
    
    fprintf('merging pair id: %s, %d of %d, %d scenes\n',unique_stripids{i},i,length(unique_stripids),length(n))
    
    %existance check
    a=dir([outdir,'/',f{n(1)}(1:48),'seg*_',res,'_dem.tif']);
    
    if ~isempty(a)
        b=dir([demdir,'/',f{n(1)}(1:48),'*_datamask.tif']);
        
        a=min([a.datenum]); %date of strip creation
        b=max([b.datenum]); %date of data mask creation
        
        % if dem date is less than strip date continue
        if b < a; fprintf('younger strip exists, skipping\n'); continue; end
        
        fprintf('old strip exists, deleting and reprocessing\n');
        eval(['!rm -f ',outdir,'/',f{n(1)}(1:48),'*']);
        
    end
    
    %% make sure all ortho's and matchtags exist, if missing, skip
    missingflag=0;
    j=1;
    for j=1:length(n)
       if ~exist(strrep([demdir,'/',f{n(j)}],'dem.tif','matchtag.tif'),'file');
           fprintf('matchtag file for %s missing, skipping this strip\n',f{n(j)});
           missingflag=1;
       end
       if ~exist(strrep([demdir,'/',f{n(j)}],'dem.tif','ortho.tif'),'file');
           fprintf('ortho file for %s missing, skipping this strip\n',f{n(j)});
           missingflag=1;
       end
    end
    if missingflag==1 ;continue; end
    %%
    seg=1;
    while ~isempty(n)
    
        fprintf('building segment %d\n',seg)

        % send n scenes to strip mosaicker
        [x,y,z,m,o,trans,rmse,scene]=scenes2strips(demdir,f(n));

        % find which files were used
        [~,IA] = intersect(f(n),scene);
        
        % remove these from the list
        n(IA) = [];
         
        if isempty(z); continue; end;
        
        
        OutDemName=[outdir,'/',scene{1}(1:48),'seg',num2str(seg),'_',res];
    
        save([OutDemName,'_trans.mat'],'scene','trans','rmse')
         
%         if exist([OutDemName,'_dem.tif'],'file'); 
%             fprintf('%s exists, skipping \n',OutDemName);
%             continue; 
%         end
%         
%         if exist([OutDemName,'_trans.mat'],'file');
%             load([OutDemName,'_trans.mat'])
%             [x,y,z,m,o,trans,rmse]=scenes2strips_v2(demdir,fn{j},trans);
%         else
%             
%             [x,y,z,m,o,trans,rmse]=scenes2strips_v2(demdir,fn{j});
%         
%             scene=fn{j};
%             save([OutDemName,'_trans.mat'],'scene','trans','rmse')
%         end
        
        d=readGeotiff([demdir,'/',scene{1}],'MapInfoOnly');
        if d.Tinfo.GeoDoubleParamsTag(1) > 0;
            projstr='polar stereo north';
        else
            projstr='polar stereo south';
        end
        
        writeGeotiff([OutDemName,'_dem.tif'],x,y,z,4,-9999,projstr);
    
        OutMatchName=[OutDemName,'_matchtag.tif'];
        writeGeotiff(OutMatchName,x,y,m,1,0,projstr);
    
        OutOrthoName=[OutDemName,'_ortho.tif'];
        writeGeotiff(OutOrthoName,x,y,o,2,0,projstr);
        
        stripmeta([OutDemName,'_dem.tif'])
            
        seg = seg + 1;
    
    end
    
end
