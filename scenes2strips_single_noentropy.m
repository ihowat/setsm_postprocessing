function scenes2strips_single_noentropy(demdir,stripid,res,outdir)

% change this
res=[res,'m'];
fprintf('source: %s\n',demdir)
fprintf('res: %s\n',res)

if ~exist(outdir,'dir'); mkdir(outdir); end

f = dir([demdir,'/*_dem.tif']);

stripids=char(f.name);
s=find(stripids(1,:)=='_');
stripids = stripids(:,s(2)+1:s(4)-1);
stripids =cellstr(stripids);

f={f.name};

%% find scenes with this strip pair id
n = find(strcmp(stripid,stripids));

fprintf('merging pair id: %s, %d scenes\n',stripid,length(n))

%existance check
a=dir([outdir,'/',f{n(1)}(1:48),'seg*_dem.tif']);
if ~isempty(a); return; end

% make sure all ortho's and matchtags exist, if missing, skip
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
if missingflag==1; return; end

%% run filtering
mask_strip_noentropy(demdir,stripid);

%% build strip segments
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
    fprintf('DEM: %s\n',OutDemName)
    save([OutDemName,'_trans.mat'],'scene','trans','rmse')

    d=readGeotiff([demdir,'/',scene{1}],'MapInfoOnly');
    if isfield(d.Tinfo,'GeoDoubleParamsTag')
            
        if d.Tinfo.GeoDoubleParamsTag(1) > 0
            projstr='polar stereo north';
        else
            projstr='polar stereo south';
        end

    else

        projstr=d.Tinfo.GeoAsciiParamsTag;
        a=findstr( projstr, 'Zone');
        b=findstr( projstr, ',');
        c = findstr( projstr,'Northern Hemisphere');

        if ~isempty(c)
            projstr=[projstr(a+4:b-1),' North'];
        else
            projstr=[projstr(a+4:b-1),' South'];
        end
    end

    writeGeotiff([OutDemName,'_dem.tif'],x,y,z,4,-9999,projstr);

    OutMatchName=[OutDemName,'_matchtag.tif'];
    writeGeotiff(OutMatchName,x,y,m,1,0,projstr);

    OutOrthoName=[OutDemName,'_ortho.tif'];
    writeGeotiff(OutOrthoName,x,y,o,2,0,projstr);

    stripmeta([OutDemName,'_dem.tif'])

    seg = seg + 1;
    
    tempfile = [OutDemName,'_dem_temp.tif'];
    hillshade = [OutDemName,'_dem_10m_shade_masked.tif'];
    system(['gdal_translate -q -tr 10 10 -r bilinear -co bigtiff=if_safer -a_nodata -9999 ',...
        [OutDemName,'_dem.tif'],' ', tempfile]);
    
    system(['gdaldem hillshade -q -z 3 -compute_edges -of GTiff -co TILED=YES -co BIGTIFF=IF_SAFER -co COMPRESS=LZW ',...
        tempfile,' ', hillshade]);
    
    delete(tempfile);

end
