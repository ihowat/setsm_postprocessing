function mask_strip(demdir,stripid)

f = dir([demdir,'/*',stripid,'*_matchtag.tif']);
fdate=[f.datenum];
f={f.name};

%% Update Mode - will only reprocess masks older than the matchtag file
fedge = dir([demdir,'/*_edgemask.tif']);

if ~isempty(fedge)
    fedgeDate=[fedge.datenum];
    fedge={fedge.name};
    [~,IA,IB] = intersect(f,strrep(fedge,'edgemask.tif','matchtag.tif'));
    n= fedgeDate(IB) - fdate(IA) >= -6.9444e-04;
    f(IA(n))=[];
    
    clear fdate fedge fedgeDate
end

i=1;
for i=1:1:length(f)
    
    matchFile = [demdir,'/',f{i}];
    orthoFile = strrep(matchFile,'matchtag.tif','ortho.tif');
    
    fprintf('processing %d of %d: %s \n',i,length(f),matchFile)
    OutEdgeMaskName= strrep(matchFile,'matchtag.tif','edgemask.tif');
    OutDataMaskName= strrep(matchFile,'matchtag.tif','datamask.tif');
    
    m=readGeotiff(matchFile);
    
    fprintf('making %s \n',OutEdgeMaskName)
    
    %find SETSM version
    metaFileName=strrep(matchFile,'matchtag.tif','meta.txt');
    if ~exist(metaFileName,'file')
        disp('no meta file, assumimg SETSM version > 2.0');
        setsmVersion=3;
    else
        c=textread(metaFileName,'%s','delimiter','\n');
        r=find(~cellfun(@isempty,strfind(c,'SETSM Version='))); 
        if ~isempty(r)
		setsmVersion=deblank(strrep(c{r(1)},'SETSM Version=',''));
        else; setsmVersion='2.03082016'; end
	fprintf('Using settings for SETSM Version = %s\n',setsmVersion)
        setsmVersion=str2num(setsmVersion);
    end
    
    if setsmVersion < 2.01292016
        n= floor(21*2/m.info.map_info.dx);
        Pmin=0.8;
        Amin=2000/m.info.map_info.dx;
        crop=n;
    else
        n= floor(101*2/m.info.map_info.dx);
        Pmin=0.99;
        Amin=2000/m.info.map_info.dx;
        crop=n;
    end
    
    % get water mask
    M0 = entropyMask(orthoFile);
 
    m1.z=edgeMask(m.z,'n',n,'Pmin',Pmin,'Amin',Amin,'crop',crop,'m0',M0);
    
    if isfield(m.Tinfo,'GeoDoubleParamsTag')
        
        if m.Tinfo.GeoDoubleParamsTag(1) > 0
            projstr='polar stereo north';
        else
            projstr='polar stereo south';
        end
        
    else
        
        projstr=m.Tinfo.GeoAsciiParamsTag;
        a=findstr( projstr, 'Zone');  
        b=findstr( projstr, ',');
        c = findstr( projstr,'Northern Hemisphere');
      
        if ~isempty(c)
            projstr=[projstr(a+4:b-1),' North']; 
        else
            projstr=[projstr(a+4:b-1),' South'];
        end
    end
    
    writeGeotiff(OutEdgeMaskName,m.x,m.y,m1.z,1,0,projstr)
    
    m.z(~m1.z)=uint8(0);
    
    clear m1

    fprintf('making %s \n',OutDataMaskName)
    
    if setsmVersion <= 2.0
        
        n= floor(21*2/m.info.map_info.dx);
        Pmin=0.3;
        Amin=1000;
        Amax=10000;
        
    else
        n= floor(101*2/m.info.map_info.dx);
        Pmin=0.90; 
        Amin=1000;
        Amax=1000;
        
    end
    
    m.z=DataDensityMask(m.z,'n',n,'Pmin',Pmin,'Amin',Amin,'Amax',Amax);
   
    writeGeotiff(OutDataMaskName,m.x,m.y,m.z,1,0,projstr)
    
end


