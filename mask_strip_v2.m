function mask_strip_v2(demDir,stripid)

demFiles = dir([demDir,'/*',stripid,'*_dem.tif']);
demDates = [demFiles.datenum];
demFiles = {demFiles.name};
demFiles = cellfun( @(x) [demDir,'/',x], demFiles, 'uniformoutput', false);

%% Update Mode - will only reprocess masks older than the matchtag file
maskFiles = dir([demDir,'/*',stripid,'*mask.tif']);

if ~isempty(maskFiles)
    
    maskDates=[maskFiles.datenum];
    maskFiles={maskFiles.name};
    maskFiles = cellfun( @(x) [demDir,'/',x], maskFiles, 'uniformoutput', false);
    [~,IA,IB] = intersect(demFiles,strrep(maskFiles,'mask.tif','dem.tif'));
    n= maskDates(IB) - demDates(IA) >= -6.9444e-04;
    demFiles(IA(n))=[];
    
    clear demDates maskFiles maskDates
end

i=1;
for i=1:1:length(demFiles)
    
    demFile = demFiles{i};
    OutMaskName = strrep(demFile,'dem.tif','mask.tif');
    fprintf('processing %d of %d: %s \n',i,length(demFiles),demFile)
   
    %read meta file
    metaFile= strrep(demFile,'dem.tif','meta.txt');
    meta=readSceneMeta(metaFile);
    m = mask(demFile,meta);
    
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
    
    writeGeotiff(OutMaskName,m.x,m.y,m.z,1,0,projstr)
    
end

