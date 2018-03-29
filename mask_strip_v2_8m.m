function mask_strip_v2_8m(demDir,stripid)

demFiles = dir([demDir,'/*',stripid,'*_dem.tif']);
demFiles = {demFiles.name};
demFiles = cellfun( @(x) [demDir,'/',x], demFiles, 'uniformoutput', false);

i=1;
for i=1:1:length(demFiles)
    
    demFile = demFiles{i};
    OutMaskName = strrep(demFile,'dem.tif','mask.tif');
    fprintf('processing %d of %d: %s \n',i,length(demFiles),demFile)
    if exist(OutMaskName,'file');  fprintf('mask exists, skipping\n'); continue; end
    
    m = mask8m(demFile);
    
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

