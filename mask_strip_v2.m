function mask_strip_v2(demDir,stripid)

demFiles = dir([demDir,'/*',stripid,'*_dem.tif']);
demFiles = {demFiles.name};
demFiles = cellfun( @(x) [demDir,'/',x], demFiles, 'uniformoutput', false);

i=1;
for i=1:1:length(demFiles)
    
    demFile = demFiles{i};
    OutMaskName = strrep(demFile,'dem.tif','mask.tif');
    fprintf('processing %d of %d: %s \n',i,length(demFiles),demFile)
    if exist(OutMaskName,'file');  fprintf('mask exists, skipping\n'); continue; end
   
    %read meta file
    metaFile= strrep(demFile,'dem.tif','meta.txt');
    meta=readSceneMeta(metaFile);
    
    % get satid and check for new cross track naming convention
    [~,satID]=fileparts(meta.image_1);
    satID = upper(satID(1:4));
    
    if strcmp(satID(1:2),'W1'); satID = 'WV01'; end 
    if strcmp(satID(1:2),'W2'); satID = 'WV02'; end
    if strcmp(satID(1:2),'W3'); satID = 'WV03'; end
    if strcmp(satID(1:2),'G1'); satID = 'GE01'; end 
    if strcmp(satID(1:2),'Q1'); satID = 'QB01'; end
    if strcmp(satID(1:2),'Q2'); satID = 'QB02'; end
    if strcmp(satID(1:2),'I1'); satID = 'IK01'; end
    
    meta.image_1_satID = satID;
    
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

