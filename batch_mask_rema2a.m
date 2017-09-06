% batch_mask_rema2a: batch script for priducing edgemask and datamask files 
% from DEM scene pairs.

% Ian Howat, ihowa@gmail.com, Ohio State University

% DEM directory
demDir='/data3/REMA/region_04_ronne_shelf/tif_results/8m';

% start on file number
strt=1;

% and process every
inc=1; 


fprintf('processing %s\n',demDir)
fprintf('processing every %d files starting with %d\n',inc,strt);

% make filename list
demFiles = dir([demDir,'/*dem.tif']);
demFiles = {demFiles.name};
demFiles = cellfun( @(x) [demDir,'/',x], demFiles, 'uniformoutput', false);

% processing list
I = strt:inc:length(demFiles);

i=1;
for i=1:1:length(I)
    
    demFile = demFiles{I(i)};
    OutMaskName = strrep(demFile,'dem.tif','mask2a.tif');
    
    fprintf('processing %d of %d: %s \n',i,length(I),demFile)
    
    if exist(OutMaskName,'file');  fprintf('mask exists, skipping\n'); continue; end
    
    m = remaMask2a(demFile);
        
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