% batch_mask: batch script for priducing edgemask and datamask files 
% from DEM scene pairs.

% Ian Howat, ihowa@gmail.com, Ohio State University

%% ArcticDEM input directory settings
% updir='/data3/REMA'; %location of region directory
% regionnum='19'; % ArcticDEM region #
% res='8'; % DEM resolution

updir='/data2/ArcticDEM'; %location of region directory
regionnum='19'; % ArcticDEM region #
res='2'; % DEM resolution


%% load file names
demDir=dir([updir,'/region_',regionnum,'*']);
demDir=[updir,'/',demDir(1).name,'/tif_results/',res,'m'];

fprintf('working: %s\n',demDir);

demFiles = dir([demDir,'/*_dem.tif']);
demDates = [demFiles.datenum];
demFiles = {demFiles.name};
demFiles = cellfun( @(x) [demDir,'/',x], demFiles, 'uniformoutput', false);

%% Update Mode - will only reprocess masks older than the matchtag file
maskFiles = dir([demdir,'/*_mask.tif']);

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
