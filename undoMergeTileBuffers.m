function undoMergeTileBuffers(fileNames)
% undoMergeTileBuffers reverts the tile edge buffers to the original values
%
% undoMergeTileBuffers(fileNames) replaces the merged tile buffers, created
% by mergeTileBuffer, with the original values, setting buffer flags to
% false. 

% make cell if not
if ~iscell(fileNames)
    fileNames={fileNames};
end

%%
i=1;
for i =1:length(fileNames)
    fprintf('De-buffering %s\n',fileNames{i})
    m=matfile(fileNames{i});
    
    m.Properties.Writable = true;
    
    zbuff = m.zbuff;
    
    if m.mergedTop
        m.z(1:21,:) = zbuff{1};
        m.mergedTop = false;
    end
    
    if m.mergedBottom
        m.z(end-20:end,:) = zbuff{2};
        m.mergedBottom=false;
    end
    
    if m.mergedLeft
        m.z(:,1:21) = zbuff{3};
        m.mergedLeft=false;
    end
    
    if m.mergedRight
        m.z(:,end-20:end) = zbuff{4};
        m.mergedRight=false;
    end
    
    m.zbuff = cell(4,1);
    
end

        
        
        