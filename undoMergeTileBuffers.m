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
    m=matfile(fileNames{i});

    varlist = who(m);
    if ~any(strcmp(varlist,'zbuff'))
        fprintf('Tile has not been merged, skipping %s\n',fileNames{i})
        continue
    end

    fprintf('De-buffering %s\n',fileNames{i})
    
    m.Properties.Writable = true;
    
    zbuff = m.zbuff;
    
    if m.mergedTop
        [nrows,~]=size(zbuff{1});
        m.z(1:nrows,:) = zbuff{1};
        m.mergedTop = false;
    end
    
    if m.mergedBottom
        [nrows,~]=size(zbuff{2});
        m.z(end-(nrows-1):end,:) = zbuff{2};
        m.mergedBottom=false;
    end
    
    if m.mergedLeft
        [~,ncols]=size(zbuff{3});
        m.z(:,1:ncols) = zbuff{3};
        m.mergedLeft=false;
    end
    
    if m.mergedRight
        [~,ncols]=size(zbuff{4});
        m.z(:,end-(ncols-1):end) = zbuff{4};
        m.mergedRight=false;
    end
    
    m.zbuff = cell(4,1);
    
end

        
        
