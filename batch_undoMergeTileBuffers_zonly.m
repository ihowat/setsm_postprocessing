function undoMergeTileBuffers_zonly(fileNames)
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
    check_file = fileNames{i};
    if ~isfile(check_file)
        error('Source file does not exist: %s\n', check_file)
    end

    m = matfile(check_file);

    varlist = who(m);
    if ~any(strcmp(varlist,'zbuff'))
        fprintf('Tile has not been merged, skipping %s\n', check_file)
        continue
    end

    zbuff = m.zbuff;

    if all(cellfun(@(x) isempty(x), zbuff))
        fprintf('Tile zbuff is empty, skipping %s\n', check_file)
        continue
    end

    fprintf('De-buffering %s\n',check_file)

    m.Properties.Writable = true;

    if any(strcmp(varlist,'mergedTop')) && m.mergedTop
        [nrows,~]=size(zbuff{1});
        m.z(1:nrows,:) = zbuff{1};
        m.mergedTop = false;
    end

    if any(strcmp(varlist,'mergedBottom')) && m.mergedBottom
        [nrows,~]=size(zbuff{2});
        m.z(end-(nrows-1):end,:) = zbuff{2};
        m.mergedBottom=false;
    end

    if any(strcmp(varlist,'mergedLeft')) && m.mergedLeft
        [~,ncols]=size(zbuff{3});
        m.z(:,1:ncols) = zbuff{3};
        m.mergedLeft=false;
    end

    if any(strcmp(varlist,'mergedRight')) && m.mergedRight
        [~,ncols]=size(zbuff{4});
        m.z(:,end-(ncols-1):end) = zbuff{4};
        m.mergedRight=false;
    end

    m.zbuff = cell(4,1);

end



