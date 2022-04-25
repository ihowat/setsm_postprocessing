function cropTile(tileFile_or_list)

%tileFile='/Users/ihowat/data4/REMA/test_2m_tiles/tiles/37_10_1_1_2m.mat';

if iscell(tileFile_or_list)
    tileFile_list=tileFile_or_list;
else
    tileFile_list={tileFile_or_list};
end

for tileFile_idx=1:length(tileFile_list)
    tileFile = tileFile_list{tileFile_idx};
    if endsWith(tileFile, '_cropped.mat')
        % we don't want to attempt to crop these!
        continue;
    end
    tileFile_cropped = strrep(tileFile, '.mat', '_cropped.mat');

    fprintf("\nWorking on tile file: %s\n", tileFile);
    if isfile(tileFile_cropped)
        fprintf("Cropped tile already exists, skipping: %s\n", tileFile_cropped);
        continue;
    end

    m = matfile(tileFile);
    res = m.x(1,2) - m.x(1,1);

    % number of buffer pixels each side is supposed to have - this is the
    % number of pixels on the outside of edge of the quad tile/subtile
    if res == 10
        NPixBuff = 10;
    elseif res == 2
        NPixBuff = 50;
    else
        fprintf("ERROR: Tile resolution (%gm) is not supported, skipping this tile\n", res);
        continue;
    end

    %% Undo buffer merge if necessary.
    %% May want to do this only if tile needs to be cropped.
    fprintf("Undoing tile buffer merge\n")
    undoMergeTileBuffersSingle(tileFile, 'all', 'ignoreMergedVars',true)
%    undoMergeTileBuffersSingle(tileFile, 'all', 'ignoreMergedVars',false)

    % index of coordinate at the edge of the tile/quadtile boundary
    tileBoundaryIndX = find(mod(m.x,50e3) == 0);
    tileBoundaryIndY = find(mod(m.y,50e3) == 0);
    boundary_detection_issue = false;
    if res == 10
        if length(tileBoundaryIndX) ~= 3 && length(tileBoundaryIndY) ~= 3
            fprintf(strjoin(...
                ["ERROR: Expected mod(m.x,50e3) and mod(m.y,50e3) for 10m tile to produce three coordinate results,",...
                 "but found %g for x and %g for y\n"]...
            ), length(tileBoundaryIndX), length(tileBoundaryIndY));
            boundary_detection_issue = true;
        end
    elseif res == 2
        if length(tileBoundaryIndX) ~= 2 && length(tileBoundaryIndY) ~= 2
            fprintf(strjoin(...
                ["ERROR: Expected mod(m.x,50e3) and mod(m.y,50e3) for 2m tile to produce two coordinate results,",...
                 "but found %g for x and %g for y\n"]...
            ), length(tileBoundaryIndX), length(tileBoundaryIndY));
            boundary_detection_issue = true;
        end
    end
    if boundary_detection_issue
        fprintf("x range = [%g, %g]\n", m.x(1,1), m.x(1,end));
        fprintf("y range = [%g, %g]\n", m.y(1,end), m.y(1,1));
        fprintf("Skipping this tile\n");
        continue;
    end
    tileBoundaryIndX = [tileBoundaryIndX(1), tileBoundaryIndX(end)];
    tileBoundaryIndY = [tileBoundaryIndY(1), tileBoundaryIndY(end)];

    %number of buffer pixels on each side
    NPixBuffLeft  = tileBoundaryIndX(1) - 1;
    NPixBuffRight = length(m.x) - tileBoundaryIndX(2);
    NPixBuffTop  = tileBoundaryIndY(1) - 1;
    NPixBuffBot = length(m.y) - tileBoundaryIndY(2);

    % check if any buffers are larger than they should be
    if ~any([NPixBuffLeft NPixBuffRight NPixBuffTop NPixBuffBot] > NPixBuff)
        fprintf("Tile does not need to be cropped\n");
    else
        fprintf("Tile needs to be cropped\n");

%        %% Undo buffer merge if necessary.
%        %% May want to do this even if tile doesn't need to be cropped.
%        fprintf('Undoing tile buffer merge: %s\n', tileFile)
%        undoMergeTileBuffersSingle(tileFile, 'all', 'ignoreMergedVars',true)
%%        undoMergeTileBuffersSingle(tileFile, 'all', 'ignoreMergedVars',false)

        m_varlist = who(m);

        mergedTop = ismember('mergedTop', m_varlist) && m.mergedTop;
        mergedBottom = ismember('mergedBottom', m_varlist) && m.mergedBottom;
        mergedLeft = ismember('mergedLeft', m_varlist) && m.mergedLeft;
        mergedRight = ismember('mergedRight', m_varlist) && m.mergedRight;

        if mergedTop || mergedBottom || mergedLeft || mergedRight
            error("Tile 'merged{Right|Left|Top|Bottom}=true', indicating failure of undo tile buffer merge");
        end

        % get index ranges of correct buffers
        IndX = tileBoundaryIndX + [-NPixBuff NPixBuff];
        IndY = tileBoundaryIndY + [-NPixBuff NPixBuff];

        % get list of variables in file
        flds=fields(m);

        % remove Properties (skip) and x and y variables (handled differently)
        flds(strcmp(flds,'Properties')) = [];
        flds(strcmp(flds,'x')) = [];
        flds(strcmp(flds,'y')) = [];

        % remove buffer variables which will be invalid after crop
        flds(endsWith(flds,'buff')) = [];

        % create new matfile
        fprintf("Writing cropped tile file: %s\n", tileFile_cropped);
        m1 = matfile(tileFile_cropped);

        % write cropped x and y vectors
        m1.x = m.x(1,IndX(1):IndX(2));
        m1.y = m.y(1,IndY(1):IndY(2));

        % get size of z array for checking variables
        sz = size(m,'z');

        % loop through variables
        i=1;
        for i=1:length(flds)

            % check size
            if all(size(m,flds{i}) == sz)
                % if same dims as z, write cropped version to new matfile
                eval(['m1.',flds{i},' = m.',flds{i},...
                    '(IndY(1):IndY(2),IndX(1):IndX(2));'])
            else
                % if not same size, just xfer directly
                eval(['m1.',flds{i},' = m.',flds{i},';'])
            end
        end

    end
end
