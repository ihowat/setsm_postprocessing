
tileFile='/Users/ihowat/data4/REMA/test_2m_tiles/tiles/37_10_1_1_2m.mat';

m = matfile(tileFile);

% number of buffer pixels each side is supposed to have - this is the
% number of pixels on the outside of edge of the quad tile/subtile
NPixBuff = 50;

% index of coordinate at the edge of the tile/quadtile boundary
tileBoundaryIndX = find(mod(m.x,50e3) == 0);
tileBoundaryIndY = find(mod(m.y,50e3) == 0);

%number of buffer pixels on each side
NPixBuffLeft  = tileBoundaryIndX(1) - 1;
NPixBuffRight = length(m.x) - tileBoundaryIndX(2);
NPixBuffTop  = tileBoundaryIndY(1) - 1;
NPixBuffBot = length(m.y) - tileBoundaryIndY(2);

% check if any buffers are larger than they should be
if any([NPixBuffLeft NPixBuffRight NPixBuffTop NPixBuffBot] > NPixBuff)

    % get index ranges of correct buffers
    IndX = tileBoundaryIndX + [-NPixBuff NPixBuff];
    IndY = tileBoundaryIndY + [-NPixBuff NPixBuff];

    % get list of variables in file
    flds=fields(m);

    % remove Properties (skip) and x and y variables (handled differently)
    flds(contains(flds,'Properties')) = [];
    flds(contains(flds,'x')) = [];
    flds(contains(flds,'y')) = [];

    %  create new matfile
    m1 = matfile(strrep(tileFile,'.mat','_corrected.mat'));

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
