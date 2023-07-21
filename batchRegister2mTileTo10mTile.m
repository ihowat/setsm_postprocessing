function batchRegister2mTileTo10mTile(tileDir,varargin)
% set paths

%     %   tileDir='/Users/ihowat/project/howat.4/rema_mosaic/rema_mosaic_v13e';
%     is2DirOrFile='/Users/ihowat/data4/REMA/altimetryByTile';

%     %    tileDir='/fs/project/howat.4/howat.4/rema_mosaic_v13e/rema_19_victoria_land';
%     is2DirOrFile='/fs/byo/howat-data4/REMA/altimetryByTile';

strt=1;
inc=1;
resolution='2m';
overwrite=false;
overwrite_flag='';

if length(varargin) == 2 && ~isnan(str2double(varargin{1})) && ~isnan(str2double(varargin{2}))
    strt=varargin{1};
    inc= varargin{2};

elseif ~isempty(varargin)
    n=find(strcmpi(varargin,'start'));
    if ~isempty(n)
        strt=varargin{n+1};
    end

    n=find(strcmpi(varargin,'inc'));
    if ~isempty(n)
        inc=varargin{n+1};
    end

    if any(strcmpi(varargin,'overwrite'))
        overwrite=true;
        overwrite_flag='overwrite';
    end
end

fprintf('Processing every %d tile, starting at %d\n',inc,strt)

tileFiles = dir(sprintf('%s/*_%s.mat', tileDir, resolution));
%tileFiles = dir([tileDir,'/*.mat']);
tileNames = {tileFiles.name};
tileFiles = fullfile({tileFiles.folder}, tileNames);
tileNames10m = cellfun(@(x) get10mTileName(x), tileNames, 'UniformOutput',false);

tileFiles=tileFiles(:);
tileNames10m=tileNames10m(:);

i=strt;
for i=strt:inc:length(tileFiles)

    % get this tile name and set aux tile file names from it
    tileFile=tileFiles{i};
    if endsWith(tileFile, '_reg.mat')
        continue
    end
    tileName10m = tileNames10m{i};

    tileFile10m = [tileDir,'/',tileName10m,'_10m.mat'];
    tileFile10m_reg = [tileDir,'/',tileName10m,'_10m_reg.mat'];
    if exist(tileFile10m_reg,'file') == 2
        tileFile10m = tileFile10m_reg;
    end

    fprintf('Working on tile %d of %d: %s\n',i,length(tileFiles),tileFile);

    unregTileFile=strrep(tileFile,'.mat','_unreg.mat.bak');
    regTempFile=strrep(tileFile,'.mat','_reg.mat.temp');
    regTileFile=strrep(tileFile,'.mat','_reg.mat');
    regFailFile=strrep(tileFile,'.mat','_reg.mat.regfail');

    if exist(unregTileFile,'file') ~= 2
        fprintf('Creating backup copy of unaltered tile matfile: %s\n',unregTileFile);
        eval(['!cp --preserve=timestamps ',tileFile,' ',unregTileFile]);
    end

    if exist(regTileFile,'file') == 2
        if overwrite
            fprintf('reg matfile already exists, will overwrite: %s\n',regTileFile);
        else
            fprintf('reg matfile already exists, skipping: %s\n',regTileFile);
            continue;
        end
    end

    missing_aux_file = false;
    if exist(tileFile10m,'file') ~= 2
        fprintf('10m DEM file does not exist: %s\n',tileFile10m);
        missing_aux_file = true;
    end
    if missing_aux_file
        fprintf('skipping registration: %s\n',unregTileFile);
        continue;
    end

    if exist(regTileFile,'file') ~= 2 || overwrite
        try
            if exist(regTileFile, 'file') == 2
                fprintf('Removing existing reg.mat file: %s\n', regTileFile);
                delete(regTileFile);
            end

            fprintf('Registering 2m to 10m: %s, %s\n', tileFile, tileFile10m);
            regFail = register2mTileTo10mTile(tileFile, tileFile10m);

            if regFail
                if exist(regTileFile, 'file') == 2
                    fprintf('Removing failure reg.mat file: %s\n', regTileFile);
                    delete(regTileFile);
                end
                fprintf('Tile cannot be registered, writing semaphore file: %s\n', regFailFile);
                fid = fopen(regFailFile, 'w');
                fclose(fid);

            elseif exist(regTileFile,'file') ~= 2
                error('Expected output reg.mat file does not exist: %s\n', regTileFile);
            end
        catch ME
            if exist(regTileFile, 'file') == 2
                fprintf('Detected error, removing reg.mat file: %s\n',regTileFile);
                delete(regTileFile);
            end
            rethrow(ME);
        end
    end

end
