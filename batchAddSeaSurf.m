function batchAddSeaSurf(tileDir,projstr,varargin)
% set paths

%     %   tileDir='/Users/ihowat/project/howat.4/rema_mosaic/rema_mosaic_v13e';
%     is2DirOrFile='/Users/ihowat/data4/REMA/altimetryByTile';

%     %    tileDir='/fs/project/howat.4/howat.4/rema_mosaic_v13e/rema_19_victoria_land';
%     is2DirOrFile='/fs/byo/howat-data4/REMA/altimetryByTile';

strt=1;
inc=1;
resolution='10m';
overwrite=false;
overwrite_flag='';

if strcmpi(projstr, 'polar stereo north')
    addSeaSurface_epsg = 3413;
elseif strcmpi(projstr, 'polar stereo south')
    addSeaSurface_epsg = 3031;
elseif contains(projstr,'UTM','IgnoreCase',true)
    [addSeaSurface_epsg,~,~] = getProjstrInfo(projstr);
else
    error("Cannot determine EPSG code for argument projstr: '%s'", projstr);
end

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

    resolution_choices = {'10m', '2m'};
    n = find(strcmpi('resolution', varargin));
    if ~isempty(n)
        resolution = varargin{n+1};
    end
    if ~ismember(resolution, resolution_choices)
        error("'resolution' must be one of the following, but was '%s': {'%s'}", resolution, strjoin(resolution_choices, "', '"));
    end

    if any(strcmpi(varargin,'overwrite'))
        overwrite=true;
        overwrite_flag='overwrite';
    end
end

fprintf('Processing every %d tile, starting at %d\n',inc,strt)

regFiles = dir(sprintf('%s/*_%s_reg.mat', tileDir, resolution));
unregFiles = dir(sprintf('%s/*_%s.mat', tileDir, resolution));
regFiles = fullfile({regFiles.folder}, {regFiles.name});
unregFiles = fullfile({unregFiles.folder}, {unregFiles.name});
unregFiles = setdiff(unregFiles, strrep(regFiles, '_reg.mat', '.mat'));
tileFiles = [unregFiles regFiles];

tileNames10m = cellfun(@(x) get10mTileName(x), tileFiles, 'UniformOutput',false);

tileFiles=tileFiles(:);
tileNames10m=tileNames10m(:);

i=strt;
for i=strt:inc:length(tileFiles)

    % get this tile name and set aux tile file names from it
    tileFile=tileFiles{i};
    tileName10m = tileNames10m{i};
    
    fprintf('Working on tile %d of %d: %s\n',i,length(tileFiles),tileFile);

    unfillTileFile=strrep(tileFile,'.mat','_unfill.mat.bak');
    fillTempFile=strrep(tileFile,'.mat','_fill.mat.temp');
    fillTileFile=strrep(tileFile,'.mat','_fill.mat');

    if exist(unfillTileFile,'file') ~= 2
        fprintf('Creating backup copy of unaltered tile matfile: %s\n',unfillTileFile);
        eval(['!cp --preserve=timestamps ',tileFile,' ',unfillTileFile]);
    end
    
    if exist(fillTileFile,'file') == 2
        if overwrite
            fprintf('fill matfile already exists, will overwrite: %s\n',fillTileFile);
        else
            fprintf('fill matfile already exists, skipping: %s\n',fillTileFile);
            continue;
        end
    end

    if exist(fillTileFile,'file') ~= 2 || overwrite
        try
            if exist(fillTileFile, 'file') == 2
                fprintf('Removing existing fill.mat file: %s\n', fillTileFile);
                delete(fillTileFile);
            end
            if exist(fillTempFile, 'file') == 2
                fprintf('Removing existing fill.mat.tmp file: %s\n', fillTempFile);
                delete(fillTempFile);
            end

            fprintf('Calculating waterfill\n');
            m = matfile(tileFile);
            m_varlist = who(m);
            if ismember('reg_extrap_nodata', m_varlist)
                [z_fill,waterfill_mask] = addSeaSurfaceHeight(m.x,m.y,m.z,m.land,'epsg',addSeaSurface_epsg,'adaptCoastline','doNotFillBorder',m.reg_extrap_nodata);
            else
                [z_fill,waterfill_mask] = addSeaSurfaceHeight(m.x,m.y,m.z,m.land,'epsg',addSeaSurface_epsg,'adaptCoastline');
            end

            fprintf('Copying unfill.mat file to new fill.mat.tmp file: %s\n',fillTempFile);
            eval(['!cp ',unfillTileFile,' ',fillTempFile]);

            fprintf('Writing filled z to fill.mat.tmp file\n');
            m1=matfile(fillTempFile);
            m1.Properties.Writable = true;
            m1.z = z_fill;
            m1.waterFillMask = waterfill_mask;
            m1.waterFilled = true;
            m1.aux_waterFill_waterMaskTif = "'land' variable";

            fprintf('Renaming fill.mat.tmp file to fill.mat: %s\n',fillTileFile);
            eval(['!mv ',fillTempFile,' ',fillTileFile]);
        catch ME
            if exist(fillTileFile, 'file') == 2
                fprintf('Detected error, removing fill.mat file: %s\n',fillTileFile);
                delete(fillTileFile);
            end
            if exist(fillTempFile, 'file') == 2
                fprintf('Detected error, removing fill.mat.tmp file: %s\n',fillTempFile);
                delete(fillTempFile);
            end
            rethrow(ME);
        end
    end

end
