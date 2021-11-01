function batch_applyRegistration(varargin)
%batch_applyRegistration initializes the calls to applyRegistration, which
%shifts and interpolates the tile mosaics (all rasters) within mat files.
% Currently set to operate on tiles within a region directory.
% Makes a new .mat file called *_reg.mat
%
% Compatible with v4 mosaics.
%
% Should be resolution agnostic.

%% Paths/files
% upper level directory where region directories sit
%remaTileDir='~/Desktop';
%remaTileDir='/fs/project/howat.4/howat.4/automosaic13e';
remaTileDir='~/project/howat.4/rema_mosaic_v13e';

% region name (which is the name of the region directory)
regionName='rema_21_mbl_north';
resstr='10m';

% File with list of registration parameters, indexed by tile name
%registrationFile=  /fs/byo/howat-data4/REMA/altimetryByTile/rema_{region number}_{region name}_is2reg.mat
registrationFile= '~/data4/REMA/altimetryByTile/rema_21_mbl_north_is2reg.mat';

%% Read file list
demMatFiles=dir([remaTileDir,'/',regionName,'/*_',resstr,'.mat']);
demMatFiles=cellfun( @(x) [remaTileDir,'/',regionName,'/',x],...
    {demMatFiles.name}, 'uniformoutput',0);

if nargin == 3
    % If arguments are supplied:
    %  Arg 1: tile dir (/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/results/output_tiles_symlink_by_region)
    %  Arg 2: region name (rema_21_mbl_north)
    %  Arg 3: resolution string (2m)
    remaTileDir=varargin{1};
    regionName=varargin{2};
    registrationFile=['/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/from_unity/altimetryByTile/',...
        regionName, '_is2reg.mat'];
    resstr=varargin{3};

    %% Read file list
    demMatFiles=dir([remaTileDir,'/',regionName,'/*/*_',resstr,'.mat']);
    demMatFiles=cellfun( @(x) [remaTileDir,'/',regionName,'/',x(1:5),'/',x],...
        {demMatFiles.name}, 'uniformoutput',0);
end

%% loop through tiles
i=1;
for i=1:length(demMatFiles)
    
    % File to register
    demMatFile=demMatFiles{i};
    
    if ~exist(demMatFile,'file')
        warning('%s doesnt exist, skipping',demMatFile)
        continue
    end
    
    % get the region and tile names from the mat file name
    S=strsplit(demMatFile,'/');
    %regionName=S{end-1}; already set this above
    tileName=S{end}(1:5);

    fprintf('tile %s\n', S{end});
    
    %find this tile
    regData=load(registrationFile);
    n = find(contains(regData.tileNames,tileName));

    if isempty(n)
        warning('could not find registration data for %s, skipping',tileName)
        continue
    end

    % extract data for this tile from registration data
    regData = rmfield(regData,{'tileFiles','tileNames'});
    regData.reg = structfun( @(x) x(n,:), regData.reg,'uniformoutput',0);
    regData.unreg = structfun( @(x) x(n,:), regData.unreg,'uniformoutput',0);
    
    if isnan(regData.reg.p(1))
        warning('NaN offset for %s, skipping',tileName)
        continue
    end
    
    % load demMatFile into a matfile object
    m = matfile(demMatFile);
    
    % get grid size
    [nrows,ncols] = size(m,'z');
    
    % make coregistration cluster mask - not currently used
    C = true(nrows,ncols);
    
    % make out mat name to write to
    outName=strrep(demMatFile,'.mat','_reg.mat');

    if exist(outName)
        fprintf("reg matfile already exists\n")

    else
        % if only vertical registration only z is changed.
        if regData.reg.p(2) == 0 && regData.reg.p(3) == 0

            % copy the matfile to the outname
    %        eval(['cp ',demMatFile,' ',outName]);
            copyfile(demMatFile,outName);

            % subtract the vertical offset from the dem
            z = m.z - regData.reg.p(1);

            % append the z (will overwrite) and the regdata to the out matfile
            save(outName,'-append','z','regData')

        % else due interpolation on each var
        else
            % load x and y coordintes
            x = m.x;
            y = m.y;
            stripList = m.stripList;
            version = m.version;

            % save to matfile
            save(outName,'regData','x','y','stripList','version','-v7.3');

            % register z and add to mat file
            z = applyRegistration(regData.reg.p,m,C,'gridvar','z','subsetSize',5000);
            save(outName,'-append','z');
            clear z

            % register z_mad and add to mat file
            z_mad =  applyRegistration(regData.reg.p,m,C,'gridvar','z_mad','subsetSize',5000);
            save(outName,'-append','z_mad');
            clear z_mad

            % register N and add to mat file
            N =  applyRegistration(regData.reg.p,m,C,'gridvar','N','subsetSize',5000);
            save(outName,'-append','N');
            clear N

            % register Nmt and add to mat file
            Nmt =  applyRegistration(regData.reg.p,m,C,'gridvar','Nmt','subsetSize',5000);
            save(outName,'-append','Nmt');
            clear Nmt

            % register tmin and add to mat file
            tmin =  applyRegistration(regData.reg.p,m,C,'gridvar','tmin','subsetSize',5000);
            save(outName,'-append','tmin');
            clear tmin

            % register tmax and add to mat file
            tmax = applyRegistration(regData.reg.p,m,C,'gridvar','tmax','subsetSize',5000);
            save(outName,'-append','tmax');
            clear tmax

        end
    end
end




