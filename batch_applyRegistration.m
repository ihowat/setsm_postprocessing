function batch_applyRegistration(remaTileDir,registrationFile,varargin)
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
%remaTileDir='~/project/howat.4/rema_mosaic_v13e';

%remaTileDir='/fs/project/howat.4/howat.4/rema_mosaic/rema_mosaic_v13e';

% region name (which is the name of the region directory). If empty, does
% not use region directory (i.e. if all files are in one folder)
% regionName='rema_21_mbl_north';
%regionName='';

% File with list of registration parameters, indexed by tile name
% if empty, uses reg variable in tile mat
%registrationFile=  /fs/byo/howat-data4/REMA/altimetryByTile/rema_{region number}_{region name}_is2reg.mat
%registrationFile= '~/data4/REMA/altimetryByTile/rema_21_mbl_north_is2reg.mat';
%registrationFile=[];

strt=1;
inc=1;
if length(varargin) == 2
    strt=varargin{1};
    inc= varargin{2};
end
fprintf('processing every %d tile, starting at %d\n',inc,strt)


%% Read file list
demMatFiles=dir([remaTileDir,'/','*_10m.mat']);
demMatFiles=cellfun( @(x) [remaTileDir,'/',x],...
    {demMatFiles.name}, 'uniformoutput',0);

fprintf('%d tiles found\n',length(demMatFiles))

%% loop through tiles
for i=strt:inc:length(demMatFiles)
    
    % File to register
    demMatFile=demMatFiles{i};
    
    applyRegistration(demMatFile,registrationFile)
    
end




