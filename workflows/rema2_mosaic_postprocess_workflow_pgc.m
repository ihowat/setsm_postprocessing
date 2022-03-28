
%% Undo registration and get back to cleanest state of tile .mat files (as if they're fresh out of MST step)

% Remove all *.tif and *meta.txt mosaic results files

% Remove all *_reg.mat mosaic tile files (we will be redoing registration)

% Undo tile buffer merge on all existing 10m and 2m (non-reg) tile .mat files
file_listing = dir(['output_tiles_nga-asap/*/','/*.mat'])
file_paths_cellarr = fullfile({file_listing.folder}.', {file_listing.name}.');
undoMergeTileBuffers(file_paths_cellarr)



%% Should only need to be done once for each new tile out of MST, could be added into MST step

% Add land mask rasters to only the 10m tile .mat files
addLandMask2REMATiles('output_tiles_testing_redo-reg/*/', '/mnt/pgc/data/projects/earthdem/tiledef_files/rema_tile_definitions_plus_sgssi.mat')



%% Register tiles to ICESat-2 and calcualte tile boundary adjustments

% Create *_reg.mat copies of the 10m tile .mat files with registration to ICESat-2.
% - Registration info is calculated and stored in the "reg" variable of both registered and unregistrered .mat tile files.
% - *_reg.mat files are created by the "applyRegistration" function with the 3D vector translation from registration applied to all data arrays.
% - The "fit2is2" function adds quadratic surface "sf" (fit equation) and "dzfit" (mesh grid) variables of point GCP offsets from ICESat-2, and applies offset to the "z" array (DEM data).
batchRegisterTiles('output_tiles_testing_redo-reg/*/', '/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/from_unity/altimetryByTile/')

% Create "tileNeighborIndex_10m.mat" index of 10m tile .mat files (*_reg.mat if exists, else unreg .mat), used in boundary/buffer adjustment steps
tileNeighborIndex('output_tiles_testing_redo-reg', 'resolution', '10m')

% Calculate DEM offset between tile and all its neighbors for all 10m tiles, save as 'dz0' variable in 10m tile .mat files
batch_boundaryAdjustCalc('output_tiles_testing_redo-reg/tileNeighborIndex_10m.mat')



%% Finish 10m tile boundary adjustments

% Simply apply 'z = z - dz0' to 10m tile .mat files ('dz0' calculated in the previous step)
batch_boundaryAdjustApply('output_tiles_testing_redo-reg/tileNeighborIndex_10m.mat')

%%%%%%%%%%%%%%%%%%%%%%%
%%%% LEFT OFF HERE %%%%
%%%%%%%%%%%%%%%%%%%%%%%

% Feather-merge 10m tile buffers
batchMergeTileBuffer('output_tiles_testing_redo-reg/tileNeighborIndex_10m.mat')



%% Finish 2m tile boundary adjustments

% Apply same ICESat-2 registration vector from 10m tile to all 2m quad tiles
applyRegTo2m('output_tiles_testing_redo-reg/tileNeighborIndex_10m.mat')

% Apply both 'dzfit' (ICESat-2 reg residual surface) and 'dz0' (tile neighbor offsets) from 10m tile to 2m quad tiles
applyDzfitTo2m('output_tiles_testing_redo-reg/tileNeighborIndex_10m.mat')

% Create "tileNeighborIndex_2m.mat" index of 2m tile .mat files (*_reg.mat if exists, else unreg .mat)
tileNeighborIndex('output_tiles_testing_redo-reg', 'resolution', '2m')

% Feather-merge 2m tile buffers
batchMergeTileBuffer('output_tiles_testing_redo-reg/tileNeighborIndex_2m.mat')


