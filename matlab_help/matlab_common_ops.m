

>> file_listing = dir([sourceDir,'/*/*.mat'])

fileNames =

  128x1 struct array with fields:

    name
    folder
    date
    bytes
    isdir
    datenum


>> file_paths_cellarr = fullfile({file_listing.folder}.', {file_listing.name}.')

file_paths_cellarr =

  128x1 cell array

    {'/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/results/output_tiles_testing_redo-reg/31_31/31_31_10m.mat'   }
    {'/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/results/output_tiles_testing_redo-reg/31_31/31_31_2_2_2m.mat'}
    {'/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/results/output_tiles_testing_redo-reg/31_32/31_32_10m.mat'   }
    {'/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/results/output_tiles_testing_redo-reg/31_32/31_32_1_2_2m.mat'}
    {'/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/results/output_tiles_testing_redo-reg/31_32/31_32_2_1_2m.mat'}
    {'/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/results/output_tiles_testing_redo-reg/31_32/31_32_2_2_2m.mat'}
    ...



