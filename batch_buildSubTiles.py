import os, string, sys, argparse, glob, subprocess
matlab_scripts = '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing4'

project_choices = [
    'arcticdem',
    'rema',
    'earthdem',
]

earthdem_tileprefix_key = '<tileprefix>'
earthdem_ref_dem_template = '/mnt/pgc/data/elev/dem/tandem-x/90m/TanDEM-X_UTM_90m/TDX_UTM_Mosaic_{}_90m.tif'.format(earthdem_tileprefix_key)

tileDefFile_utm_north = 'PGC_UTM_Mosaic_Tiles_North.mat'
tileDefFile_utm_south = 'PGC_UTM_Mosaic_Tiles_South.mat'
tileDefFile_utm_options = "{} or {}".format(tileDefFile_utm_north, tileDefFile_utm_south)
project_tileDefFile_dict = {
    'arcticdem': 'PGC_Imagery_Mosaic_Tiles_Arctic.mat',
    'rema': 'PGC_Imagery_Mosaic_Tiles_Antarctic.mat',
    'earthdem': tileDefFile_utm_options,
}

project_databaseFile_dict = {
    'arcticdem': 'arcticDEMdatabase4_2m_v4_20200806.mat',
    'rema': 'REMAdatabase4_2m_v4_20200806.mat',
    'earthdem': 'EarthDEMdatabase4_2m_v4_20201012_combined_fixed.mat',
}
waterTileDir_dict = {
    'arcticdem': '/mnt/pgc/data/projects/arcticdem/watermasks/global_surface_water/tiled_watermasks/',
    'rema':      '/mnt/pgc/data/projects/rema/watermasks/global_surface_water/tiled_watermasks/',
    'earthdem':  '/mnt/pgc/data/projects/earthdem/watermasks/global_surface_water/tiled_watermasks/',
}

def main():

    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("dstdir", help="target directory (tile subfolders will be created)")
    parser.add_argument("tiles", help="list of mosaic tiles, comma delimited")

    parser.add_argument("--ref-dem", default=None, help="reference DEM (required for ArcticDEM & REMA, automatically selected "
                        "for EarthDEM by file path template {})".format(earthdem_ref_dem_template))
    parser.add_argument("--project", default=None, choices=project_choices,
                        help="sets the default value of project-specific arguments")
    parser.add_argument("--tile-def", default=None,
                        help="mosaic tile definition mat file (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_tileDefFile_dict.items()])
                        ))
    parser.add_argument("--strip-db", default=None,
                        help="strip database mat file (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_databaseFile_dict.items()])
                        ))
    parser.add_argument("--water-tile-dir", default=None,
                        help="directory of water tifs (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in waterTileDir_dict.items()])
                        ))
    parser.add_argument("--lib-path", default=matlab_scripts,
                        help="path to referenced Matlab functions (default={})".format(matlab_scripts))

    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--rerun", action='store_true', default=False,
            help="rerun tile, behavior determined by redoFlag in Matlab code")
    parser.add_argument("--rerun-without-cleanup", action='store_true', default=False,
            help="rerun tile without attempting to clean up potentially corrupted subtiles")
    parser.add_argument("--sort-fix", action="store_true", default=False,
            help="run tile with buildSubTilesSortFix script")
    parser.add_argument("--make-10m-only", action='store_true', default=False,
            help="do not give 'make2m' argument in call to buildSubTiles.m, breaks --rerun capability")
    parser.add_argument("--qsubscript",
            help="qsub script to use in PBS submission (default is qsub_buildSubTiles.sh in script root folder)")
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')

    args = parser.parse_args()

    tiles = args.tiles.split(',')
    dstdir = os.path.abspath(args.dstdir)
    scriptdir = os.path.abspath(os.path.dirname(sys.argv[0]))

    matlab_script = 'buildSubTiles'
    if args.sort_fix:
        matlab_script = 'buildSubTilesSortFix'

    make2m_arg = 'false' if args.make_10m_only else 'true'

    ## Set default arguments by project setting
    if args.project is None and True in [arg is None for arg in [args.tile_def, args.strip_db, args.water_tile_dir]]:
        parser.error("--project arg must be provided if one of the following arguments is not provided: {}".format(
            ' '.join(["--tile-def", "--strip-db", "--water-tile-dir"])
        ))
    if args.ref_dem is None and args.project != 'earthdem':
        parser.error("--ref-dem argument must be provided if not --project=earthdem")
    if args.tile_def is None:
        args.tile_def = project_tileDefFile_dict[args.project]
    if args.strip_db is None:
        args.strip_db = project_databaseFile_dict[args.project]
    if args.water_tile_dir is None:
        args.water_tile_dir = waterTileDir_dict[args.project]

    ## Verify path arguments
    if not os.path.isdir(dstdir):
        parser.error("srcdir does not exist: {}".format(dstdir))
    if args.ref_dem is not None and not os.path.isfile(args.ref_dem):
        parser.error("--ref-dem does not exist: {}".format(args.ref_dem))
    if args.project == 'earthdem':
        pass
    else:
        tile_def_abs = os.path.join(scriptdir, args.tile_def)
        if not os.path.isfile(tile_def_abs):
            parser.error("--tile-def file does not exit: {}".format(tile_def_abs))
    strip_db_abs = os.path.join(scriptdir, args.strip_db)
    if not os.path.isfile(strip_db_abs):
        parser.error("--strip-db does not exist: {}".format(strip_db_abs))
    if not os.path.isdir(args.water_tile_dir):
        parser.error("--water-tile-dir does not exist: {}".format(args.water_tile_dir))
    if not os.path.isdir(args.lib_path):
        parser.error("--lib-path does not exist: {}".format(args.lib_path))

    if args.rerun and args.rerun_without_cleanup:
        parser.error("--rerun and --rerun-without-cleanup are mutually exclusive")

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,'qsub_buildSubTiles.sh')
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)

    i=0
    if len(tiles) > 0:

        for tile in tiles:

            ref_dem = args.ref_dem
            tile_def = args.tile_def

            if ref_dem is None or tile_def == tileDefFile_utm_options:
                assert args.project == 'earthdem'

                utm_tilename_parts = tile.split('_')
                utm_tilename_prefix = utm_tilename_parts[0]
                if not utm_tilename_prefix.startswith('utm'):
                    parser.error("Expected only UTM tile names (e.g. 'utm10n_01_01'), but got '{}'".format(tile))

                if ref_dem is None:
                    ref_dem = earthdem_ref_dem_template.replace(earthdem_tileprefix_key, utm_tilename_prefix)

                if tile_def == tileDefFile_utm_options:
                    if utm_tilename_prefix.endswith('n'):
                        tile_def = tileDefFile_utm_north
                    elif utm_tilename_prefix.endswith('s'):
                        tile_def = tileDefFile_utm_south
                    else:
                        parser.error("UTM tile name prefix does not end with 'n' or 's' (e.g. 'utm10n'): {}".format(tile))

                tile_def_abs = os.path.join(scriptdir, tile_def)
                if not os.path.isfile(tile_def_abs):
                    parser.error("tile def file does not exit: {}".format(tile_def_abs))

            ## if output does not exist, add to task list
            tile_dstdir = os.path.join(dstdir,tile,'subtiles')
            if not os.path.isdir(tile_dstdir):
                if not args.dryrun:
                    os.makedirs(tile_dstdir)
            dstfps = glob.glob(os.path.join(tile_dstdir,'{}_*m.mat'.format(tile)))

            if args.make_10m_only:
                final_subtile_fp = os.path.join(tile_dstdir,'{}_10000_10m.mat'.format(tile))
            else:
                final_subtile_fp = os.path.join(tile_dstdir,'{}_10000_2m.mat'.format(tile))
            finfile = final_subtile_fp.replace('.mat', '.fin')

            run_tile = True
            if args.sort_fix or args.rerun_without_cleanup:
                pass
            elif args.rerun:
                print('Verifying tile {} before rerun'.format(tile))

                if os.path.isfile(final_subtile_fp):
                    print('Tile seems complete ({} exists)'.format(os.path.basename(final_subtile_fp)))
                    run_tile = False
                elif os.path.isfile(finfile):
                    print('Tile seems complete ({} exists)'.format(os.path.basename(finfile)))
                    run_tile = False

                if not args.make_10m_only:
                    ## clean up subtiles with only 10m version
                    dstfps_10m = glob.glob(os.path.join(tile_dstdir,'{}_*10m.mat'.format(tile)))
                    for dstfp_10m in dstfps_10m:
                        dstfp_2m = dstfp_10m.replace('10m.mat','2m.mat')
                        if not os.path.isfile(dstfp_2m):
                            print('Removing 10m subtile missing 2m component: {}'.format(os.path.basename(dstfp_10m)))
                            run_tile = True
                            if not args.dryrun:
                                os.remove(dstfp_10m)
            
            elif len(dstfps) > 0:
                print('{} subtiles exists, skipping'.format(tile))
                run_tile = False
            
            if run_tile:
                ## if pbs, submit to scheduler
                i+=1
                if args.pbs:
                    job_name = 'bst_{}'.format(tile)
                    cmd = r'qsub -N {1} -v p1={2},p2={3},p3={4},p4={5},p5={6},p6={7},p7={8},p8={9},p9={10},p10={11},p11={12} {0}'.format(
                        qsubpath,
                        job_name,
                        scriptdir,  #p1
                        args.lib_path,  #p2
                        matlab_script, #p3
                        tile,  #p4
                        tile_dstdir, #p5
                        tile_def,  #p6
                        args.strip_db,  #p7
                        args.water_tile_dir,  #p8
                        ref_dem,  #p9
                        make2m_arg,  #p10
                        finfile,  #p11
                    )
                    print(cmd)
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)

                ## else run matlab
                else:
                    cmd = """try; matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); {2}('{3}','{4}','{5}','{6}','landTile','{7}','refDemFile','{8}','make2m',{9}); catch e; disp(getReport(e)); exit(1); end; exit(0);" """.format(
                        scriptdir,
                        args.lib_path,
                        matlab_script,
                        tile,
                        tile_dstdir,
                        tile_def,
                        args.strip_db,
                        args.water_tile_dir,
                        ref_dem,
                        make2m_arg,
                    )
                    print("{}, {}".format(i, cmd))
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()
