import os, string, sys, argparse, glob, subprocess

SCRIPT_FILE = os.path.abspath(os.path.realpath(__file__))
SCRIPT_FNAME = os.path.basename(SCRIPT_FILE)
SCRIPT_NAME, SCRIPT_EXT = os.path.splitext(SCRIPT_FNAME)
SCRIPT_DIR = os.path.dirname(SCRIPT_FILE)

matlab_scripts = os.path.join(SCRIPT_DIR, '../setsm_postprocessing4')

project_choices = [
    'arcticdem',
    'rema',
    'earthdem',
]

epsg_projstr_dict = {
    3413: 'polar stereo north',
    3031: 'polar stereo south',
}
project_epsg_dict = {
    'arcticdem': 3413,
    'earthdem': None,
    'rema': 3031,
}

earthdem_tileprefix_key = '<tileprefix>'
earthdem_ref_dem_template = '/mnt/pgc/data/elev/dem/tandem-x/90m/TanDEM-X_UTM_90m/TDX_UTM_Mosaic_{}_90m.tif'.format(earthdem_tileprefix_key)

project_refDemFile_dict = {
    'arcticdem': None,
    'earthdem': earthdem_ref_dem_template,
    'rema': '/mnt/pgc/data/elev/dem/tandem-x/90m/TanDEM-X_Antarctica_90m/TanDEM_Antarctica_Mosaic.tif',
}

tileDefFile_utm_north = 'PGC_UTM_Mosaic_Tiles_North.mat'
tileDefFile_utm_south = 'PGC_UTM_Mosaic_Tiles_South.mat'
tileDefFile_utm_options = "{} or {}".format(tileDefFile_utm_north, tileDefFile_utm_south)
project_tileDefFile_dict = {
    'arcticdem': 'PGC_Imagery_Mosaic_Tiles_Arctic.mat',
    # 'rema': 'PGC_Imagery_Mosaic_Tiles_Antarctic.mat',
    'rema': 'rema_tile_definitions.mat',
    'earthdem': tileDefFile_utm_options,
}

project_databaseFile_dict = {
    'arcticdem': '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing_pgc/ArcticDEMdatabase4_2m_v4_20201218.mat',
    # 'rema': '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing_pgc/REMAdatabase4_2m_v4_20200806.mat',
    'rema': '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing_pgc/rema_strips_v13e.shp',
    'earthdem': '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing_pgc/EarthDEMdatabase4_2m_v4_20210101.mat',
}
project_waterTileDir_dict = {
    'arcticdem': '/mnt/pgc/data/projects/arcticdem/watermasks/global_surface_water/tiled_watermasks/',
    'rema':      '',
    'earthdem':  '/mnt/pgc/data/projects/earthdem/watermasks/global_surface_water/tiled_watermasks/',
}
project_stripsDirectory_dict = {
    'arcticdem': '/mnt/pgc/data/elev/dem/setsm/ArcticDEM/region',
    'rema':      '/mnt/pgc/data/elev/dem/setsm/REMA/region',
    'earthdem':  '/mnt/pgc/data/elev/dem/setsm/EarthDEM/region',
}
project_tileqcDir_dict = {
    'arcticdem': '',
    'rema':      '/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/tile_qc_v13e',
    'earthdem':  '',
}
project_tileParamList_dict = {
    'arcticdem': '',
    'rema':      '/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/tileParamList_v13e.txt',
    'earthdem':  '',
}


def main():

    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("dstdir", help="target directory (tile subfolders will be created)")
    parser.add_argument("tiles", help="list of mosaic tiles, comma delimited")

    parser.add_argument("--project", default=None, choices=project_choices,
                        help="sets the default value of project-specific arguments")
    parser.add_argument("--ref-dem", default=None, help="reference DEM (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_refDemFile_dict.items()])
                        ))
    parser.add_argument("--epsg", type=int, default=None, choices=list(epsg_projstr_dict.keys()),
                        help="output mosaic tile projection EPSG code (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_epsg_dict.items()])
                        ))
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
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_waterTileDir_dict.items()])
                        ))
    parser.add_argument("--tileqc-dir", default=None,
                        help="directory of tile qc mat files (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_tileqcDir_dict.items()])
                        ))
    parser.add_argument("--lib-path", default=matlab_scripts,
                        help="path to referenced Matlab functions (default={})".format(matlab_scripts))

    # parser.add_argument('--require-finfiles', action='store_true', default=False,
    #         help="let existence of finfiles dictate reruns")
    parser.add_argument('--bypass-finfile-req', action='store_true', default=False,
            help="upon rerun, deem tiles complete when the 10,000-th subtile exists even though finfile does not exist")

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
    scriptdir = SCRIPT_DIR

    matlab_script = 'buildSubTiles'
    if args.sort_fix:
        matlab_script = 'buildSubTilesSortFix'

    make2m_arg = 'false' if args.make_10m_only else 'true'

    ## Set default arguments by project setting
    if args.project is None and True in [arg is None for arg in
            [args.ref_dem, args.epsg, args.tile_def, args.strip_db, args.water_tile_dir]]:
        parser.error("--project arg must be provided if one of the following arguments is not provided: {}".format(
            ' '.join(["--ref-dem", "--epsg", "--tile-def", "--strip-db", "--water-tile-dir"])
        ))
    if args.epsg is None:
        args.epsg = project_epsg_dict[args.project]
    if args.ref_dem is None:
        args.ref_dem = project_refDemFile_dict[args.project]
        if args.ref_dem is None:
            parser.error("--ref-dem argument must be provided if --project={}".format(args.project))
    if args.tile_def is None:
        args.tile_def = project_tileDefFile_dict[args.project]
    if args.strip_db is None:
        args.strip_db = project_databaseFile_dict[args.project]
    if args.water_tile_dir is None:
        args.water_tile_dir = project_waterTileDir_dict[args.project]
    if args.project is not None:
        if args.tileqc_dir is None:
            args.tileqc_dir = project_tileqcDir_dict[args.project]
        strips_dir = project_stripsDirectory_dict[args.project]
        tileparam_file = project_tileParamList_dict[args.project]
    else:
        args.tileqc_dir = ''
        strips_dir = ''
        tileparam_file = ''

    ## Convert argument paths to absolute paths
    if args.water_tile_dir != '':
        args.water_tile_dir = os.path.abspath(args.water_tile_dir)
    if strips_dir != '':
        strips_dir = os.path.abspath(strips_dir)
    if args.tileqc_dir != '':
        args.tileqc_dir = os.path.abspath(args.tileqc_dir)
    if tileparam_file != '':
        tileparam_file = os.path.abspath(tileparam_file)

    ## Verify path arguments
    if not os.path.isdir(dstdir):
        parser.error("srcdir does not exist: {}".format(dstdir))
    if args.project == 'earthdem':
        projection_string = None
    else:
        if args.ref_dem is not None and not os.path.isfile(args.ref_dem):
            parser.error("--ref-dem does not exist: {}".format(args.ref_dem))
        projection_string = epsg_projstr_dict[args.epsg]
        tile_def_abs = os.path.abspath(args.tile_def if os.path.isfile(args.tile_def) else os.path.join(scriptdir, args.tile_def))
        if not os.path.isfile(tile_def_abs):
            parser.error("--tile-def file does not exit: {}".format(tile_def_abs))
    strip_db_abs = os.path.abspath(args.strip_db if os.path.isfile(args.strip_db) else os.path.join(scriptdir, args.strip_db))
    if not os.path.isfile(strip_db_abs):
        parser.error("--strip-db does not exist: {}".format(strip_db_abs))
    if args.water_tile_dir != '' and not os.path.isdir(args.water_tile_dir):
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

    num_tiles_to_run = 0

    i=0
    if len(tiles) > 0:

        for tile in tiles:

            ref_dem = args.ref_dem
            tile_projstr = projection_string
            tile_def = args.tile_def

            if ref_dem == earthdem_ref_dem_template or tile_projstr is None or tile_def == tileDefFile_utm_options:
                assert args.project == 'earthdem'

                utm_tilename_parts = tile.split('_')
                utm_tilename_prefix = utm_tilename_parts[0]
                if not utm_tilename_prefix.startswith('utm'):
                    parser.error("Expected only UTM tile names (e.g. 'utm10n_01_01'), but got '{}'".format(tile))

                if tile_projstr is None:
                    tile_projstr = utm_tilename_prefix

                if ref_dem == earthdem_ref_dem_template:
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

            final_subtile_fp_10m = os.path.join(tile_dstdir,'{}_10000_10m.mat'.format(tile))
            final_subtile_fp_2m = os.path.join(tile_dstdir,'{}_10000_2m.mat'.format(tile))
            if args.make_10m_only:
                final_subtile_fp = final_subtile_fp_10m
            else:
                final_subtile_fp = final_subtile_fp_2m
            finfile = final_subtile_fp.replace('.mat', '.fin')
            finfile_2m = final_subtile_fp_2m.replace('.mat', '.fin')

            run_tile = True
            if args.sort_fix or args.rerun_without_cleanup:
                pass
            elif args.rerun:
                print('Verifying tile {} before rerun'.format(tile))

                # if os.path.isfile(final_subtile_fp) and not args.require_finfiles:
                if os.path.isfile(final_subtile_fp) and args.bypass_finfile_req:
                    print('Tile seems complete ({} exists)'.format(os.path.basename(final_subtile_fp)))
                    run_tile = False
                elif os.path.isfile(finfile):
                    print('Tile seems complete ({} exists)'.format(os.path.basename(finfile)))
                    run_tile = False
                elif os.path.isfile(finfile_2m):
                    print('Tile seems complete (2m finfile {} exists)'.format(os.path.basename(finfile_2m)))
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
                num_tiles_to_run += 1

                matlab_cmd = "matlab -nojvm -nodisplay -nosplash -r \\\"try; " \
                    "addpath('{0}'); addpath('{1}'); " \
                    "[meta,landtile] = initializeMosaic('','{3}',''," \
                        "'projection','{10}','tileDefFile','{5}'," \
                        "'stripDatabaseFile','{6}','stripsDirectory','{11}'," \
                        "'waterTileDir','{7}','refDemFile','{8}'," \
                        "'tileqcDir','{12}','tileParamListFile','{13}'," \
                        "'returnMetaOnly',true); " \
                    "{2}('{3}','{4}','{5}',meta," \
                        "'landTile',landtile," \
                        "'refDemFile','{8}'," \
                        "'make2m',{9}," \
                        "'projection','{10}'); " \
                    "catch e; disp(getReport(e)); exit(1); end; exit(0);\\\"".format(
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
                    tile_projstr,
                    strips_dir,
                    args.tileqc_dir,
                    tileparam_file,
                )

                ## if pbs, submit to scheduler
                i+=1
                if args.pbs:
                    job_name = 'bst_{}'.format(tile)
                    cmd = r"""qsub -N {1} -v task_cmd="{2}",finfile="{3}" {0}""".format(
                        qsubpath,
                        job_name,
                        matlab_cmd.replace(',', '|COMMA|'),
                        finfile,
                    )
                    print(cmd)
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)

                ## else run matlab
                else:
                    print("{}, {}".format(i, matlab_cmd))
                    if not args.dryrun:
                        subprocess.call(matlab_cmd, shell=True)

    print("Running {} tiles".format(num_tiles_to_run))

    print("Done")


if __name__ == '__main__':
    main()
