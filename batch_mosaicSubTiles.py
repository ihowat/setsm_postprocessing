import os, string, sys, argparse, glob, subprocess
from collections import namedtuple

SCRIPT_FILE = os.path.abspath(os.path.realpath(__file__))
SCRIPT_FNAME = os.path.basename(SCRIPT_FILE)
SCRIPT_NAME, SCRIPT_EXT = os.path.splitext(SCRIPT_FNAME)
SCRIPT_DIR = os.path.dirname(SCRIPT_FILE)

matlab_scripts = os.path.join(SCRIPT_DIR, '../setsm_postprocessing4')

quads = ['1_1','1_2','2_1','2_2']

Task = namedtuple('Task', 't st')

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

tileDefFile_utm_north = 'PGC_UTM_Mosaic_Tiles_North.mat'
tileDefFile_utm_south = 'PGC_UTM_Mosaic_Tiles_South.mat'
tileDefFile_utm_options = "{} or {}".format(tileDefFile_utm_north, tileDefFile_utm_south)
project_tileDefFile_dict = {
    'arcticdem': 'PGC_Imagery_Mosaic_Tiles_Arctic.mat',
    # 'rema': 'PGC_Imagery_Mosaic_Tiles_Antarctic.mat',
    'rema': 'rema_tile_definitions.mat',
    'earthdem': tileDefFile_utm_options,
}
project_version_dict = {
    'arcticdem': 'ArcticDEM|4.1',
    'rema': 'REMA|2.0',
    'earthdem': 'EarthDEM|1.0',
}

def main():

    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("srcdir", help="source root dir (level above tile name dir)")
    parser.add_argument("tiles", help="list of tiles, comma delimited")
    parser.add_argument("res", type=int, choices=[2, 10], help="resolution (2 or 10)")

    parser.add_argument("--lib-path", default=matlab_scripts,
                        help="path to referenced Matlab functions (default={}".format(matlab_scripts))
    parser.add_argument("--project", default=None, choices=project_choices,
                        help="sets the default value of project-specific arguments")
    parser.add_argument("--epsg", type=int, default=None, choices=list(epsg_projstr_dict.keys()),
                        help="output mosaic tile projection EPSG code (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_epsg_dict.items()])
                        ))
    parser.add_argument("--tile-def", default=None,
                        help="mosaic tile definition mat file (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_tileDefFile_dict.items()])
                        ))
    parser.add_argument("--version", default=None,
                        help="mosaic version (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_version_dict.items()])
                        ))
    parser.add_argument('--quads', action='store_true', default=False,
            help="build into quad subtiles")

    parser.add_argument('--bypass-bst-finfile-req', action='store_true', default=False,
            help="do not require BST finfiles exist before mosaicking tiles")
    parser.add_argument('--relax-bst-finfile-req', action='store_true', default=False,
            help="allow mosaicking tiles with no BST finfile if 10,000-th subtile exists")
    # parser.add_argument('--require-mst-finfiles', action='store_true', default=False,
    #         help="let existence of MST finfiles dictate reruns")
    parser.add_argument('--bypass-mst-finfile-req', action='store_true', default=False,
            help="upon rerun, deem tiles complete when MST results exist even though MST finfile does not exist")
    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--qsubscript",
            help="qsub script to use in PBS submission (default is qsub_mosaicSubTiles.sh in script root folder)")
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')

    args = parser.parse_args()

    if os.path.isfile(args.tiles):
        tilelist_file = args.tiles
        with open(tilelist_file, 'r') as tilelist_fp:
            tiles = tilelist_fp.read().splitlines()
    else:
        tiles = args.tiles.split(',')
    tiles = sorted(list(set(tiles)))

    srcdir = os.path.abspath(args.srcdir)
    scriptdir = SCRIPT_DIR

    matlab_script = 'mosaicSubTiles'

    ## Set default arguments by project setting
    if args.project is None and True in [arg is None for arg in [args.epsg, args.tile_def, args.version]]:
        parser.error("--project arg must be provided if one of the following arguments is not provided: {}".format(
            ' '.join(["--epsg", "--tile-def", "--version"])
        ))
    if args.epsg is None:
        args.epsg = project_epsg_dict[args.project]
    if args.tile_def is None:
        args.tile_def = project_tileDefFile_dict[args.project]
    if args.version is None:
        args.version = project_version_dict[args.project]

    ## Verify path arguments
    if not os.path.isdir(srcdir):
        parser.error("srcdir does not exist: {}".format(srcdir))
    if not os.path.isdir(args.lib_path):
        parser.error("--lib-path does not exist: {}".format(args.lib_path))
    if args.project == 'earthdem':
        projection_string = None
    else:
        projection_string = epsg_projstr_dict[args.epsg]
        tile_def_abs = os.path.join(scriptdir, args.tile_def)
        if not os.path.isfile(tile_def_abs):
            parser.error("--tile-def file does not exit: {}".format(tile_def_abs))

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,'qsub_mosaicSubTiles.sh')
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)

    if args.bypass_bst_finfile_req and args.relax_bst_finfile_req:
        parser.error("--bypass-bst-finfile-req and --relax-bst-finfile-req arguments are mutually exclusive")

    tasks = []
    error_messages = []
    supertile_num_nodata_dict = dict()

    i=0
    if len(tiles) > 0:

        for tile in tiles:
            if args.quads:
                for quad in quads:
                    tasks.append(Task(tile, quad))
            else:
                tasks.append(Task(tile, 'null'))

    print("{} tasks found".format(len(tasks)))
    num_tiles_to_run = 0

    if len(tasks) > 0:
        for task in tasks:

            tile = task.t

            tile_projstr = projection_string
            tile_def = args.tile_def

            if tile_projstr is None or tile_def == tileDefFile_utm_options:
                assert args.project == 'earthdem'

                utm_tilename_parts = tile.split('_')
                utm_tilename_prefix = utm_tilename_parts[0]
                if not utm_tilename_prefix.startswith('utm'):
                    parser.error("Expected only UTM tile names (e.g. 'utm10n_01_01'), but got '{}'".format(tile))

                if tile_projstr is None:
                    tile_projstr = utm_tilename_prefix

                if tile_def == tileDefFile_utm_options:
                    if utm_tilename_prefix.endswith('n'):
                        tile_def = tileDefFile_utm_north
                    elif utm_tilename_prefix.endswith('s'):
                        tile_def = tileDefFile_utm_south
                    else:
                        parser.error("UTM tile name prefix does not end with 'n' or 's' (e.g. 'utm10n'): {}".format(tile))

                tile_def_abs = os.path.join(scriptdir, tile_def)
                if not os.path.isfile(tile_def_abs):
                    parser.error("tile def file does not exist: {}".format(tile_def_abs))

            if task.st == 'null':
                dstfn = "{}_{}m.mat".format(task.t,args.res)
            else:
                dstfn = "{}_{}_{}m.mat".format(task.t,task.st,args.res)
            dstfp = os.path.join(srcdir, task.t, dstfn)
            finfile = os.path.join(srcdir, task.t, dstfn.replace('.mat','.fin'))
            subtile_dir = os.path.join(srcdir,task.t,'subtiles')

            if not os.path.isdir(subtile_dir):
                message = "ERROR! Subtile directory ({}) does not exist, skipping {}".format(subtile_dir, dstfn)
                print(message)
                error_messages.append(message)
                continue

            run_tile = True
            removing_existing_output = False

            mst_finfile = finfile
            bst_final_subtile_fp = os.path.join(subtile_dir, '{}_10000_{}m.mat'.format(task.t, args.res))
            bst_finfile_10m = "{}_10m.fin".format(subtile_dir)
            bst_finfile_2m = "{}_2m.fin".format(subtile_dir)
            if args.res == 10:
                bst_finfile = bst_finfile_10m
            elif args.res == 2:
                bst_finfile = bst_finfile_2m

            if (not args.bypass_bst_finfile_req) and (not any([os.path.isfile(f) for f in [bst_finfile, bst_finfile_2m]])):
                if args.relax_bst_finfile_req and os.path.isfile(bst_final_subtile_fp):
                    message = "WARNING! BST finfile ({}) does not exist for tile {}, but 10,000-th subtile exists so may run".format(bst_finfile, dstfn)
                    print(message)
                else:
                    message = "ERROR! BST finfile ({}) does not exist, skipping {}".format(bst_finfile, dstfn)
                    print(message)
                    error_messages.append(message)
                    if os.path.isfile(bst_final_subtile_fp):
                        message = "  (but 10,000-th subtile exists; can provide --relax-bst-finfile-req argument to run this tile anyways)"
                        print(message)
                        error_messages.append(message)
                    run_tile = False
            else:
                for bst_finfile_temp in list({bst_finfile, bst_finfile_2m, bst_final_subtile_fp}):
                    if os.path.isfile(bst_finfile_temp):
                        message = None
                        if os.path.isfile(mst_finfile) and (os.path.getmtime(bst_finfile_temp) > os.path.getmtime(mst_finfile)):
                            message = "BST finfile ({}) is newer than MST finfile ({}), so existing results will be removed".format(bst_finfile_temp, mst_finfile)
                            removing_existing_output = True
                        elif os.path.isfile(dstfp) and (os.path.getmtime(bst_finfile_temp) > os.path.getmtime(dstfp)):
                            message = "BST finfile ({}) is newer than MST output ({}), so existing results will be removed".format(bst_finfile_temp, dstfp)
                            removing_existing_output = True
                        if message is not None:
                            print(message)
                            error_messages.append(message)

            # if os.path.isfile(dstfp) and args.require_mst_finfiles:
            if os.path.isfile(dstfp) and not args.bypass_mst_finfile_req:
                if not os.path.isfile(finfile):
                    message = "WARNING! MST finfile ({}) does not exist, so existing results will be removed".format(finfile)
                    print(message)
                    error_messages.append(message)
                    removing_existing_output = True

            if removing_existing_output:
                dstfps_old_pattern = dstfp.replace('.mat', '*')
                dstfps_old = glob.glob(dstfps_old_pattern)
                if dstfps_old:
                    print("{}Removing old MST results matching {}".format('(dryrun) ' if args.dryrun else '', dstfps_old_pattern))
                    if not args.dryrun:
                        for dstfp_old in dstfps_old:
                            os.remove(dstfp_old)
                run_tile = True

            else:
                # if os.path.isfile(dstfp) and not args.require_mst_finfiles:
                if os.path.isfile(dstfp) and args.bypass_mst_finfile_req:
                    print("Output exists, skipping {}".format(dstfn))
                    run_tile = False
                elif os.path.isfile(finfile):
                    print("finfile exists, skipping {}".format(dstfn))
                    run_tile = False
                    if not os.path.isfile(dstfp):
                        message = "WARNING! MST finfile exists ({}) but expected output does not exist ({}) for tile {}".format(
                            finfile, dstfp, dstfn
                        )
                        # print(message)
                        if task.t not in supertile_num_nodata_dict:
                            supertile_num_nodata_dict[task.t] = 0
                        supertile_num_nodata_dict[task.t] += 1

            if run_tile:
                num_tiles_to_run += 1

                ## if pbs, submit to scheduler
                i+=1
                if args.pbs:
                    job_name = 'mst_{}'.format(task.t)
                    cmd = r"""qsub -N {1} -v p1={2},p2={3},p3={4},p4={5},p5={6},p6={7},p7={8},p8={9},p9={10},p10={11},p11='{12}',p12='{13}' {0}""".format(
                        qsubpath,
                        job_name,
                        scriptdir,
                        args.lib_path,
                        matlab_script,
                        subtile_dir,
                        args.res,
                        dstfp,
                        task.t,
                        tile_def,
                        task.st,
                        finfile,
                        tile_projstr,
                        args.version,
                    )
                    print(cmd)
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)

                ## else run matlab
                else:
                    if task.st == 'null':
                        cmd = """matlab -nojvm -nodisplay -nosplash -r "try; addpath('{0}'); addpath('{1}'); [x0,x1,y0,y1]=getTileExtents('{6}','{7}'); {2}('{3}',{4},'{5}','projection','{8}','version','{9}','extent',[x0,x1,y0,y1]); catch e; disp(getReport(e)); exit(1); end; exit(0);" """.format(
                            scriptdir,
                            args.lib_path,
                            matlab_script,
                            subtile_dir,
                            args.res,
                            dstfp,
                            task.t,
                            tile_def,
                            tile_projstr,
                            args.version,
                        )
                    else:
                        cmd = """matlab -nojvm -nodisplay -nosplash -r "try; addpath('{0}'); addpath('{1}'); [x0,x1,y0,y1]=getTileExtents('{7}','{8}','quadrant','{6}'); {2}('{3}',{4},'{5}','projection','{9}','version','{10}','quadrant','{6}','extent',[x0,x1,y0,y1]); catch e; disp(getReport(e)); exit(1); end; exit(0);" """.format(
                            scriptdir,
                            args.lib_path,
                            matlab_script,
                            subtile_dir,
                            args.res,
                            dstfp,
                            task.st,
                            task.t,
                            tile_def,
                            tile_projstr,
                            args.version,
                        )
                    print("{}, {}".format(i, cmd))
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)

    if error_messages:
        print('----')
        print("The following error messages were received")
        print('----')
        for errmsg in error_messages:
            print(errmsg)
        print('----')

    inspect_tiles = []
    for tile, num_nodata in supertile_num_nodata_dict.items():
        if (args.quads and num_nodata == 4) or (not args.quads and num_nodata > 0):
            inspect_tiles.append(tile)

    if inspect_tiles:
        print('-----')
        print("{} tiles have all MST finfiles but no output mosaic results!".format(len(inspect_tiles)))
        print("Please investigate why these tiles produce no results!!")
        print(','.join(inspect_tiles))
        print('-----')
        print("Checking those {} super-tiles for existence of subtile results...".format(len(inspect_tiles)))
        for tile in inspect_tiles:
            subtile_dir = os.path.join(srcdir,tile,'subtiles')
            if not glob.glob(os.path.join(subtile_dir, '{}_*{}m.mat'.format(tile, args.res))):
                print("ERROR! No {}m results exist in subtile directory for tile {}: {}".format(args.res, tile, subtile_dir))
        print('-----')

    print("Running {} {}tiles".format(num_tiles_to_run, 'quad-' if args.quads else ''))

    print("Done")



if __name__ == '__main__':
    main()
