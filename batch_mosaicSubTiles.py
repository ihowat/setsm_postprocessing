import os, string, sys, argparse, glob, subprocess
from collections import namedtuple
matlab_scripts = '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing4'
quads = ['1_1','1_2','2_1','2_2']

Task = namedtuple('Task', 't st')

project_choices = [
    'arcticdem',
    'rema',
    'earthdem',
]

tileDefFile_utm_north = 'PGC_UTM_Mosaic_Tiles_North.mat'
tileDefFile_utm_south = 'PGC_UTM_Mosaic_Tiles_South.mat'
tileDefFile_utm_options = "{} or {}".format(tileDefFile_utm_north, tileDefFile_utm_south)
project_tileDefFile_dict = {
    'arcticdem': 'PGC_Imagery_Mosaic_Tiles_Arctic.mat',
    'rema': 'PGC_Imagery_Mosaic_Tiles_Antarctic.mat',
    'earthdem': tileDefFile_utm_options,
}

def main():

    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("srcdir", help="source root dir (level above tile name dir)")
    parser.add_argument("tiles", help="list of tiles, comma delimited")
    parser.add_argument("res", choices=['2','10'], help="resolution (2 or 10)")

    parser.add_argument("--lib-path", default=matlab_scripts,
            help="path to referenced Matlab functions (default={}".format(matlab_scripts))
    parser.add_argument("--project", default=None, choices=project_choices,
                        help="sets the default value of project-specific arguments")
    parser.add_argument("--tile-def", default=None,
                        help="mosaic tile definition mat file (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_tileDefFile_dict.items()])
                        ))
    parser.add_argument('--quads', action='store_true', default=False,
            help="build into quad subtiles")
    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--qsubscript",
            help="qsub script to use in PBS submission (default is qsub_mosaicSubTiles.sh in script root folder)")
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')

    args = parser.parse_args()

    tiles = args.tiles.split(',')
    srcdir = os.path.abspath(args.srcdir)
    scriptdir = os.path.abspath(os.path.dirname(sys.argv[0]))

    matlab_script = 'mosaicSubTiles'

    ## Set default arguments by project setting
    if args.project is None and True in [arg is None for arg in [args.tile_def]]:
        parser.error("--project arg must be provided if one of the following arguments is not provided: {}".format(
            ' '.join(["--tile-def"])
        ))
    if args.tile_def is None:
        args.tile_def = project_tileDefFile_dict[args.project]

    ## Verify path arguments
    if not os.path.isdir(srcdir):
        parser.error("srcdir does not exist: {}".format(srcdir))
    if not os.path.isdir(args.lib_path):
        parser.error("--lib-path does not exist: {}".format(args.lib_path))
    if args.project == 'earthdem':
        pass
    else:
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

    tasks = []

    i=0
    if len(tiles) > 0:

        for tile in tiles:
            if args.quads:
                for quad in quads:
                    tasks.append(Task(tile, quad))
            else:
                tasks.append(Task(tile, ''))


    print("{} tasks found".format(len(tasks)))

    if len(tasks) > 0:
        for task in tasks:

            tile = task.t

            tile_def = args.tile_def

            if tile_def == tileDefFile_utm_options:
                assert args.project == 'earthdem'

                utm_tilename_parts = tile.split('_')
                utm_tilename_prefix = utm_tilename_parts[0]
                if not utm_tilename_prefix.startswith('utm'):
                    parser.error("Expected only UTM tile names (e.g. 'utm10n_01_01'), but got '{}'".format(tile))

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

            if task.st == '':
                dstfn = "{}_{}m.mat".format(task.t,args.res)
            else:
                dstfn = "{}_{}_{}m.mat".format(task.t,task.st,args.res)
            dstfp = os.path.join(srcdir, task.t, dstfn)
            sem = os.path.join(srcdir, task.t, dstfn.replace('.mat','_empty.txt'))
            subtile_dir = os.path.join(srcdir,task.t,'subtiles')

            if os.path.isfile(dstfp):
                print('Output exists, skipping {}'.format(dstfn))
            elif os.path.isfile(sem):
                print('N array was empty on last run, skipping {}'.format(dstfn))

            else:
                ## if pbs, submit to scheduler
                i+=1
                if args.pbs:
                    job_name = 'mst_{}'.format(task.t)
                    cmd = r'qsub -N {1} -v p1={2},p2={3},p3={4},p4={5},p5={6},p6={7},p7={8},p8={9},p9={10} {0}'.format(
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
                    )
                    print(cmd)
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)

                ## else run matlab
                else:
                    if task.st == '':
                        cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); [x0,x1,y0,y1]=getTileExtents('{6}','{7}'); projstr=getTileProjection('{7}'); {2}('{3}',{4},'{5}','projection',projstr,'extent',[x0,x1,y0,y1]); exit" """.format(
                            scriptdir,
                            args.lib_path,
                            matlab_script,
                            subtile_dir,
                            args.res,
                            dstfp,
                            task.t,
                            tile_def
                        )
                    else:
                        cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); [x0,x1,y0,y1]=getTileExtents('{7}','{8}','quadrant','{6}'); projstr=getTileProjection('{8}'); {2}('{3}',{4},'{5}','projection',projstr,'quadrant','{6}','extent',[x0,x1,y0,y1]); exit" """.format(
                            scriptdir,
                            args.lib_path,
                            matlab_script,
                            subtile_dir,
                            args.res,
                            dstfp,
                            task.st,
                            task.t,
                            tile_def
                        )
                    print("{}, {}".format(i, cmd))
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()
