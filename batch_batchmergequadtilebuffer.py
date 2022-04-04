import os, string, sys, argparse, glob, subprocess

SCRIPT_FILE = os.path.abspath(os.path.realpath(__file__))
SCRIPT_FNAME = os.path.basename(SCRIPT_FILE)
SCRIPT_NAME, SCRIPT_EXT = os.path.splitext(SCRIPT_FNAME)
SCRIPT_DIR = os.path.dirname(SCRIPT_FILE)

matlab_scripts = os.path.join(SCRIPT_DIR, '../setsm_postprocessing4')
quadnames = ('1_1','1_2','2_1','2_2')
qsub_default = 'qsub_mergetilebuffer.sh'
# RESOLUTIONS = ['2','10']
RESOLUTIONS = ['2']

def main():

    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("tiledir", help="target directory")
    parser.add_argument("tiles",
        help=' '.join([
            "list of mosaic supertiles; either specified on command line (comma delimited),",
            "or a text file list (each tile on separate line)"
        ])
    )
    parser.add_argument("res", choices=RESOLUTIONS, help="resolution ({})".format(','.join(RESOLUTIONS)))
    parser.add_argument("dimension", choices=['row','column'], help="dimension on which to group tiles for merging")
    parser.add_argument("--lib-path", default=matlab_scripts,
            help="path to referenced Matlab functions (default={}".format(matlab_scripts))
    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--qsubscript",
            help="qsub script to use in PBS submission (default is {} in script root folder)".format(qsub_default))
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')

    args = parser.parse_args()

    if args.tiles.lower().endswith(('.txt', '.csv')) or os.path.isfile(args.tiles):
        tilelist_file = args.tiles
        if not os.path.isfile(args.tiles):
            parser.error("'tiles' argument tilelist file does not exist: {}".format(tilelist_file))
        with open(tilelist_file, 'r') as tilelist_fp:
            supertile_list = [line for line in tilelist_fp.read().splitlines() if line != '']
    else:
        supertile_list = args.tiles.split(',')
    supertile_list = sorted(list(set(supertile_list)))

    res_name = '{}m'.format(args.res)

    root_tiledir = os.path.abspath(args.tiledir)
    scriptdir = SCRIPT_DIR

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,qsub_default)
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)

    if not os.path.isdir(root_tiledir):
        parser.error("tiledir does not exist: {}".format(root_tiledir))

    # Test tiles exist and group into mosaic groups
    mosaic_groups = {}
    for supertile in supertile_list:
        supertile_parts = supertile.split('_')
        if len(supertile_parts) == 2:
            mos = 'polar'
            tnum = supertile
        elif len(supertile_parts) == 3:
            mos = supertile_parts[0]
            tnum = '_'.join(supertile_parts[1:2])
            assert mos.startswith('utm')
        else:
            print("Tile name does not match a known pattern: {}".format(supertile))
            sys.exit(1)

        num_quads_missing_mat = 0
        for q in quadnames:
            tq = "{}_{}".format(supertile, q)
            filename = "{}/{}/{}_{}*.mat".format(root_tiledir, supertile, tq, res_name)
            matfiles = glob.glob(filename)
            if len(matfiles) == 0:
                print("Tile {0} {1} .mat and _reg.mat do not exist: {2}".format(tq, res_name, filename))
                num_quads_missing_mat += 1
            else:
                if not mos in mosaic_groups:
                    mosaic_groups[mos] = []
                mosaic_groups[mos].append(tq)

        dstfps_old_pattern = [
            "{0}/{1}/{1}_*_{2}*.tif".format(root_tiledir, supertile, res_name),
            "{0}/{1}/{1}_*_{2}_meta.txt".format(root_tiledir, supertile, res_name)
        ]
        dstfps_old = [fp for pat in dstfps_old_pattern for fp in glob.glob(pat)]
        if dstfps_old:
            if num_quads_missing_mat == 4:
                print("ERROR! No quad mat files exist, but other MST results exist matching {}".format(dstfps_old_pattern))
                continue
            print("{}Removing old MST results matching {}".format('(dryrun) ' if args.dryrun else '', dstfps_old_pattern))
            if not args.dryrun:
                for dstfp_old in dstfps_old:
                    os.remove(dstfp_old)

    # group tiles by dimension
    groups = {}
    for mos in mosaic_groups:
        existing_tiles = mosaic_groups[mos]
        for quadtile in existing_tiles:
            quadtile_parts = quadtile.split('_')

            if mos == 'polar' and len(quadtile_parts) == 4:
                super_row, super_col, quad_row, quad_col = quadtile_parts
            elif mos.startswith('utm') and len(quadtile) == 5:
                quad_mos, super_row, super_col, quad_row, quad_col = quadtile_parts
                assert quad_mos == mos
            else:
                print("Failed to parse row/col from quad tile name: {}".format(quadtile))
                sys.exit(1)

            row = '_'.join([super_row, quad_row])
            col = '_'.join([super_col, quad_col])

            if args.dimension == 'row':
                key = 'row_{}'.format(row)
            else:
                key = 'col_{}'.format(col)
            key = '{}_{}'.format(mos, key)
    
            if key not in groups:
                groups[key] = []
            groups[key].append(quadtile)

    i=0
    if len(groups) > 0:
        keys = list(groups.keys())
        keys.sort()

        for key in keys:

            print("Submitting tile group from {} {}".format(args.dimension,key))
            quads = groups[key]

            if len(quads) < 2:
                print("Tile group {} has only 1 member: {}. Skipping".format(key, quads))
            else:
                tile_str = ";".join(quads)

                ## if pbs, submit to scheduler
                i+=1
                if args.pbs:
                    job_name = 'tbm2m_{}'.format(key)
                    cmd = r'qsub -N {} -v p1={},p2={},p3="{}",p4={},p5={} {}'.format(
                        job_name,
                        scriptdir,
                        root_tiledir,
                        tile_str,
                        res_name,
                        args.lib_path,
                        qsubpath
                    )
                    print(cmd)
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)

                ## else run matlab
                else:
                    cmd = """matlab -nojvm -nodisplay -nosplash -r "try; addpath('{}'); addpath('{}'); batch_batchMergeTileBuffer('{}',{{'{}'}},'{}'); catch e; disp(getReport(e)); exit(1); end; exit(0)" """.format(
                        scriptdir,
                        args.lib_path,
                        root_tiledir,
                        tile_str.replace(";","','"),
                        res_name
                    )
                    print("{}, {}".format(i, cmd))
                    if not args.dryrun:
                        subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()
