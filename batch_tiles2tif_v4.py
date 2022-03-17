import os, string, sys, argparse, glob, subprocess

SCRIPT_FILE = os.path.abspath(os.path.realpath(__file__))
SCRIPT_FNAME = os.path.basename(SCRIPT_FILE)
SCRIPT_NAME, SCRIPT_EXT = os.path.splitext(SCRIPT_FNAME)
SCRIPT_DIR = os.path.dirname(SCRIPT_FILE)

matlab_scripts = os.path.join(SCRIPT_DIR, '../setsm_postprocessing4')
quadnames = ('1_1','1_2','2_1','2_2')
RESOLUTIONS = ['2','10']
REGIONS = ['arctic','antarctic','earthdem']
default_qsub = 'qsub_tiles2tif_v4.sh'

#### TODO add projstring to passed args
def main():
    global quadnames

    ## args
    parser = argparse.ArgumentParser()
    parser.add_argument("dstdir", help="target directory (tile subfolders will be created)")
    parser.add_argument("tiles",
        help=' '.join([
            "list of mosaic tiles; either specified on command line (comma delimited),",
            "or a text file list (each tile on separate line)"
        ])
    )
    parser.add_argument("res", choices=RESOLUTIONS, help="resolution ({})".format(','.join(RESOLUTIONS)))
    parser.add_argument("region", choices=REGIONS, help="region ({})".format(','.join(REGIONS)))

    parser.add_argument("--meta-only", action='store_true', default=False,
                        help="build meta files only")
    parser.add_argument("--rerun", action='store_true', default=False,
            help="run script even if target dem already exists")
    parser.add_argument("--lib-path", default=matlab_scripts,
            help="path to referenced Matlab functions (default={}".format(matlab_scripts))
    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--hold", action='store_true', default=False,
            help="when submitting PBS jobs, submit them with held (H) status")
    parser.add_argument("--qsubscript",
            help="qsub script to use in PBS submission (default is {} in script root folder)".format(default_qsub))
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')

    args = parser.parse_args()

    if args.tiles.lower().endswith(('.txt', '.csv')) or os.path.isfile(args.tiles):
        tilelist_file = args.tiles
        if not os.path.isfile(args.tiles):
            parser.error("'tiles' argument tilelist file does not exist: {}".format(tilelist_file))
        with open(tilelist_file, 'r') as tilelist_fp:
            tiles = [line for line in tilelist_fp.read().splitlines() if line != '']
    else:
        tiles = args.tiles.split(',')
    tiles = sorted(list(set(tiles)))

    dstdir = os.path.abspath(args.dstdir)
    scriptdir = SCRIPT_DIR

    ## Verify qsubscript
    if args.qsubscript is None:
        qsubpath = os.path.join(scriptdir,default_qsub)
    else:
        qsubpath = os.path.abspath(args.qsubscript)
    if not os.path.isfile(qsubpath):
        parser.error("qsub script path is not valid: %s" %qsubpath)

    if not os.path.isdir(dstdir):
        parser.error("dstdir does not exist: {}".format(dstdir))

    if args.region == 'arctic':
        projstr = 'polar stereo north'
    elif args.region == 'antarctic':
        projstr = 'polar stereo south'
    elif args.region == 'earthdem':
        projstr = None

    if args.res == '10':
        quadnames = ['']

    i=0
    if len(tiles) > 0:

        for tile in tiles:

            tile_projstr = projstr

            if tile_projstr is None:
                assert args.region == 'earthdem'

                utm_tilename_parts = tile.split('_')
                utm_tilename_prefix = utm_tilename_parts[0]
                if not utm_tilename_prefix.startswith('utm'):
                    parser.error("Expected only UTM tile names (e.g. 'utm10n_01_01'), but got '{}'".format(tile))

                tile_projstr = utm_tilename_prefix

            for q in quadnames:
                tq = "{}_{}".format(tile,q) if args.res == '2' else tile
                dstfp = os.path.join(dstdir,tile,'{}_{}m_dem.tif'.format(tq, args.res))
                metafp = os.path.join(dstdir,tile,'{}_{}m_meta.txt'.format(tq, args.res))
                unregmatfile = os.path.join(dstdir,tile,'{}_{}m.mat'.format(tq, args.res))
                regmatfile = os.path.join(dstdir,tile,'{}_{}m_reg.mat'.format(tq, args.res))
                matfile = regmatfile if os.path.isfile(regmatfile) else unregmatfile

                run_tile = True

                if not os.path.isfile(unregmatfile) and not os.path.isfile(regmatfile):
                    print("Tile {} {}m mat and reg.mat files do not exist: {}".format(tq,args.res,matfile))
                    run_tile = False

                elif not args.meta_only and os.path.isfile(dstfp):
                    if args.rerun or not os.path.isfile(metafp):
                        dstfps_old_pattern = [
                            dstfp.replace('_dem.tif', '*.tif'),
                            metafp
                        ]
                        dstfps_old = [fp for pat in dstfps_old_pattern for fp in glob.glob(pat)]
                        if dstfps_old:
                            print("{}Removing existing tif tile results matching {}".format('(dryrun) ' if args.dryrun else '', dstfps_old_pattern))
                            if not args.dryrun:
                                for dstfp_old in dstfps_old:
                                    os.remove(dstfp_old)
                    else:
                        print('{} exists, skipping'.format(dstfp))
                        run_tile = False
                
                elif args.meta_only and os.path.isfile(metafp):
                    if args.rerun:
                        print('Removing existing meta file: {}'.format(metafp))
                        if not args.dryrun:
                            os.remove(metafp)
                    else:
                        print('{} exists, skipping'.format(metafp))
                        run_tile = False

                if run_tile:
                    ## if pbs, submit to scheduler
                    i+=1
                    if args.pbs:
                        job_name = 't2t_{}'.format(tq)
                        cmd = r'qsub {} -N {} -v p1={},p2={},p3="{}",p4={},p5={} {}'.format(
                            '-h' if args.hold else '',
                            job_name,
                            scriptdir,
                            matfile,
                            tile_projstr,
                            args.lib_path,
                            'true' if args.meta_only else 'false',
                            qsubpath,
                        )
                        print(cmd)
                        if not args.dryrun:
                            subprocess.call(cmd, shell=True)

                    ## else run matlab
                    else:
                        if args.meta_only:
                            cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); tileMetav4('{2}'); exit" """.format(
                                scriptdir,
                                args.lib_path,
                                matfile,
                            )
                        else:
                            cmd = """matlab -nojvm -nodisplay -nosplash -r "addpath('{0}'); addpath('{1}'); writeTileToTifv4('{2}','{3}','outRasterType','COG'); tileMetav4('{2}'); exit" """.format(
                                scriptdir,
                                args.lib_path,
                                matfile,
                                tile_projstr,
                            )
                        print("{}, {}".format(i, cmd))
                        if not args.dryrun:
                            subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()
