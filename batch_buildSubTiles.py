#!/usr/bin/env python

import argparse
import glob
import os
import subprocess
import warnings
from datetime import datetime

SCRIPT_FILE = os.path.abspath(os.path.realpath(__file__))
SCRIPT_FNAME = os.path.basename(SCRIPT_FILE)
SCRIPT_NAME, SCRIPT_EXT = os.path.splitext(SCRIPT_FNAME)
SCRIPT_DIR = os.path.dirname(SCRIPT_FILE)


## General argument defaults and settings
default_matlab_scriptdir = os.path.join(SCRIPT_DIR, '../setsm_postprocessing4')
default_jobscript = os.path.join(SCRIPT_DIR, 'qsub_buildSubTiles.sh')
default_tempdir = os.path.join(SCRIPT_DIR, 'temp')
default_logdir = '../logs/'  # relative path appended to argument dstdir
swift_program = os.path.abspath('/projects/sciteam/bazu/tools/swift-2/bin/swift')
swift_rundir = '/scratch/sciteam/GS_bazu/user/{}/swiftruns'.format(os.environ['USER'])
swift_config = os.path.join(SCRIPT_DIR, 'swift.conf')
swift_script = os.path.join(SCRIPT_DIR, 'mosaic.swift')
swift_site = 'bst'


## System-specific settings
hostname = os.environ['HOSTNAME'].lower()
if hostname.startswith('h2o'):
    system_name = 'bw'
    sched_presubmit_cmd = 'export NOAPRUNWARN=1'
    sched_addl_envvars = "CRAY_ROOTFS=SHIFTER,UDI='ubuntu:xenial'"
    sched_specify_outerr_paths = True
    sched_addl_vars = "-l nodes=1:ppn=32:xe,gres=shifter,walltime=96:00:00 -m n -q high"
elif hostname.startswith('nunatak'):
    system_name = 'pgc'
    sched_presubmit_cmd = ''
    sched_addl_envvars = ''
    # sched_specify_outerr_paths = True
    # sched_addl_vars = "-l walltime=200:00:00,nodes=1:ppn=16,mem=64gb -m n -q batch"
    sched_specify_outerr_paths = False
    sched_addl_vars = "-l walltime=200:00:00,nodes=1:ppn=16,mem=64gb -m n -k oe -j oe -q batch"
else:
    warnings.warn("Hostname '{}' not recognized. System-specific settings will not be applied.".format(hostname))
    system_name = ''
    sched_presubmit_cmd = ''
    sched_addl_envvars = ''
    sched_specify_outerr_paths = False
    sched_addl_vars = ''


## Argument defaults by 'project'

project_choices = [
    'arcticdem',
    'rema',
    'earthdem',
]

epsg_projstr_dict = {
    None: '',
    3413: 'polar stereo north',
    3031: 'polar stereo south',
}
project_epsg_dict = {
    'arcticdem': 3413,
    'earthdem': None,
    'rema': 3031,
}

earthdem_tileprefix_key = '<tileprefix>'
project_refDemFile_dict = {
    'arcticdem': None,
    'earthdem': '/mnt/pgc/data/elev/dem/tandem-x/90m/TanDEM-X_UTM_90m/TDX_UTM_Mosaic_{}_90m.tif'.format(earthdem_tileprefix_key),
    'rema': '/mnt/pgc/data/elev/dem/tandem-x/90m/TanDEM-X_Antarctica_90m/TanDEM_Antarctica_Mosaic.tif',
}

earthdem_hemisphere_key = '<hemisphere>'
project_tileDefFile_dict = {
    'arcticdem': 'PGC_Imagery_Mosaic_Tiles_Arctic.mat',
    # 'rema': 'PGC_Imagery_Mosaic_Tiles_Antarctic.mat',
    'rema': 'rema_tile_definitions.mat',
    'earthdem': 'PGC_UTM_Mosaic_Tiles_{}.mat'.format(earthdem_hemisphere_key),
}

project_databaseFile_dict = {
    'arcticdem': 'ArcticDEMdatabase4_2m_v4_20201218.mat',
    # 'rema': 'REMAdatabase4_2m_v4_20200806.mat',
    'rema': 'rema_strips_v13e.shp',
    'earthdem': 'EarthDEMdatabase4_2m_v4_20210101_terrnva-paths.mat',
}
project_waterTileDir_dict = {
    'arcticdem': '/mnt/pgc/data/projects/arcticdem/watermasks/global_surface_water/tiled_watermasks/',
    'rema':      '',
    'earthdem':  '/mnt/pgc/data/projects/earthdem/watermasks/global_surface_water/tiled_watermasks/',
}
project_stripsDirectory_dict = {
    'arcticdem': '/mnt/pgc/terrnva_data/elev/dem/setsm/ArcticDEM/region',
    'rema':      '/mnt/pgc/terrnva_data/elev/dem/setsm/REMA/region',
    'earthdem':  '/mnt/pgc/terrnva_data/elev/dem/setsm/EarthDEM/region',
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

    parser = argparse.ArgumentParser()
    parser.add_argument("dstdir", help="target directory (tile subfolders will be created)")
    parser.add_argument("tiles",
        help=' '.join([
            "list of mosaic tiles; either specified on command line (comma delimited),",
            "or a text file list (each tile on separate line)"
        ])
    )

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
    parser.add_argument("--strip-db", default=None,
                        help="strip database mat file (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_databaseFile_dict.items()])
                        ))
    parser.add_argument("--strips-dir", default=None,
                        help="root directory of source strip dems (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_stripsDirectory_dict.items()])
                        ))
    parser.add_argument("--ref-dem", default=None, help="reference DEM (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_refDemFile_dict.items()])
                        ))
    parser.add_argument("--water-tile-dir", default=None,
                        help="directory of water tifs (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_waterTileDir_dict.items()])
                        ))
    parser.add_argument("--tileqc-dir", default=None,
                        help="directory of tile qc mat files (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_tileqcDir_dict.items()])
                        ))
    parser.add_argument("--tileparam-list", default=None,
                        help="tile parameters text file (default is {})".format(
                            ', '.join(["{} if --project={}".format(val, dom) for dom, val in project_tileParamList_dict.items()])
                        ))
    
    parser.add_argument("--libdir", default=default_matlab_scriptdir,
                        help="directory of referenced Matlab functions (default={})".format(default_matlab_scriptdir))

    parser.add_argument("--make-10m-only", action='store_true', default=False,
            help="do not give 'make2m' argument in call to buildSubTiles.m, breaks --rerun capability")
    parser.add_argument("--rerun", action='store_true', default=False,
            help="rerun tile, behavior determined by redoFlag in Matlab code")
    parser.add_argument("--rerun-without-cleanup", action='store_true', default=False,
            help="rerun tile without attempting to clean up potentially corrupted subtiles")
    # parser.add_argument("--sort-fix", action="store_true", default=False,
    #         help="run tile with buildSubTilesSortFix script")

    # parser.add_argument('--require-finfiles', action='store_true', default=False,
    #         help="let existence of finfiles dictate reruns")
    parser.add_argument('--bypass-finfile-req', action='store_true', default=False,
            help="upon rerun, deem tiles complete when the 10,000-th subtile exists even though finfile does not exist")

    parser.add_argument("--jobscript", default=default_jobscript,
            help="jobscript used in task submission (default={})".format(default_jobscript))
    parser.add_argument("--tempdir", default=default_tempdir,
            help="directory where filled-out running copy of jobscript is created (default={})".format(default_tempdir))
    parser.add_argument("--logdir", default=default_logdir,
            help="directory where logfiles for Matlab tile processing are created (default is {} from dstdir)".format(default_logdir))

    parser.add_argument("--pbs", action='store_true', default=False,
            help="submit tasks to PBS")
    parser.add_argument("--slurm", action='store_true', default=False,
            help="submit tasks to SLURM")
    parser.add_argument("--swift", action='store_true', default=False,
            help="submit tasks with Swift")
    parser.add_argument("--swift-program", default=swift_program,
            help="path to Swift executable (default={})".format(swift_program))
    parser.add_argument("--swift-logdir", default=None,
            help="directory where output job logs from Swift will go (default is ../logs/ from dstdir)")

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

    make2m_arg = 'false' if args.make_10m_only else 'true'
    target_res = '10m' if args.make_10m_only else '2m'

    ## Set default arguments by project setting
    if args.project is None and True in [arg is None for arg in
            [args.epsg, args.tile_def, args.strip_db, args.ref_dem, args.water_tile_dir]]:
        parser.error("--project arg must be provided if one of the following arguments is not provided: {}".format(
            ' '.join(["--epsg", "--tile-def", "--strip-db", "--ref-dem", "--water-tile-dir"])
        ))
    if args.epsg is None:
        args.epsg = project_epsg_dict[args.project]
    projection_string = epsg_projstr_dict[args.epsg]
    if args.tile_def is None:
        args.tile_def = project_tileDefFile_dict[args.project]
    if args.strip_db is None:
        args.strip_db = project_databaseFile_dict[args.project]
    if args.strips_dir is None:
        args.strips_dir = project_stripsDirectory_dict[args.project]
    if args.ref_dem is None:
        args.ref_dem = project_refDemFile_dict[args.project]
        if args.ref_dem is None:
            parser.error("--ref-dem argument must be provided if --project={}".format(args.project))
    if args.water_tile_dir is None:
        args.water_tile_dir = project_waterTileDir_dict[args.project]
    if args.tileqc_dir is None:
        args.tileqc_dir = project_tileqcDir_dict[args.project]
    if args.tileparam_list is None:
        args.tileparam_list = project_tileParamList_dict[args.project]

    ## Convert argument paths to absolute paths
    args.dstdir = os.path.abspath(args.dstdir)
    args.tile_def = os.path.abspath(args.tile_def if os.path.isfile(args.tile_def) else os.path.join(SCRIPT_DIR, args.tile_def))
    args.strip_db = os.path.abspath(args.strip_db if os.path.isfile(args.strip_db) else os.path.join(SCRIPT_DIR, args.strip_db))
    if args.strips_dir != '':
        args.strips_dir = os.path.abspath(args.strips_dir)
    if args.ref_dem != '':
        args.ref_dem = os.path.abspath(args.ref_dem)
    if args.water_tile_dir != '':
        args.water_tile_dir = os.path.abspath(args.water_tile_dir)
    if args.tileqc_dir != '':
        args.tileqc_dir = os.path.abspath(args.tileqc_dir)
    if args.tileparam_list != '':
        args.tileparam_list = os.path.abspath(args.tileparam_list)
    args.libdir = os.path.abspath(args.libdir)
    args.tempdir = os.path.abspath(args.tempdir)
    args.logdir = os.path.abspath(os.path.join(args.dstdir, args.logdir) if args.logdir.startswith('../') else args.logdir)

    ## Verify path arguments
    if not os.path.isdir(args.dstdir):
        parser.error("dstdir does not exist: {}".format(args.dstdir))
    if args.project != 'earthdem':
        if not os.path.isfile(args.tile_def):
            parser.error("--tile-def file does not exit: {}".format(args.tile_def))
        if args.ref_dem is not None and not os.path.isfile(args.ref_dem):
            parser.error("--ref-dem does not exist: {}".format(args.ref_dem))
    if not os.path.isfile(args.strip_db):
        parser.error("--strip-db does not exist: {}".format(args.strip_db))
    if args.strips_dir != '' and not os.path.isdir(args.strips_dir):
        parser.error("--strips-dir does not exist: {}".format(args.strips_dir))
    if args.water_tile_dir != '' and not os.path.isdir(args.water_tile_dir):
        parser.error("--water-tile-dir does not exist: {}".format(args.water_tile_dir))
    if args.tileqc_dir != '' and not os.path.isdir(args.tileqc_dir):
        parser.error("--tileqc-dir does not exist: {}".format(args.tileqc_dir))
    if args.tileparam_list != '' and not os.path.isfile(args.tileparam_list):
        parser.error("--tileqc-file does not exist: {}".format(args.tileparam_list))
    if not os.path.isdir(args.libdir):
        parser.error("--libdir does not exist: {}".format(args.libdir))
    if not os.path.isfile(args.jobscript):
        parser.error("--jobscript does not exist: {}".format(args.jobscript))

    ## Verify other arguments
    if [args.pbs, args.slurm, args.swift].count(True) > 1:
        parser.error("--pbs --slurm --swift are mutually exclusive")
    if args.rerun and args.rerun_without_cleanup:
        parser.error("--rerun and --rerun-without-cleanup are mutually exclusive")

    bst_logdir = os.path.join(args.logdir, target_res, 'bst')
    pbs_logdir = os.path.join(args.logdir, 'pbs', target_res, 'bst')
    swift_logrootdir = os.path.join(args.logdir, 'swift')
    swift_logdir = os.path.join(swift_logrootdir, target_res, 'bst')
    swift_tasklist_dir = os.path.join(swift_logrootdir, 'tasklist')
    swift_tasklist_fname = 'bst_{}_tasklist_{}.txt'.format(target_res, datetime.now().strftime("%Y%m%d%H%M%S"))
    swift_tasklist_file = os.path.join(swift_tasklist_dir, swift_tasklist_fname)

    ## Create output directories
    outdir_list = [args.tempdir, bst_logdir]
    if args.pbs:
        outdir_list.append(pbs_logdir)
    if args.swift:
        outdir_list.extend([swift_rundir, swift_logdir, swift_tasklist_dir])
    for outdir in outdir_list:
        if not os.path.isdir(outdir):
            print("Creating output directory: {}".format(outdir))
            os.makedirs(outdir)

    ## Create temp jobscript with comment mosaicking args filled in
    tilename_key = '<tilename>'
    template_outdir = os.path.join(args.dstdir, tilename_key, 'subtiles')
    template_finfile = "{}_{}.fin".format(template_outdir, target_res)
    template_logfile = os.path.join(bst_logdir, tilename_key+'.log')
    jobscript_static_args_dict = {
        'system': system_name,
        'scriptdir': SCRIPT_DIR,
        'libdir': args.libdir,
        'outDir': template_outdir,
        'projection': projection_string,
        'tileDefFile': args.tile_def,
        'stripDatabaseFile': args.strip_db,
        'stripsDirectory': args.strips_dir,
        'waterTileDir': args.water_tile_dir,
        'refDemFile': args.ref_dem,
        'tileqcDir': args.tileqc_dir,
        'tileParamListFile': args.tileparam_list,
        'make2m': make2m_arg,
        'finfile': template_finfile,
        'logfile': template_logfile,
    }
    jobscript_fname = os.path.basename(args.jobscript)
    jobscript_temp_fname = jobscript_fname.replace(
        '.sh',
        '{}_{}.sh'.format(
            '_{}'.format(args.project) if args.project is not None else '',
            datetime.now().strftime("%Y%m%d%H%M%S")
        )
    )

    jobscript_temp = os.path.join(args.tempdir, jobscript_temp_fname)
    with open(args.jobscript, 'r') as jobscript_fp:
        jobscript_text = jobscript_fp.read()

    jobscript_temp_text = jobscript_text
    for argname, argval in jobscript_static_args_dict.items():
        jobscript_argname = '"$ARG_{}"'.format(argname).upper()
        jobscript_argval = '"{}"'.format(argval)
        jobscript_temp_text = jobscript_temp_text.replace(jobscript_argname, jobscript_argval)

    print("Writing temporary jobscript file: {}".format(jobscript_temp))
    with open(jobscript_temp, 'w') as jobscript_fp:
        jobscript_fp.write(jobscript_temp_text)


    tiles_to_run = []

    for tile in tiles:

        tile_projstr = projection_string
        tile_def = args.tile_def
        ref_dem = args.ref_dem

        if tile_projstr == '' or earthdem_hemisphere_key in tile_def or earthdem_tileprefix_key in ref_dem:
            assert args.project == 'earthdem'

            utm_tilename_parts = tile.split('_')
            utm_tilename_prefix = utm_tilename_parts[0]
            if not utm_tilename_prefix.startswith('utm'):
                parser.error("Expected only UTM tile names (e.g. 'utm10n_01_01'), but got '{}'".format(tile))

            if tile_projstr == '':
                tile_projstr = utm_tilename_prefix

            if earthdem_hemisphere_key in tile_def:
                if utm_tilename_prefix.endswith('n'):
                    hemisphere = 'North'
                elif utm_tilename_prefix.endswith('s'):
                    hemisphere = 'South'
                else:
                    parser.error("UTM tile name prefix does not end with 'n' or 's' (e.g. 'utm10n'): {}".format(tile))

                tile_def = tile_def.replace(earthdem_hemisphere_key, hemisphere)
                if not os.path.isfile(tile_def):
                    parser.error("Tile definition file does not exist: {}".format(tile_def))

            if earthdem_tileprefix_key in ref_dem:
                ref_dem = ref_dem.replace(earthdem_tileprefix_key, utm_tilename_prefix)
                if not os.path.isfile(ref_dem):
                    parser.error("Reference DEM file does not exist: {}".format(tile_def))

        ## If output does not exist, add to task list
        tile_outdir = os.path.join(args.dstdir, tile, 'subtiles')
        tile_outdir_exists = os.path.isdir(tile_outdir)
        if not os.path.isdir(tile_outdir):
            if not args.dryrun:
                os.makedirs(tile_outdir)

        final_subtile_fp_10m = os.path.join(tile_outdir, '{}_10000_10m.mat'.format(tile))
        final_subtile_fp_2m = os.path.join(tile_outdir, '{}_10000_2m.mat'.format(tile))
        finfile_10m = "{}_10m.fin".format(tile_outdir)
        finfile_2m = "{}_2m.fin".format(tile_outdir)
        if args.make_10m_only:
            final_subtile_fp = final_subtile_fp_10m
            finfile = finfile_10m
        else:
            final_subtile_fp = final_subtile_fp_2m
            finfile = finfile_2m

        run_tile = True

        if tile_outdir_exists:
            # if args.sort_fix or args.rerun_without_cleanup:
            if args.rerun_without_cleanup:
                pass

            elif args.rerun:
                print("Verifying tile {} before rerun".format(tile))

                # if os.path.isfile(final_subtile_fp) and not args.require_finfiles:
                if os.path.isfile(final_subtile_fp) and args.bypass_finfile_req:
                    print("Tile seems complete ({} exists)".format(os.path.basename(final_subtile_fp)))
                    run_tile = False
                elif os.path.isfile(finfile):
                    print("Tile seems complete ({} exists)".format(os.path.basename(finfile)))
                    run_tile = False
                elif os.path.isfile(finfile_2m):
                    print("Tile seems complete (2m finfile {} exists)".format(os.path.basename(finfile_2m)))
                    run_tile = False

                # if not args.make_10m_only:
                #     ## Remove subtiles with only 10m version
                #     tile_outfiles_10m = glob.glob(os.path.join(tile_outdir, '{}_*10m.mat'.format(tile)))
                #     for outfile_10m in tile_outfiles_10m:
                #         outfile_2m = outfile_10m.replace('10m.mat', '2m.mat')
                #         if not os.path.isfile(outfile_2m):
                #             print("Removing 10m subtile missing 2m component: {}".format(os.path.basename(outfile_10m)))
                #             run_tile = True
                #             if not args.dryrun:
                #                 os.remove(outfile_10m)

            elif any(os.scandir(tile_outdir)) > 0:
                print("Subtiles exist, skipping tile {}".format(tile))
                run_tile = False

        if run_tile:
            tiles_to_run.append(tile)


    print("Running {} tiles".format(len(tiles_to_run)))

    if args.swift:
        with open(swift_tasklist_file, 'w') as swift_tasklist_fp:
            for tile in tiles_to_run:
                swift_tasklist_fp.write(tile+'\n')
        swift_cmd = r""" "{}" -config "{}" -sites {},local "{}" -tasklist_file="{}" -jobscript="{}" -task_description="{}" -logdir="{}" """.format(
            swift_program,
            swift_config,
            swift_site,
            swift_script,
            swift_tasklist_file,
            jobscript_temp,
            "tiles to run",
            swift_logdir,
        )
        print("Running Swift command: {}".format(swift_cmd))
        if not args.dryrun:
            subprocess.call(swift_cmd, shell=True, cwd=swift_rundir)
        print("Ran {} tiles".format(len(tiles_to_run)))

    else:
        for tasknum, tile in enumerate(tiles_to_run, 1):

            job_name = 'bst_{}'.format(tile)
            job_outfile = os.path.join(pbs_logdir, tile+'.out')
            job_errfile = os.path.join(pbs_logdir, tile+'.err')

            if args.pbs:
                cmd = r""" {}qsub -N {} -v ARG_TILENAME={}{} {} {} "{}" """.format(
                    sched_presubmit_cmd+' ; ' if sched_presubmit_cmd != '' else '',
                    job_name,
                    tile,
                    ','+sched_addl_envvars if sched_addl_envvars != '' else '',
                    '-o "{}" -e "{}"'.format(job_outfile, job_errfile) if sched_specify_outerr_paths else '',
                    sched_addl_vars,
                    jobscript_temp,
                )

            elif args.slurm:
                cmd = r""" {}sbatch -J {} -v ARG_TILENAME={}{} {} {} "{}" """.format(
                    sched_presubmit_cmd+' ; ' if sched_presubmit_cmd != '' else '',
                    job_name,
                    tile,
                    ','+sched_addl_envvars if sched_addl_envvars != '' else '',
                    '-o "{}" -e "{}"'.format(job_outfile, job_errfile) if sched_specify_outerr_paths else '',
                    sched_addl_vars,
                    jobscript_temp,
                )

            else:
                cmd = r""" bash "{}" {} """.format(
                    jobscript_temp,
                    tile,
                )

            print("{}, {}".format(tasknum, cmd))
            if not args.dryrun:
                subprocess.call(cmd, shell=True, cwd=(pbs_logdir if args.pbs else None))

        if args.pbs or args.slurm:
            print("Submitted {} tiles to scheduler".format(len(tiles_to_run)))
        else:
            print("Ran {} tiles".format(len(tiles_to_run)))


    print("Done")


if __name__ == '__main__':
    main()
