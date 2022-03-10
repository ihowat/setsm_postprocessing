#!/usr/bin/env python

import argparse
import glob
import math
import os
import socket
import subprocess
import sys
import time
import warnings
from datetime import datetime

SCRIPT_FILE = os.path.abspath(os.path.realpath(__file__))
SCRIPT_FNAME = os.path.basename(SCRIPT_FILE)
SCRIPT_NAME, SCRIPT_EXT = os.path.splitext(SCRIPT_FNAME)
SCRIPT_DIR = os.path.dirname(SCRIPT_FILE)


## General argument defaults and settings
mst_pyscript = os.path.join(SCRIPT_DIR, 'batch_mosaicSubTiles.py')
default_matlab_scriptdir = os.path.join(SCRIPT_DIR, '../setsm_postprocessing4')
default_bst_jobscript = os.path.join(SCRIPT_DIR, 'qsub_buildSubTiles.sh')
chain_jobscript = os.path.join(SCRIPT_DIR, 'qsub_runBstThenMst.sh')
batch_jobscript = os.path.join(SCRIPT_DIR, 'qsub_batchRunTiles.sh')
default_tempdir = os.path.join(SCRIPT_DIR, 'temp')
default_logdir = '../logs/[DSTDIR_FOLDER_NAME]'  # relative path appended to argument dstdir
swift_program = os.path.abspath('/projects/sciteam/bazu/tools/swift-2/bin/swift')
swift_rundir = '/scratch/sciteam/GS_bazu/user/{}/swiftruns'.format(os.environ['USER'])
swift_config = os.path.join(SCRIPT_DIR, 'swift.conf')
swift_script = os.path.join(SCRIPT_DIR, 'mosaic.swift')
swift_site = 'bst'


## System-specific settings
hostname = socket.gethostname().lower()
if hostname.startswith('h2o'):
    system_name = 'bw'
    sched_presubmit_cmd = 'export NOAPRUNWARN=1'
    sched_addl_envvars = "CRAY_ROOTFS=SHIFTER,UDI='ubuntu:xenial'"
    sched_specify_outerr_paths = True
    sched_addl_vars = "-l nodes=1:ppn=16:xe,gres=shifter,walltime=96:00:00 -m n"
    sched_default_queue = 'normal'
elif hostname.startswith('nunatak'):
    system_name = 'pgc'
    sched_presubmit_cmd = ''
    sched_addl_envvars = ''
    sched_specify_outerr_paths = False
    sched_addl_vars = "-l nodes=1:ppn=16,mem=64gb,walltime=200:00:00 -m n -k oe -j oe"
    sched_default_queue = 'old'
else:
    warnings.warn("Hostname '{}' not recognized. System-specific settings will not be applied.".format(hostname))
    system_name = ''
    sched_presubmit_cmd = ''
    sched_addl_envvars = ''
    sched_specify_outerr_paths = False
    sched_addl_vars = ''
    sched_default_queue = None


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
    'earthdem':  None,
    'rema':      3031,
}

earthdem_tilePrefix_key = '<tilePrefix>'
project_refDemFile_dict = {
    'arcticdem': None,
    'earthdem':  '/mnt/pgc/data/elev/dem/tandem-x/90m/TanDEM-X_UTM_90m/TDX_UTM_Mosaic_{}_90m.tif'.format(earthdem_tilePrefix_key),
    'rema':      '/mnt/pgc/data/elev/dem/tandem-x/90m/TanDEM-X_Antarctica_90m/TanDEM_Antarctica_Mosaic.tif',
}

earthdem_hemisphere_key = '<hemisphere>'
project_tileDefFile_dict = {
    'arcticdem': '/mnt/pgc/data/projects/earthdem/tiledef_files/PGC_Imagery_Mosaic_Tiles_Arctic.mat',
    'rema':      '/mnt/pgc/data/projects/earthdem/tiledef_files/rema_tile_definitions.mat',
    'earthdem':  '/mnt/pgc/data/projects/earthdem/tiledef_files/PGC_UTM_Mosaic_Tiles_{}.mat'.format(earthdem_hemisphere_key),
}

project_databaseFile_dict = {
    'arcticdem': '/scratch/sciteam/GS_bazu/mosaic_data/strip_databases/ArcticDEMdatabase4_2m_v4_20210817.mat',
    'rema':      '/mnt/pgc/data/projects/earthdem/strip_databases/rema_strips_v13e.shp',
    'earthdem':  '/scratch/sciteam/GS_bazu/mosaic_data/strip_databases/EarthDEMdatabase4_2m_v4_20211014.mat',
}
project_waterTileDir_dict = {
    'arcticdem': '/mnt/pgc/data/projects/arcticdem/watermasks/',
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
    'rema':      '/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/tile_qc/tile_qc_v13e',
    'earthdem':  '',
}
project_tileParamList_dict = {
    'arcticdem': '',
    'rema':      '/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/tile_params/tileParamList_v13e.txt',
    'earthdem':  '',
}


with open(os.path.join(SCRIPT_DIR, 'watermask_tiles_greenland.txt'), 'r') as tilelist_fp:
    watermask_tiles_greenland = [line.strip() for line in tilelist_fp.read().splitlines() if line != '']
with open(os.path.join(SCRIPT_DIR, 'watermask_tiles_needing_visnav.txt'), 'r') as tilelist_fp:
    watermask_tiles_needing_visnav = [line.strip() for line in tilelist_fp.read().splitlines() if line != '']
watermask_tiles_visnav_need_editing_file = os.path.join(SCRIPT_DIR, 'watermask_tiles_visnav_need_editing_canada.txt')
with open(watermask_tiles_visnav_need_editing_file, 'r') as tilelist_fp:
    watermask_tiles_visnav_need_editing = [line.strip() for line in tilelist_fp.read().splitlines() if line != '']

quads = ['1_1', '1_2', '2_1', '2_2']


def main():
    global sched_addl_vars

    parser = argparse.ArgumentParser()
    parser.add_argument("dstdir", help="target directory (where tile subfolders will be created)")
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
    parser.add_argument("--chain-mst", action='store_true', default=False,
            help=("Run 10m and 2m mosaicSubTiles processes for each tile after buildSubTiles completes successfully. "
                  "Delete 'subtiles' folder for the tile upon successful completion of 10m and 2m MST processes. "
                  "Only 10m MST process will be run if --make-10m-only option is provided."))
    parser.add_argument("--chain-mst-keep-subtiles", action='store_true', default=False,
            help="do not remove tiles' 'subtiles' folders after successful BST+MST steps")
    parser.add_argument("--chain-mst-no-local", action='store_true', default=False,
            help="do not use local space on compute nodes for tiles' 'subtiles' folders")

    # parser.add_argument('--require-finfiles', action='store_true', default=False,
    #         help="let existence of finfiles dictate reruns")
    parser.add_argument('--bypass-finfile-req', action='store_true', default=False,
            help="upon rerun, deem tiles complete when the 10,000-th subtile exists even though finfile does not exist")
    parser.add_argument('--skip-missing-watermasks', action='store_true', default=False,
            help="if tiles are missing watermasks, skip them and submit the rest of the tiles for processing")

    parser.add_argument("--jobscript", default=default_bst_jobscript,
            help="jobscript used in BST task submission (default={})".format(default_bst_jobscript))
    parser.add_argument("--tempdir", default=default_tempdir,
            help="directory where filled-out running copy of jobscript is created (default={})".format(default_tempdir))
    parser.add_argument("--logdir", default=default_logdir,
            help="directory where logfiles for Matlab tile processing are created (default is {} from dstdir)".format(default_logdir))

    parser.add_argument("--queue", default=sched_default_queue,
            help="queue for scheduler submission")
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

    parser.add_argument("--tasks-per-job", type=int, default=1,
            help="Number of tiles to run in a single PBS or Slurm job")
    parser.add_argument("--dryrun", action='store_true', default=False,
            help='print actions without executing')

    args = parser.parse_args()

    if args.chain_mst_no_local or args.chain_mst_keep_subtiles:
        args.chain_mst = True

    if args.tiles.lower().endswith(('.txt', '.csv')) or os.path.isfile(args.tiles):
        tilelist_file = args.tiles
        if not os.path.isfile(args.tiles):
            parser.error("'tiles' argument tilelist file does not exist: {}".format(tilelist_file))
        with open(tilelist_file, 'r') as tilelist_fp:
            tiles = [line for line in tilelist_fp.read().splitlines() if line != '']
    else:
        tiles = args.tiles.split(',')
    tiles = sorted(list(set(tiles)))

    make_10m_only = 'true' if args.make_10m_only else 'false'
    make2m_arg = 'false' if args.make_10m_only else 'true'
    target_res = '10m' if args.make_10m_only else '2m'

    if args.tasks_per_job > 1:
        if not (args.pbs or args.slurm):
            parser.error("--tasks-per-job > 1 can only be provided in addition to"
                         " --pbs or --slurm scheduler job submission options")
        batch_job_submission = True
        if args.pbs:
            if 'nodes=1:' not in sched_addl_vars:
                parser.error("'nodes=1' must be in PBS resource request for --tasks-per-job > 1")
            sched_addl_vars = sched_addl_vars.replace('nodes=1', 'nodes={}'.format(args.tasks_per_job))
        if args.slurm:
            # TODO: Implement --tasks-per-job for SLURM job submission
            pass
    else:
        batch_job_submission = False

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
    auto_select_arcticdem_water_tile_dir = False
    if args.water_tile_dir is None:
        args.water_tile_dir = project_waterTileDir_dict[args.project]
        if args.project == 'arcticdem':
            auto_select_arcticdem_water_tile_dir = True
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
    if '[DSTDIR_FOLDER_NAME]' in args.logdir:
        args.logdir = args.logdir.replace('[DSTDIR_FOLDER_NAME]', os.path.basename(args.dstdir))

    ## Verify path arguments
    if not os.path.isdir(args.dstdir):
        parser.error("dstdir does not exist: {}".format(args.dstdir))
    if args.project != 'earthdem':
        if not os.path.isfile(args.tile_def):
            parser.error("--tile-def file does not exist: {}".format(args.tile_def))
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
    if args.chain_mst and args.project is None:
        parser.error('--chain-mst option requires --project arg is also provided')

    log_time = datetime.now().strftime("%Y%m%d%H%M%S")
    bst_logdir = os.path.join(args.logdir, 'matlab', 'bst', target_res)
    pbs_logdir = os.path.join(args.logdir, 'pbs', 'bst', target_res)
    if batch_job_submission:
        pbs_logdir = '{}_batch_{}_{}-tiles'.format(pbs_logdir, log_time, len(tiles))
    swift_logrootdir = os.path.join(args.logdir, 'swift')
    swift_logdir = os.path.join(swift_logrootdir, 'bst', target_res)
    swift_tasklist_dir = os.path.join(swift_logrootdir, 'tasklist')
    swift_tasklist_fname = 'bst_{}_tasklist_{}.txt'.format(target_res, datetime.now().strftime("%Y%m%d%H%M%S"))
    swift_tasklist_file = os.path.join(swift_tasklist_dir, swift_tasklist_fname)

    ## Create output directories
    outdir_list = [args.tempdir, bst_logdir]
    if args.pbs:
        outdir_list.append(pbs_logdir)
    if args.swift:
        outdir_list.extend([swift_rundir, swift_logdir, swift_tasklist_dir])
    if not args.dryrun:
        for outdir in outdir_list:
            if not os.path.isdir(outdir):
                print("Creating output directory: {}".format(outdir))
                os.makedirs(outdir)
        if batch_job_submission:
            pbs_batch_tilelist = os.path.join(pbs_logdir, 'tilelist.txt')
            with open(pbs_batch_tilelist, 'w') as tilelist_fp:
                for tile in tiles:
                    tilelist_fp.write(tile+'\n')


    ## Create temp jobscript with mosaicking args filled in
    tilename_key = '<tileName>'
    template_outdir = os.path.join(args.dstdir, tilename_key, 'subtiles')
    template_finfile = "{}_{}.fin".format(template_outdir, target_res)
    template_logfile = os.path.join(bst_logdir, tilename_key+'.log')
    template_runscript = "{}_startmatlab_{}.sh".format(template_outdir, target_res)
    jobscript_static_args_dict = {
        'system': system_name,
        'scriptdir': SCRIPT_DIR,
        'libdir': args.libdir,
        'outDir': template_outdir,
        'projection': projection_string,
        'tileDefFile': args.tile_def,
        'stripDatabaseFile': args.strip_db,
        'stripsDirectory': args.strips_dir,
        'refDemFile': args.ref_dem,
        'tileqcDir': args.tileqc_dir,
        'tileParamListFile': args.tileparam_list,
        'make2m': make2m_arg,
        'finfile': template_finfile,
        'logfile': template_logfile,
        'runscript': template_runscript,
    }
    if not auto_select_arcticdem_water_tile_dir:
        jobscript_static_args_dict['waterTileDir'] = args.water_tile_dir
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

    print("Writing temporary BST jobscript file: {}".format(jobscript_temp))
    with open(jobscript_temp, 'w') as jobscript_fp:
        jobscript_fp.write(jobscript_temp_text)
    bst_jobscript = jobscript_temp


    if args.chain_mst:

        ## Create temp MST jobscript
        mst_jobscript = None
        mst_args = ['python', mst_pyscript, args.dstdir, tiles[0], '10', '--project', args.project, '--write-jobscript-and-exit']
        if not args.make_10m_only:
            mst_args.append('--make-2m-logdirs')
        print("Creating MST temporary jobscript file with the following command:\n    {}".format(' '.join(mst_args)))
        mst_proc = subprocess.Popen(mst_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        mst_stdout, mst_stderr = mst_proc.communicate()
        mst_stdout = mst_stdout.decode('utf-8').strip()
        mst_stderr = mst_stderr.decode('utf-8').strip()
        if mst_stdout != '':
            print("STDOUT:\n'''\n{}\n'''".format(mst_stdout))
        if mst_stderr != '':
            print("STDERR:\n'''\n{}\n'''".format(mst_stderr))
        mst_out_lines = mst_stdout.splitlines()
        for line in mst_out_lines:
            if "Writing temporary jobscript file:" in line:
                mst_jobscript = line.split(':')[1].strip()
        if mst_jobscript is None:
            parser.error("Failed to parse path to MST temporary jobscript from stdout")
            
        ## Create temp chain jobscript
        jobscript_static_args_dict = {
            'system': system_name,
            'bst_jobscript': bst_jobscript,
            'mst_jobscript': mst_jobscript,
            'mst_pyscript': mst_pyscript,
            'output_tiles_dir': args.dstdir,
            'project': args.project,
            'make_10m_only': make_10m_only,
            'keep_subtiles': str(args.chain_mst_keep_subtiles).lower(),
            'use_local': str(not args.chain_mst_no_local).lower(),
        }
        if not auto_select_arcticdem_water_tile_dir:
            jobscript_static_args_dict['waterTileDir'] = args.water_tile_dir
        jobscript_fname = os.path.basename(chain_jobscript)
        jobscript_temp_fname = jobscript_fname.replace(
            '.sh',
            '{}_{}.sh'.format(
                '_{}'.format(args.project) if args.project is not None else '',
                datetime.now().strftime("%Y%m%d%H%M%S")
            )
        )

        jobscript_temp = os.path.join(args.tempdir, jobscript_temp_fname)
        with open(chain_jobscript, 'r') as jobscript_fp:
            jobscript_text = jobscript_fp.read()

        jobscript_temp_text = jobscript_text
        for argname, argval in jobscript_static_args_dict.items():
            jobscript_argname = '"$ARG_{}"'.format(argname).upper()
            jobscript_argval = '"{}"'.format(argval)
            jobscript_temp_text = jobscript_temp_text.replace(jobscript_argname, jobscript_argval)

        print("Writing temporary chain jobscript file: {}".format(jobscript_temp))
        with open(jobscript_temp, 'w') as jobscript_fp:
            jobscript_fp.write(jobscript_temp_text)

    tilerun_jobscript = jobscript_temp


    number_incomplete_tiles = 0
    tiles_to_run = []
    tiles_missing_watermask = []
    tile_waterTileDir_dict = dict()
    waterTileDir_tilesToRun_dict = dict()
    incomplete_tiles_missing_watermask = []

    for tile in tiles:

        tile_projstr = projection_string
        tile_def = args.tile_def
        ref_dem = args.ref_dem
        water_tile_dir = args.water_tile_dir
        tile_missing_watermask = False

        if tile_projstr == '' or earthdem_hemisphere_key in tile_def or earthdem_tilePrefix_key in ref_dem:
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

            if earthdem_tilePrefix_key in ref_dem:
                ref_dem = ref_dem.replace(earthdem_tilePrefix_key, utm_tilename_prefix)
                if not os.path.isfile(ref_dem):
                    parser.error("Reference DEM file does not exist: {}".format(tile_def))

        if water_tile_dir != '':
            if auto_select_arcticdem_water_tile_dir:
                if tile in watermask_tiles_visnav_need_editing:
                    parser.error("Tile '{}' is in tilelist indicating the visnav watermask "
                                 "still needs fixing: {}".format(tile, watermask_tiles_visnav_need_editing_file))
                if tile in watermask_tiles_greenland:
                    water_tile_dir = os.path.join(water_tile_dir, 'howat_greenland/tiled_watermasks/')
                elif tile in watermask_tiles_needing_visnav:
                    water_tile_dir = os.path.join(water_tile_dir, 'visnav/tiled_watermasks/')
                else:
                    water_tile_dir = os.path.join(water_tile_dir, 'global_surface_water/tiled_watermasks/')
                tile_waterTileDir_dict[tile] = water_tile_dir
            tile_land_file = os.path.join(water_tile_dir, '{}_land.tif'.format(tile))
            tile_ice_file = os.path.join(water_tile_dir, '{}_ice.tif'.format(tile))
            tile_water_file = os.path.join(water_tile_dir, '{}_water.tif'.format(tile))
            if not (os.path.isfile(tile_land_file) or os.path.isfile(tile_water_file)):
                print("ERROR: Tile land/water mask file does not exist: {}".format(tile_water_file))
                tile_missing_watermask = True
                tiles_missing_watermask.append(tile)

        ## If output does not exist, add to task list
        tile_outdir = os.path.join(args.dstdir, tile)
        tile_stdir = os.path.join(args.dstdir, tile, 'subtiles')
        tile_stdir_exists = os.path.isdir(tile_stdir)
        if not args.chain_mst and not os.path.isdir(tile_stdir):
            if not args.dryrun:
                os.makedirs(tile_stdir)

        final_subtile_fp_10m = os.path.join(tile_stdir, '{}_10000_10m.mat'.format(tile))
        final_subtile_fp_2m = os.path.join(tile_stdir, '{}_10000_2m.mat'.format(tile))
        finfile_10m = "{}_10m.fin".format(tile_stdir)
        finfile_2m = "{}_2m.fin".format(tile_stdir)
        if args.make_10m_only:
            final_subtile_fp = final_subtile_fp_10m
            finfile = finfile_10m
        else:
            final_subtile_fp = final_subtile_fp_2m
            finfile = finfile_2m

        run_tile = True

        mst_complete = False
        if args.chain_mst:
            mst_complete = True
            mst_finfile_10m = os.path.join(tile_outdir, '{}_10m.fin'.format(tile))
            if not os.path.isfile(mst_finfile_10m):
                mst_complete = False
            if not args.make_10m_only:
                for q in quads:
                    mst_finfile_2m = os.path.join(tile_outdir, '{}_{}_2m.fin'.format(tile, q))
                    if not os.path.isfile(mst_finfile_2m):
                        mst_complete = False
                        break
            if mst_complete:
                if args.make_10m_only:
                    print("Tile {} seems complete (MST 10m finfile exists)".format(tile))
                else:
                    print("Tile {} seems complete (MST 10m and 2m quad finfiles exist)".format(tile))
                run_tile = False

        if run_tile:
            # if args.sort_fix or args.rerun_without_cleanup:
            if args.rerun_without_cleanup:
                pass

            elif args.rerun:
                print("Verifying tile {} BST results before rerun".format(tile))

                bst_complete = False
                # if os.path.isfile(final_subtile_fp) and not args.require_finfiles:
                if os.path.isfile(final_subtile_fp) and args.bypass_finfile_req:
                    print("Tile seems complete with BST step ({} exists)".format(os.path.basename(final_subtile_fp)))
                    bst_complete = True
                elif os.path.isfile(finfile):
                    print("Tile seems complete with BST step ({} exists)".format(os.path.basename(finfile)))
                    bst_complete = True
                elif os.path.isfile(finfile_2m):
                    print("Tile seems complete with BST step (2m finfile {} exists)".format(os.path.basename(finfile_2m)))
                    bst_complete = True

                # if tile_stdir_exists and not args.make_10m_only:
                #     ## Remove subtiles with only 10m version
                #     tile_outfiles_10m = glob.glob(os.path.join(tile_stdir, '{}_*10m.mat'.format(tile)))
                #     for outfile_10m in tile_outfiles_10m:
                #         outfile_2m = outfile_10m.replace('10m.mat', '2m.mat')
                #         if not os.path.isfile(outfile_2m):
                #             print("Removing 10m subtile missing 2m component: {}".format(os.path.basename(outfile_10m)))
                #             bst_complete = False
                #             if not args.dryrun:
                #                 os.remove(outfile_10m)

                if args.chain_mst and not mst_complete:
                    if bst_complete:
                        print("... but MST step is not complete, so tile will be run")
                elif bst_complete:
                    run_tile = False

            elif tile_stdir_exists and any(os.scandir(tile_stdir)):
                print("Subtiles exist, skipping BST step for tile {}".format(tile))
                run_tile = False

        if run_tile:
            number_incomplete_tiles += 1
            if tile_missing_watermask:
                incomplete_tiles_missing_watermask.append(tile)
            else:
                tiles_to_run.append(tile)
                if water_tile_dir not in waterTileDir_tilesToRun_dict:
                    waterTileDir_tilesToRun_dict[water_tile_dir] = []
                waterTileDir_tilesToRun_dict[water_tile_dir].append(tile)


    if len(tiles_to_run) != number_incomplete_tiles:
        number_run_tiles_text = "{} of {} incomplete".format(len(tiles_to_run), number_incomplete_tiles)
    else:
        number_run_tiles_text = "{}".format(len(tiles_to_run))

    tiles_missing_watermask.sort()
    if len(tiles_missing_watermask) > 0:
        print("\n!!! ERROR !!! {} tiles are missing watermasks:\n{}\n".format(
            len(tiles_missing_watermask), '\n'.join(tiles_missing_watermask))
        )
        if len(incomplete_tiles_missing_watermask) > 0:
            print("!!! ERROR !!! {} INCOMPLETE tiles are missing watermasks:\n{}\n".format(
                len(incomplete_tiles_missing_watermask), '\n'.join(incomplete_tiles_missing_watermask))
            )
        if not args.skip_missing_watermasks:
            if len(tiles_to_run) > 0:
                print("{} total incomplete tiles, {} incomplete tiles missing watermasks".format(
                    number_incomplete_tiles, len(incomplete_tiles_missing_watermask)
                ))
                print("Provide the --skip-missing-watermasks argument to skip those tiles missing watermasks and run the rest")
                print("\nCould run {} tiles that have watermasks".format(number_run_tiles_text))
            sys.exit(1)


    print("Running {} tiles".format(number_run_tiles_text))
    if len(tiles_to_run) == 0:
        sys.exit(0)
    sleep_seconds = 10
    print("Sleeping {} seconds before submission".format(sleep_seconds))
    time.sleep(sleep_seconds)

    if not args.dryrun and batch_job_submission:
        pbs_ran_tilelist = os.path.join(pbs_logdir, 'ran_tiles.txt')
        with open(pbs_ran_tilelist, 'w') as tilelist_fp:
            for tile in tiles:
                tilelist_fp.write(tile+'\n')


    task_success_rc = 3 if len(tiles_to_run) > 0 else -1
    matlab_cmd_success = None

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

        jobnum_total = int(math.ceil(len(tiles_to_run) / float(args.tasks_per_job)))
        jobnum_fmt = '{:0>'+str(min(3, len(str(jobnum_total))))+'}'

        for water_tile_dir, tiles_to_run in waterTileDir_tilesToRun_dict.items():
            for jobnum, task_start_idx in enumerate(range(0, len(tiles_to_run), args.tasks_per_job), 1):

                tile_bundle = tiles_to_run[task_start_idx:task_start_idx+args.tasks_per_job]
                if batch_job_submission:
                    tile = None
                    task_name = jobnum_fmt.format(jobnum)
                else:
                    tile = tiles_to_run[task_start_idx]
                    task_name = tile
                job_name = 'bst_{}'.format(task_name)
                job_outfile = os.path.join(pbs_logdir, task_name+'.out')
                job_errfile = os.path.join(pbs_logdir, task_name+'.err')

                arg_waterTileDir = None
                if auto_select_arcticdem_water_tile_dir:
                    arg_waterTileDir = water_tile_dir

                if batch_job_submission and len(tile_bundle) < args.tasks_per_job:
                    sched_addl_vars_inst = sched_addl_vars.replace(
                        'nodes={}'.format(args.tasks_per_job),
                        'nodes={}'.format(len(tile_bundle))
                    )
                else:
                    sched_addl_vars_inst = sched_addl_vars

                if args.pbs:
                    cmd = r""" {}qsub -N {} -v {}{}ARG_TILENAME={}{}{} {} {} {} "{}" """.format(
                        sched_presubmit_cmd+' ; ' if sched_presubmit_cmd != '' else '',
                        job_name,
                        'TILERUN_JOBSCRIPT="{}",'.format(tilerun_jobscript) if batch_job_submission else '',
                        'IN_PARALLEL=true,' if batch_job_submission else '',
                        '@'.join(tile_bundle) if batch_job_submission else tile,
                        ',ARG_WATERTILEDIR="{}"'.format(arg_waterTileDir) if arg_waterTileDir is not None else '',
                        ','+sched_addl_envvars if sched_addl_envvars != '' else '',
                        '-q {}'.format(args.queue) if args.queue is not None else '',
                        '-o "{}" -e "{}"'.format(job_outfile, job_errfile) if sched_specify_outerr_paths else '',
                        sched_addl_vars_inst,
                        batch_jobscript if batch_job_submission else tilerun_jobscript,
                    )

                elif args.slurm:
                    cmd = r""" {}sbatch -J {} -v {}{}ARG_TILENAME={}{}{} {} {} "{}" """.format(
                        sched_presubmit_cmd+' ; ' if sched_presubmit_cmd != '' else '',
                        job_name,
                        'TILERUN_JOBSCRIPT="{}",'.format(tilerun_jobscript) if batch_job_submission else '',
                        'IN_PARALLEL=true,' if batch_job_submission else '',
                        '@'.join(tile_bundle) if batch_job_submission else tile,
                        ',ARG_WATERTILEDIR="{}"'.format(arg_waterTileDir) if arg_waterTileDir is not None else '',
                        ','+sched_addl_envvars if sched_addl_envvars != '' else '',
                        '-o "{}" -e "{}"'.format(job_outfile, job_errfile) if sched_specify_outerr_paths else '',
                        sched_addl_vars_inst,
                        batch_jobscript if batch_job_submission else tilerun_jobscript,
                    )

                else:
                    cmd = r""" {}bash "{}" {} """.format(
                        'export ARG_WATERTILEDIR="{}"'.format(arg_waterTileDir) if arg_waterTileDir is not None else '',
                        jobscript_temp,
                        tile,
                    )

                print("{}, {}".format(jobnum, cmd))
                if not args.dryrun:
                    matlab_cmd_rc = subprocess.call(cmd, shell=True, cwd=(pbs_logdir if args.pbs else None))
                    if matlab_cmd_rc == 2 and matlab_cmd_success is not False:
                        matlab_cmd_success = True
                    else:
                        matlab_cmd_success = False

        if args.pbs or args.slurm:
            print("Submitted {} tiles to scheduler".format(number_run_tiles_text))
        else:
            print("Ran {} tiles".format(number_run_tiles_text))

    if len(tiles_missing_watermask) > 0:
        print("\n!!! WARNING !!! {} tiles are missing watermasks:\n{}\n".format(
            len(tiles_missing_watermask), '\n'.join(tiles_missing_watermask))
        )
        if len(incomplete_tiles_missing_watermask):
            print("!!! WARNING !!! {} INCOMPLETE tiles are missing watermasks:\n{}\n".format(
                len(incomplete_tiles_missing_watermask), '\n'.join(incomplete_tiles_missing_watermask))
            )

    if matlab_cmd_success is True:
        task_success_rc = 2
    elif matlab_cmd_success is False:
        task_success_rc = 1

    sys.exit(task_success_rc)



if __name__ == '__main__':
    main()
