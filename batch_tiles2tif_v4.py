
import argparse
import glob
import os
import subprocess
import sys

from lib import jobscript_utils
from lib.jobscript_utils import wrap_multiline_str
from lib.jobscript_utils import UnsupportedMethodError


# This script paths
SCRIPT_FILE = os.path.abspath(os.path.realpath(__file__))
SCRIPT_FNAME = os.path.basename(SCRIPT_FILE)
SCRIPT_NAME, SCRIPT_EXT = os.path.splitext(SCRIPT_FNAME)
SCRIPT_DIR = os.path.dirname(SCRIPT_FILE)

# Paths relative to this script
MATLAB_LIBDIR = os.path.join(SCRIPT_DIR, "../setsm_postprocessing4")


# Argument defaults by 'domain'

domain_finalQcMask_dict = {
    'arcticdem': None,
    'earthdem': None,
    # 'rema': "/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/final_qc_mask/rema_v2/rema_final_mask.mat",
    'rema': "/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/final_qc_mask/rema_v2.1/rema_final_mask_v2_1.mat",
}
supertile_key = '<supertile>'
domain_refDemPath_dict = {
    'arcticdem': "/mnt/pgc/data/elev/dem/copernicus-dem-30m/mosaic/arctic_tiles_wgs84/<supertile>_10m_cop30_wgs84.tif",
    'earthdem': None,
    'rema': "/mnt/pgc/data/elev/dem/copernicus-dem-30m/mosaic/rema_tiles_wgs84/<supertile>_10m_cop30_wgs84.tif",
}
domain_coverTifPath_dict = {
    'arcticdem': "/mnt/pgc/data/thematic/landcover/esa_worldcover_2020/mosaics/arctic_tiles/<supertile>_10m_cover.tif",
    'earthdem': None,
    'rema': None,
}
domain_cop30SkipregShp_dict = {
    'arcticdem': "/mnt/pgc/data/elev/dem/setsm/ArcticDEM/mosaic/v4.1/arcticdem_v4.1_mosaic_reg-cop30_skip-reg.shp",
    'earthdem': None,
    'rema': None,
}


class RawTextArgumentDefaultsHelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter): pass

def get_arg_parser():
    
    parser = argparse.ArgumentParser(
        formatter_class=RawTextArgumentDefaultsHelpFormatter,
        description=wrap_multiline_str("""
            Export raster data from Matlab '<tilename>_<resolution>.mat'
            DEM mosaic tile files as GeoTIFF raster images, along with
            associated processing metadata in text format.
        """)
    )
    
    ## Positional args
    
    parser.add_argument(
        'tiledir',
        type=str,
        help="Root directory of tiles to be processed."
    )
    parser.add_argument(
        'tiles',
        help=wrap_multiline_str("""
            List of mosaic supertiles (10m tile names) to process,
            either specified on command line (comma delimited)
            or a text file list (each tile on separate line).
        """)
    )
    parser.add_argument(
        'domain',
        type=str,
        choices=['arcticdem', 'earthdem', 'rema'],
        help="DEM production domain of source tiles"
    )
    parser.add_argument(
        'resolution',
        type=int,
        choices=[10, 2],
        help=wrap_multiline_str("""
            Resolution class of existing tiles in source tiledir
            to be processed (meters).
        """)
    )

    ## Optional args

    parser.add_argument(
        '--tif-format',
        type=str,
        choices=['lzw', 'cog'],
        default='cog',
        help="Format of output GeoTIFF rasters."
    )
    parser.add_argument(
        '--output-set',
        type=str,
        choices=['full', 'dem-browse', 'dem', 'browse', 'meta', 'coreg-debug'],
        default='full',
        help=wrap_multiline_str(r"""
            Set of output tile result files to create. All settings include
            creation of the meta.txt file.
            \nIf 'browse' or 'meta', the default setting of --keep-old-results
            will become 'other' so that only the browse.tif or meta.txt file,
            will be built (or removed and rebuilt upon --rerun).
        """)
    )
    parser.add_argument(
        '--tile-nocrop',
        action='store_true',
        help=wrap_multiline_str("""
            Export full tile data arrays, without cropping to
            nearest multiple of 100km/50km (10m/2m tiles) in x/y
            coordinate values.
        """)
    )
    parser.add_argument(
        '--tile-buffer-meters',
        type=int,
        default=100,
        help=wrap_multiline_str("""
            Size of tile overlap buffer for exported rasters (meters),
            where the buffer is applied after cropping to nearest
            multiple of 100km/50km (10m/2m tiles) in x/y coordinate values.
        """)
    )
    parser.add_argument(
        '--ref-dem-path',
        type=str,
        default=None,
        help=wrap_multiline_str(r"""
            Path to reference DEM used for filtering with the --apply-ref-filter argument.
            DEM should cover all input tiles, or include the '{}' substring in path to be
            replaced with supertile name.
            \n(default is {})
        """.format(
            supertile_key,
            ', '.join(["{} if domain={}".format(val, dom) for dom, val in domain_refDemPath_dict.items()]
        )))
    )
    parser.add_argument(
        '--cover-tif-path',
        type=str,
        default=None,
        help=wrap_multiline_str(r"""
            Path to ESA WorldCover landcover raster used for filling water bodies with the reference DEM
            if the --apply-water-fill option is provided. The path should include the '{}'
            substring to be replaced with supertile name.
            \n(default is {})
        """.format(
            supertile_key,
            ', '.join(["{} if domain={}".format(val, dom) for dom, val in domain_coverTifPath_dict.items()]
        )))
    )
    parser.add_argument(
        '--apply-ref-filter',
        type=str,
        choices=['true', 'false'],
        default='false',
        help=wrap_multiline_str("""
            Apply Ian's slopeDifferenceFilter to tile DEM data in memory before tif export,
            using the reference DEM and landcover raster from the --ref-dem-path
            and --cover-tif-path arguments.
        """)
    )
    parser.add_argument(
        '--apply-topo-filter',
        type=str,
        choices=['true', 'false'],
        default='false',
        help=wrap_multiline_str("""
            Apply Ian's Residual Topography Fractional Difference Filter to tile
            DEM data in memory before tif export using the reference DEM and landcover
            raster from the --ref-dem-path and --cover-tif-path arguments.
        """)
    )
    parser.add_argument(
        '--apply-water-fill',
        type=str,
        choices=['true', 'false'],
        default='false',
        help=wrap_multiline_str("""
            Apply Ian's fillWater to tile DEM data in memory before tif export,
            using the reference DEM and landcover raster from the --ref-dem-path
            and --cover-tif-path arguments.
        """)
    )
    parser.add_argument(
        '--fill-water-interp-method',
        type=int,
        choices=list(range(6)),
        default=2,
        help=wrap_multiline_str("""
            Interpolation method to use for fillWater's inpaint_nans call.
        """)
    )

    parser.add_argument(
        '--register-to-ref',
        type=str,
        choices=['none', 'tiles', 'blobs', 'reportOffsetOnly'],
        default='none',
        help=wrap_multiline_str("""
            Register DEM to reference DEM. Requires both reference DEM and
            landcover raster from the --ref-dem-path and --cover-tif-path arguments.
        """)
    )
    parser.add_argument(
        '--register-to-ref-debug',
        action='store_true',
        help=wrap_multiline_str("""
            (Only applicable when --register-to-ref value other than 'none' is provided)
            Only write out a debug version of the registered DEM at 32m resolution,
            with no filter or fill applied, and then exit.
        """)
    )
    parser.add_argument(
        '--register-cop30-skipreg-shp',
        type=str,
        default=None,
        help=wrap_multiline_str(r"""
            (Only applicable when --register-to-ref value other than 'none' is provided)
            \nPath to shapefile of polygon features with 'skipreg' field value set to 1
            where, in the COP30 registration process, DEM data blobs should not be
            vertically adjusted.
            \n(default is {})
        """.format(
            supertile_key,
            ', '.join(["{} if domain={}".format(val, dom) for dom, val in domain_cop30SkipregShp_dict.items()]
        )))
    )
    
    parser.add_argument(
        '--add-sea-surface-height',
        type=str,
        choices=['true', 'false'],
        default='true',
        help=wrap_multiline_str("""
            In tile dem.tif area where ocean pixels ('land' mask is false) would
            normally be set to NoData, set elevation values to EGM2008 height above the
            WGS84 ellipsoid. 
        """)
    )
    parser.add_argument(
        '--final-qc-mask',
        type=str,
        default=None,
        help=wrap_multiline_str(r"""
            Domain-wide QC mask matfile in a format created by Ian for the
            REMA v2 final QC mask application step. Mask polygons in the matfile
            are designated for filling with either sea surface height or NoData.
            \n(default is {})
        """.format(
            ', '.join(["{} if domain={}".format(val, dom) for dom, val in domain_finalQcMask_dict.items()]
        )))
    )
    parser.add_argument(
        '--use-final-qc-mask',
        type=str,
        choices=['true', 'false'],
        default='true',
        help=wrap_multiline_str("""
            Whether or not to do the final QC mask application step, using the
            QC mask matfile provided through argument `--final-qc-mask`.
        """)
    )
    
    parser.add_argument(
        '--rerun',
        action='store_true',
        help=wrap_multiline_str("""
            Submit processing jobs even if output-set results files
            already exist.
        """)
    )
    parser.add_argument(
        '--keep-old-results',
        type=str,
        choices=['output-set', 'other', 'meta', 'all', 'none'],
        nargs='+',
        default=[],
        help=wrap_multiline_str(r"""
            Do not remove these classes of existing results files before
            submitting processing jobs. 'other' refers to tif results files
            outside of the provided --output-set setting.
            \nIf --output-set is 'browse' or 'meta', the default setting
            of --keep-old-results will become 'other' so that only the
            browse.tif or meta.txt file will be built
            (or removed and rebuilt upon --rerun).
        """)
    )

    parser.add_argument(
        '--tile-org',
        type=str,
        choices=['pgc', 'osu'],
        default='pgc',
        help="Tile file organization scheme."
    )
    parser.add_argument(
        '--process-by',
        type=str,
        choices=['supertile-dir', 'tile-file'],
        default='tile-file',
        help=wrap_multiline_str(r"""
            \nIf 'supertile-dir', if one tile within the supertile folder
            needs to be processed, hand whole supertile folder to Matlab script
            to be processed as a single task.
            \nIf 'tile-file', send each tile to be processed individually
            to Matlab script as a single task.
            \n
        """)
    )
    
    parser.add_argument(
        '--matlib',
        type=str,
        default=MATLAB_LIBDIR,
        help=wrap_multiline_str("""
            Path to directory containing necessary Matlab scripts
            and functions.
        """)
    )
    
    parser.add_argument(
        '--dryrun',
        action='store_true',
        help="Print actions without executing."
    )
    
    ## Argument groups
    
    jobscript_utils.argparse_add_job_scheduler_group(parser, config_group=SCRIPT_FNAME)

    return parser


def main():

    ## Parse and adjust arguments

    arg_parser = get_arg_parser()
    script_args = arg_parser.parse_args()

    root_tiledir = os.path.abspath(script_args.tiledir)

    if script_args.tiles.lower().endswith(('.txt', '.csv')) or os.path.isfile(script_args.tiles):
        tilelist_file = script_args.tiles
        if not os.path.isfile(script_args.tiles):
            arg_parser.error("'tiles' argument tilelist file does not exist: {}".format(tilelist_file))
        with open(tilelist_file, 'r') as tilelist_fp:
            supertile_list = [line for line in tilelist_fp.read().splitlines() if line != '']
    else:
        supertile_list = script_args.tiles.split(',')
    supertile_list = sorted(list(set(supertile_list)))

    if script_args.domain == 'arcticdem':
        global_projstr = 'polar stereo north'
    elif script_args.domain == 'earthdem':
        global_projstr = None
    elif script_args.domain == 'rema':
        global_projstr = 'polar stereo south'
    else:
        raise UnsupportedMethodError(
            "No projstr mapping for argument 'domain': {}".format(script_args.domain)
        )

    if script_args.ref_dem_path is None:
        script_args.ref_dem_path = domain_refDemPath_dict[script_args.domain]
    if script_args.cover_tif_path is None:
        script_args.cover_tif_path = domain_coverTifPath_dict[script_args.domain]
    if script_args.final_qc_mask is None:
        script_args.final_qc_mask = domain_finalQcMask_dict[script_args.domain]
    if script_args.register_cop30_skipreg_shp is None:
        script_args.register_cop30_skipreg_shp = domain_cop30SkipregShp_dict[script_args.domain]

    if script_args.output_set in ('browse', 'meta') and not script_args.keep_old_results:
        script_args.keep_old_results = ['other']

    res_name = '{}m'.format(script_args.resolution)

    if script_args.resolution == 2:
        quadname_list = ['1_1', '1_2', '2_1', '2_2']
    else:
        quadname_list = ['']

    output_set_to_matscript_arg_dict = {
        'full': 'full',
        'dem-browse': 'demAndBrowse',
        'dem': 'demOnly',
        'browse': 'browseOnly',
        'meta': 'metaOnly',
        'coreg-debug': 'demAndBrowse',
    }

    single_t2t_args = wrap_multiline_str(f"""
        , 'resolution','{res_name}'
        , 'outFormat','{script_args.tif_format.upper()}'
        , 'outSet','{output_set_to_matscript_arg_dict[script_args.output_set]}'
        , 'bufferMeters',{script_args.tile_buffer_meters}
        , 'registerToRef','{script_args.register_to_ref}'
    """)
    if script_args.tile_nocrop:
        single_t2t_args += ", 'noCrop'"
    if script_args.add_sea_surface_height == 'true':
        single_t2t_args += ", 'addSeaSurface'"
    if script_args.use_final_qc_mask == 'true' and script_args.final_qc_mask is not None:
        single_t2t_args += ", 'qcMaskFile','{}'".format(script_args.final_qc_mask)
    if script_args.register_to_ref_debug:
        single_t2t_args += ", 'registerToRefDebug'"
    if script_args.register_cop30_skipreg_shp:
        single_t2t_args += ", 'registerBlobsSkipregShp','{}'".format(script_args.register_cop30_skipreg_shp)
    single_t2t_args += ", 'fillWaterInterpMethod',{}".format(script_args.fill_water_interp_method)

    batch_t2t_args = single_t2t_args
    if script_args.output_set == 'meta':
        batch_t2t_args += ", 'metaOnly'"

    jobscript_utils.adjust_args(script_args, arg_parser)
    jobscript_utils.create_dirs(script_args, arg_parser)


    ## Verify arguments

    if not os.path.isdir(root_tiledir):
        arg_parser.error("Argument 'tiledir' is not an existing directory: {}".format(root_tiledir))
    if script_args.apply_ref_filter == 'true' and not script_args.ref_dem_path:
        arg_parser.error("--ref-dem-path cannot be empty when --apply-ref-filter=true")
    if script_args.apply_topo_filter == 'true' and not script_args.ref_dem_path:
        arg_parser.error("--ref-dem-path cannot be empty when --apply-topo-filter=true")
    if script_args.apply_water_fill == 'true' and (not script_args.ref_dem_path or not script_args.cover_tif_path):
        arg_parser.error("--ref-dem-path and --cover-tif-path cannot be empty when --apply-water-fill=true")
    if script_args.register_to_ref == 'none' and script_args.register_to_ref_debug:
        arg_parser.error("--register-to-ref-debug cannot be provided when --register-to-ref='none'")
    if script_args.register_to_ref_debug and script_args.output_set != 'coreg-debug':
        arg_parser.error("--output-set must be 'coreg-debug' if --register-to-ref-debug is used")
    if script_args.tile_org == 'osu' and script_args.process_by == 'supertile-dir':
        arg_parser.error("--process-by must be set to to 'tile-file' when --tile-org='osu'")

    if script_args.apply_ref_filter == 'true' or script_args.apply_water_fill == 'true' or script_args.apply_topo_filter == 'true':
        check_dir, _, _ = script_args.ref_dem_path.partition(supertile_key)
        check_dir = os.path.dirname(check_dir)
        if not os.path.isdir(check_dir):
            arg_parser.error("--ref-dem-path directory does not exist: {}".format(check_dir))
        if script_args.apply_water_fill:
            check_dir, _, _ = script_args.cover_tif_path.partition(supertile_key)
            check_dir = os.path.dirname(check_dir)
            if not os.path.isdir(check_dir):
                arg_parser.error("--cover-tif-path directory does not exist: {}".format(check_dir))

    jobscript_utils.verify_args(script_args, arg_parser)


    ## Main processing loop

    num_supertiles = len(supertile_list)
    if script_args.resolution == 2 and script_args.process_by == 'tile-file':
        est_num_tasks = num_supertiles * 4
    else:
        est_num_tasks = num_supertiles

    job_handler = jobscript_utils.JobHandler(
        script_args,
        est_num_tasks,
        num_tasks_is_estimate=True,
        init_env_requests='matlab gdal'
    )

    run_tile_list_all = []

    for supertile in supertile_list:

        if supertile.startswith('utm'):
            if script_args.domain != 'earthdem':
                arg_parser.error("domain should be 'earthdem' when 'utm*' prefix tilenames are provided")
        else:
            if script_args.domain == 'earthdem':
                arg_parser.error("domain should NOT be 'earthdem' when tilenames do not have 'utm*' prefix")

        tile_projstr = global_projstr
        if tile_projstr is None:
            assert script_args.domain == 'earthdem'

            utm_tilename_parts = supertile.split('_')
            utm_tilename_prefix = utm_tilename_parts[0]
            if not utm_tilename_prefix.startswith('utm'):
                arg_parser.error("Expected only UTM tile names (e.g. 'utm10n_01_01'), but got '{}'".format(supertile))

            tile_projstr = utm_tilename_prefix

        supertile_args = ''
        if script_args.ref_dem_path:
            ref_dem_file = script_args.ref_dem_path.replace(supertile_key, supertile)
            supertile_args += ", 'refDemFile','{}'".format(ref_dem_file)
        if script_args.cover_tif_path:
            cover_tif_file = script_args.cover_tif_path.replace(supertile_key, supertile)
            supertile_args += ", 'waterMaskFile','{}'".format(cover_tif_file)
        if script_args.apply_ref_filter == 'true':
            supertile_args += ", 'applySlopeDiffFilt'"
        if script_args.apply_topo_filter == 'true':
            supertile_args += ", 'applyResidualTopographyFractionalDifferenceFilter'"
        if script_args.apply_water_fill == 'true':
            supertile_args += ", 'applyWaterFill'"

        run_tile_matlist = []

        for quadname in quadname_list:
            tile_name = '{}_{}'.format(supertile, quadname) if quadname != '' else supertile
            tile_rootpath = os.path.join(
                root_tiledir,
                supertile if script_args.tile_org == 'pgc' else '',
                '{}_{}'.format(tile_name, res_name)
            )

            unregmatfile    = '{}.mat'.format(tile_rootpath)
            regmatfile      = '{}_reg.mat'.format(tile_rootpath)
            finfp           = '{}.fin'.format(tile_rootpath)
            demfp           = '{}_dem.tif'.format(tile_rootpath)
            browsefp        = '{}_browse.tif'.format(tile_rootpath)
            metafp          = '{}_meta.txt'.format(tile_rootpath)
            matfile         = regmatfile if os.path.isfile(regmatfile) else unregmatfile

            fp_mad      = '{}_mad.tif'.format(tile_rootpath)
            fp_count    = '{}_count.tif'.format(tile_rootpath)
            fp_countmt  = '{}_countmt.tif'.format(tile_rootpath)
            fp_mindate  = '{}_mindate.tif'.format(tile_rootpath)
            fp_maxdate  = '{}_maxdate.tif'.format(tile_rootpath)
            fp_datamask = '{}_datamask.tif'.format(tile_rootpath)

            coreg_debug_offset = '{}_dem_debug-reg_{}_offset.tif'.format(tile_rootpath, script_args.register_to_ref)
            coreg_debug_dem = '{}_dem_debug-reg_{}.tif'.format(tile_rootpath, script_args.register_to_ref)
            coreg_debug_demdiff = '{}_demdiff_debug-reg_{}.tif'.format(tile_rootpath, script_args.register_to_ref)

            if not os.path.isfile(unregmatfile) and not os.path.isfile(regmatfile):
                print("Tile {} {}m mat and reg.mat files do not exist{}: {}".format(
                    tile_name, script_args.resolution,
                    " (AND .fin file also does not exist!!)" if not os.path.isfile(finfp) else '',
                    matfile
                ))
                run_tile = False

            else:
                output_set_results_files_dict = {
                    'full':         [metafp, demfp, browsefp, fp_mad, fp_count, fp_countmt, fp_mindate, fp_maxdate, fp_datamask],
                    'dem-browse':   [metafp, demfp, browsefp],
                    'dem':          [metafp, demfp],
                    'browse':       [metafp, browsefp],
                    'meta':         [metafp],
                    'coreg-debug':  [coreg_debug_dem, coreg_debug_offset, coreg_debug_demdiff]
                }
                results_fp_list = output_set_results_files_dict[script_args.output_set]
                results_fp_exist_set = set([fp for fp in results_fp_list if os.path.isfile(fp)])

                if not script_args.rerun and len(results_fp_list) == len(results_fp_exist_set):
                    print("Tile {}: All results files exist, skipping".format(matfile))
                    run_tile = False

                else:
                    run_tile = True
                    results_fp_remove_set = set()

                    if 'all' in script_args.keep_old_results:
                        pass
                    else:
                        results_fp_remove_set.update(results_fp_exist_set)
                        if 'other' not in script_args.keep_old_results:
                            dstfps_old_pattern = [
                                metafp,
                                demfp.replace('_dem.tif', '*.tif'),
                            ]
                            dstfps_old = set([fp for pat in dstfps_old_pattern for fp in glob.glob(pat)])
                            results_fp_remove_set.update(dstfps_old)
                        if 'output-set' in script_args.keep_old_results:
                            results_fp_remove_set = results_fp_remove_set.difference(results_fp_exist_set)
                        if 'meta' in script_args.keep_old_results and metafp in results_fp_remove_set:
                            results_fp_remove_set.remove(metafp)

                    if results_fp_remove_set:
                        results_fp_remove_list = sorted(list(results_fp_remove_set))
                        print("Tile {}: Removing existing results files{}:\n  {}".format(
                            matfile,
                            ' (dryrun)' if script_args.dryrun else '',
                            '\n  '.join(results_fp_remove_list)
                        ))
                        if not script_args.dryrun:
                            for fp in results_fp_remove_list:
                                os.remove(fp)

            if run_tile:
                run_tile_matlist.append(matfile)

        run_tile_list = []

        if script_args.process_by == 'tile-file':
            run_tile_list = run_tile_matlist
        elif script_args.process_by == 'supertile-dir' and len(run_tile_matlist) > 0:
            supertile_dir = os.path.join(root_tiledir, supertile)
            run_tile_list = [supertile_dir]

        run_tile_list_all.extend(run_tile_list)
        # Un-comment the following line and un-indent the proceeding
        # `for` loop if you want to check & gather the full list of
        # processing tasks before looping through list for job submission.
    # for tile_path in run_tile_list_all:
        for tile_path in run_tile_list:
            job_id = None

            matlab_addpath = wrap_multiline_str(f"""
                addpath('{SCRIPT_DIR}');
                addpath('{script_args.matlib}');
            """)

            if not tile_path.endswith('.mat'):
                supertile_dir = tile_path
                job_id = os.path.basename(supertile_dir)

                task_cmd = jobscript_utils.matlab_cmdstr_to_shell_cmdstr(wrap_multiline_str(f"""
                    {matlab_addpath}
                    batch_tiles2tif_v4('{supertile_dir}', '{tile_projstr}'{batch_t2t_args}{supertile_args});
                """))

            else:
                tile_matfile = tile_path
                job_id, _ = os.path.splitext(os.path.basename(tile_matfile))

                task_cmd = jobscript_utils.matlab_cmdstr_to_shell_cmdstr(wrap_multiline_str(f"""
                    {matlab_addpath}
                    writeTileToTifv4('{tile_matfile}', '{tile_projstr}'{single_t2t_args}{supertile_args});
                    tileMetav4('{tile_matfile}');
                """))

            submit_cmd = job_handler.add_task_cmd(task_cmd, job_id)
            if submit_cmd is not None:
                print("Submitting job [{}/{}(est)]:\n  {}".format(
                    job_handler.job_num, job_handler.num_jobs, submit_cmd
                ))
                if not script_args.dryrun:
                    subprocess.call(submit_cmd, shell=True)

    if job_handler.bundle_tasks:
        submit_cmd = job_handler.get_last_bundle_submit_cmd()
        if submit_cmd is not None:
            print("Submitting job [{}/{}]:\n  {}".format(
                job_handler.job_num, job_handler.job_num, submit_cmd
            ))
            if not script_args.dryrun:
                subprocess.call(submit_cmd, shell=True)

    sys.exit(0)


if __name__ == '__main__':
    main()
