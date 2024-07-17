import shlex
import subprocess
from pathlib import Path

import click

from earthdem_mosaic.commands._utils import EXISTING_FILE
from earthdem_mosaic.config import Settings
from earthdem_mosaic.utm_zone import UtmZone


def _coreg_debug_base_cmd(settings: Settings, tiles: Path, utm_zone: UtmZone) -> str:
    tiledir = settings.WORKING_ZONES_DIR / str(utm_zone) / "10-coregistration-debug"
    ref_dem_path = (
        settings.REFERENCE_DEM_DIR / str(utm_zone) / "<supertile>_10m_cop30_wgs84.tif"
    )
    cover_tif_path = (
        settings.LANDCOVER_DIR
        / str(utm_zone)
        / "<supertile>_10m_esa_worldcover_2021.tif"
    )

    return "\n\t".join(
        [
            f"python {settings.SETSM_POSTPROCESSING_PYTHON_DIR}/batch_tiles2tif_v4.py",
            f"{tiledir} {tiles} earthdem 2",
            "--tif-format cog",
            "--output-set coreg-debug",
            "--tile-buffer-meters 100",
            f"--ref-dem-path '{ref_dem_path}'",
            f"--cover-tif-path '{cover_tif_path}'",
            "--apply-ref-filter false",
            "--apply-water-fill false",
            "--fill-water-interp-method 2",
            "--register-to-ref blobs",
            "--register-to-ref-debug",
            "--add-sea-surface-height false",
            "--use-final-qc-mask false",
            "--keep-old-results output-set",
            "--tile-org pgc",
            "--process-by tile-file",
            f"--matlib {settings.SETSM_POSTPROCESSING_MATLAB_DIR}",
        ],
    )


def _dryrun_option() -> str:
    return "--dryrun"


def _slurm_options(settings: Settings, job_name_prefix: str, utm_zone: UtmZone) -> str:
    job_config_file = (
        settings.SETSM_POSTPROCESSING_PYTHON_DIR / "lib/jobscript_config.ini"
    )
    job_script = (
        settings.SETSM_POSTPROCESSING_PYTHON_DIR / "lib/jobscript_<scheduler>.sh"
    )
    job_log_dir = settings.WORKING_ZONES_DIR / f"{utm_zone}" / "processing_logs"

    return "\n\t".join(
        [
            "--slurm",
            f"--job-config-file {job_config_file}",
            f"--job-script '{job_script}'",
            f"--job-name-prefix '{job_name_prefix}'",
            f"--job-log-dir {job_log_dir}",
            "--job-ncpus 2",
            "--job-mem 99",
            "--job-walltime 4",
            "--job-queue big_mem",
            "--tasks-per-job 1",
            "--tasks-per-job-mode serial",
            "--task-bundle-dir $HOME/scratch/task_bundles",
        ],
    )


@click.command(short_help="Create coregistration review TIFs")
@click.option("--slurm", is_flag=True, help="Submit jobs to slurm cluster")
@click.option("--dryrun", is_flag=True, help="Print actions without executing")
@click.option(
    "--show-command",
    is_flag=True,
    help="Print the generated command without executing",
)
@click.argument("utm_zone", nargs=1, type=UtmZone)
@click.argument(
    "tiles",
    nargs=1,
    type=EXISTING_FILE,
)
@click.pass_obj
def coreg_debug(
    settings: Settings,
    slurm: bool,
    dryrun: bool,
    show_command: bool,
    utm_zone: UtmZone,
    tiles: Path,
) -> None:
    """Produce coregistered DEMs and offset TIFs for review."""
    job_name_prefix = "t2t_coreg_debug"

    cmd = _coreg_debug_base_cmd(settings=settings, tiles=tiles, utm_zone=utm_zone)
    if slurm:
        slurm_options = _slurm_options(
            settings=settings, job_name_prefix=job_name_prefix, utm_zone=utm_zone
        )
        cmd = "\n\t".join([cmd, slurm_options])
    if dryrun:
        cmd = "\n\t".join([cmd, _dryrun_option()])
    if show_command:
        print(cmd)
        return

    subprocess.run(shlex.split(cmd), check=True)


def _export_tif_base_cmd(
    settings: Settings,
    tiles: Path,
    utm_zone: UtmZone,
    tiledir: Path,
    apply_slope_filter: bool,
) -> str:
    ref_dem_path = (
        settings.REFERENCE_DEM_DIR / str(utm_zone) / "<supertile>_10m_cop30_wgs84.tif"
    )
    cover_tif_path = (
        settings.LANDCOVER_DIR
        / str(utm_zone)
        / "<supertile>_10m_esa_worldcover_2021.tif"
    )

    return "\n\t".join(
        [
            f"python {settings.SETSM_POSTPROCESSING_PYTHON_DIR}/batch_tiles2tif_v4.py",
            f"{tiledir} {tiles} earthdem 2",
            "--tif-format cog",
            "--output-set full",
            "--tile-buffer-meters 100",
            f"--ref-dem-path '{ref_dem_path}'",
            f"--cover-tif-path '{cover_tif_path}'",
            f"--apply-ref-filter {str(apply_slope_filter).lower()}",
            "--apply-water-fill false",
            "--fill-water-interp-method 2",
            "--register-to-ref none",
            "--add-sea-surface-height false",
            "--use-final-qc-mask false",
            "--tile-org pgc",
            "--process-by tile-file",
            f"--matlib {settings.SETSM_POSTPROCESSING_MATLAB_DIR}",
        ],
    )


@click.command(short_help="Export final TIFs with or with slope filtering")
@click.option("--slurm", is_flag=True, help="Submit jobs to slurm cluster")
@click.option("--dryrun", is_flag=True, help="Print actions without executing")
@click.option(
    "--show-command",
    is_flag=True,
    help="Print the generated command without executing",
)
@click.option(
    "--apply-slope-filter",
    is_flag=True,
    help="Apply slope filter before exporting TIFs",
)
@click.option(
    "--rerun",
    is_flag=True,
    help="Remove existing outputs before running",
)
@click.argument("utm_zone", nargs=1, type=UtmZone)
@click.argument(
    "tiles",
    nargs=1,
    type=EXISTING_FILE,
)
@click.pass_obj
def export_final_tifs(
    settings: Settings,
    slurm: bool,
    dryrun: bool,
    show_command: bool,
    apply_slope_filter: bool,
    rerun: bool,
    utm_zone: UtmZone,
    tiles: Path,
) -> None:
    """Export final TIFs with or without slope filtering"""
    if not apply_slope_filter:
        job_name_prefix = "t2t_no_slope_filter"
        tiledir = settings.WORKING_ZONES_DIR / str(utm_zone) / "20-no-slope-filter"
    else:
        job_name_prefix = "t2t_yes_slope_filter"
        tiledir = settings.WORKING_ZONES_DIR / str(utm_zone) / "30-yes-slope-filter"

    cmd = _export_tif_base_cmd(
        settings=settings,
        tiles=tiles,
        utm_zone=utm_zone,
        tiledir=tiledir,
        apply_slope_filter=apply_slope_filter,
    )
    if rerun:
        cmd = "\n\t".join([cmd, "--rerun"])
    else:
        cmd = "\n\t".join([cmd, "--keep-old-results output-set"])
    if slurm:
        slurm_options = _slurm_options(
            settings=settings, job_name_prefix=job_name_prefix, utm_zone=utm_zone
        )
        cmd = "\n\t".join([cmd, slurm_options])
    if dryrun:
        cmd = "\n\t".join([cmd, _dryrun_option()])
    if show_command:
        print(cmd)
        return

    subprocess.run(shlex.split(cmd), check=True)
