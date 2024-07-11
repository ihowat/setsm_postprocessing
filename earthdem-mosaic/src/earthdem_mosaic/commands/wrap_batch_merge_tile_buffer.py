import shlex
import subprocess
from pathlib import Path

import click

from earthdem_mosaic.commands._utils import EXISTING_FILE
from earthdem_mosaic.config import Settings
from earthdem_mosaic.utm_zone import UtmZone


def _base_cmd(settings: Settings, tiles: Path, utm_zone: UtmZone, axis: str) -> str:
    tiledir = settings.WORKING_ZONES_DIR / str(utm_zone) / "00-matfiles"
    tile_index = (
        settings.WORKING_ZONES_DIR
        / str(utm_zone)
        / "tile_index_files/tileNeighborIndex_2m.mat"
    )

    return "\n\t".join(
        [
            f"python {settings.SETSM_POSTPROCESSING_PYTHON_DIR}/batch_mergeTileBuffer.py",
            f"{tiledir} {tiles} 2",
            f"--tile-index {tile_index}",
            "--tile-org pgc",
            f"--process-group {axis}",
            f"--matlib {settings.SETSM_POSTPROCESSING_MATLAB_DIR}",
        ],
    )


def _dryrun_option() -> str:
    return "--dryrun"


def _slurm_options(settings: Settings) -> str:
    job_config_file = (
        settings.SETSM_POSTPROCESSING_PYTHON_DIR / "lib/jobscript_config.ini"
    )
    job_script = (
        settings.SETSM_POSTPROCESSING_PYTHON_DIR / "lib/jobscript_<scheduler>.sh"
    )

    return "\n\t".join(
        [
            "--slurm",
            f"--job-config-file {job_config_file}",
            "--job-config-on",
            f"--job-script '{job_script}'",
            "--job-name-prefix 'mtb<resolution>m'",
            "--job-ncpus 4",
            "--job-mem 200",
            "--job-walltime 48",
            "--job-queue big_mem",
            "--tasks-per-job 1",
            "--tasks-per-job-mode serial",
            "--task-bundle-dir $HOME/scratch/task_bundles",
        ],
    )


@click.command(short_help="Merge overlapping matfile tile buffers")
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
@click.argument(
    "axis",
    type=click.Choice(["row", "column"]),
    nargs=1,
)
@click.pass_obj
def merge_buffers(
    settings: Settings,
    slurm: bool,
    dryrun: bool,
    show_command: bool,
    axis: str,
    utm_zone: UtmZone,
    tiles: Path,
) -> None:
    """Adjust overlapping tile buffers to merge adjacent tiles"""
    cmd = _base_cmd(settings=settings, tiles=tiles, utm_zone=utm_zone, axis=axis)
    if slurm:
        cmd = "\n\t".join([cmd, _slurm_options(settings=settings)])
    if dryrun:
        cmd = "\n\t".join([cmd, _dryrun_option()])
    if show_command:
        print(cmd)
        return

    subprocess.run(shlex.split(cmd), check=True)
