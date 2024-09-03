import shlex
import subprocess
from itertools import chain
from pathlib import Path
from typing import Self

import click

from earthdem_mosaic.commands._utils import EXISTING_FILE
from earthdem_mosaic.config import Settings
from earthdem_mosaic.utm_zone import UtmZone


class BatchMergeTileBufferBase:
    def __init__(self, settings: Settings, utm_zone: UtmZone, tiles: Path):
        # This dict stores portions of the command. The keys represent independent
        # portions of the command such as "base", "slurm" "dryrun". This allows for
        # multiple calls to the same portion of the command builder to override the
        # previous calls to it. For example: cmd.dryrun().dryrun() will not result in
        # multiple '--dryrun' flags being included in the final command.
        self._cmd_parts: dict[[str], list[str]] = {}

        self.tile_index = (
            settings.WORKING_ZONES_DIR
            / f"{utm_zone}"
            / "tile_index_files"
            / "tileNeighborIndex_2m.mat"
        )

        self.tiles = tiles
        self.script = (
            settings.SETSM_POSTPROCESSING_PYTHON_DIR / "batch_mergeTileBuffer.py"
        )
        self.matlib = settings.SETSM_POSTPROCESSING_MATLAB_DIR
        self.utm_zone_working_dir = settings.WORKING_ZONES_DIR / f"{utm_zone}"
        self.ref_dem_path = (
            settings.REFERENCE_DEM_DIR
            / f"{utm_zone}"
            / "<supertile>_10m_cop30_wgs84.tif"
        )
        self.cover_tif_path = (
            settings.LANDCOVER_DIR
            / f"{utm_zone}"
            / "<supertile>_10m_esa_worldcover_2021.tif"
        )
        self.job_config_file = (
            settings.SETSM_POSTPROCESSING_PYTHON_DIR / "lib/jobscript_config.ini"
        )
        self.job_script = (
            settings.SETSM_POSTPROCESSING_PYTHON_DIR / "lib/jobscript_<scheduler>.sh"
        )
        self.job_log_dir = (
            settings.WORKING_ZONES_DIR / f"{utm_zone}" / "processing_logs"
        )
        self.task_bundle_dir = Path.home() / "scratch" / "task_bundles"

    def slurm(
        self,
        job_name_prefix: str,
        job_ncpus: int,
        job_mem: int,
        job_walltime: int,
        job_queue: str,
    ) -> Self:
        """Prepare the command for submission to the slurm cluster."""
        key = "slurm"
        value = [
            "--slurm",
            f"--job-config-file {self.job_config_file}",
            f"--job-script '{self.job_script}'",
            f"--job-name-prefix '{job_name_prefix}'",
            f"--job-log-dir {self.job_log_dir}",
            f"--job-ncpus {job_ncpus}",
            f"--job-mem {job_mem}",
            f"--job-walltime {job_walltime}",
            f"--job-queue {job_queue}",
            "--tasks-per-job 1",
            "--tasks-per-job-mode serial",
            f"--task-bundle-dir {self.task_bundle_dir}",
        ]
        self._cmd_parts.update({key: value})
        return self

    def dryrun(self) -> Self:
        """Preform a dryrun of the command."""
        key = "dryrun"
        value = ["--dryrun"]
        self._cmd_parts.update({key: value})
        return self

    def as_str(self, *, pretty: bool = False) -> str:
        """Export a string, either formatted or not."""
        sep = "\n\t" if pretty else " "
        parts = chain(*self._cmd_parts.values())
        return sep.join(parts)

    def as_list(self) -> list[str]:
        """Export as a list of strings that can be passed to subprocess.run()."""
        parts = chain(*self._cmd_parts.values())
        return shlex.split(" ".join(parts))


class MergeRow(BatchMergeTileBufferBase):
    def __init__(self, settings: Settings, utm_zone: UtmZone, tiles: Path):
        super().__init__(settings=settings, utm_zone=utm_zone, tiles=tiles)
        self.base_command()

    def base_command(self) -> Self:
        """Merge Row base command."""
        stage_dir = "00-matfiles"
        axis = "row"
        key = "base"
        value = [
            f"python {self.script}",
            f"{self.utm_zone_working_dir / stage_dir} {self.tiles} 2",
            f"--tile-index {self.tile_index}",
            "--tile-org pgc",
            f"--process-group {axis}",
            f"--matlib {self.matlib}",
        ]
        self._cmd_parts.update({key: value})
        return self


class MergeColumn(BatchMergeTileBufferBase):
    def __init__(self, settings: Settings, utm_zone: UtmZone, tiles: Path):
        super().__init__(settings=settings, utm_zone=utm_zone, tiles=tiles)
        self.base_command()

    def base_command(self) -> Self:
        """Merge Row base command."""
        stage_dir = "00-matfiles"
        axis = "column"
        key = "base"
        value = [
            f"python {self.script}",
            f"{self.utm_zone_working_dir / stage_dir} {self.tiles} 2",
            f"--tile-index {self.tile_index}",
            "--tile-org pgc",
            f"--process-group {axis}",
            f"--matlib {self.matlib}",
        ]
        self._cmd_parts.update({key: value})
        return self


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
    if axis == "row":
        cmd = MergeRow(settings=settings, utm_zone=utm_zone, tiles=tiles)
    if axis == "column":
        cmd = MergeColumn(settings=settings, utm_zone=utm_zone, tiles=tiles)
    if slurm:
        cmd.slurm(
            job_name_prefix=f"mtb_{axis}",
            job_ncpus=4,
            job_mem=16,
            job_walltime=48,
            job_queue="big_mem",
        )
    if dryrun:
        cmd.dryrun()
    if show_command:
        print(cmd.as_str(pretty=True))
        return

    subprocess.run(cmd.as_list(), check=True)
