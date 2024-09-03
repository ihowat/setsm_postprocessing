import shlex
import subprocess
from itertools import chain
from pathlib import Path
from typing import Self

import click

from earthdem_mosaic.commands._utils import EXISTING_FILE
from earthdem_mosaic.config import Settings
from earthdem_mosaic.utm_zone import UtmZone


class BatchTiles2TifV4Base:
    def __init__(self, settings: Settings, utm_zone: UtmZone, tiles: Path):
        # This dict stores portions of the command. The keys represent independent
        # portions of the command such as "base", "slurm" "dryrun". This allows for
        # multiple calls to the same portion of the command builder to override the
        # previous calls to it. For example: cmd.dryrun().dryrun() will not result in
        # multiple '--dryrun' flags being included in the final command.
        self._cmd_parts: dict[[str], list[str]] = {}

        self.tiles = tiles
        self.script = settings.SETSM_POSTPROCESSING_PYTHON_DIR / "batch_tiles2tif_v4.py"
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

    def keep_old_results(self) -> Self:
        """Keep existing files that would be created by this command.

        Mutually exclusive with `delete_old_results()`
        """
        key = "old_results"
        value = ["--keep-old-results output-set"]
        self._cmd_parts.update({key: value})
        return self

    def delete_old_results(self) -> Self:
        """Delete existing files that would be created by this command.

        Mutually exclusive with `keep_old_results()`
        """
        key = "old_results"
        value = ["--rerun"]
        self._cmd_parts.update({key: value})
        return self

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


class CoregDebug(BatchTiles2TifV4Base):
    def __init__(self, settings: Settings, utm_zone: UtmZone, tiles: Path):
        super().__init__(settings=settings, utm_zone=utm_zone, tiles=tiles)
        self.base_command()

    def base_command(self) -> Self:
        """Coregistration Debug base command."""
        stage_dir = "10-coregistration-debug"

        key = "base"
        value = [
            f"python {self.script}",
            f"{self.utm_zone_working_dir / stage_dir} {self.tiles} earthdem 2",
            "--tif-format cog",
            "--output-set coreg-debug",
            "--tile-buffer-meters 100",
            f"--ref-dem-path '{self.ref_dem_path}'",
            f"--cover-tif-path '{self.cover_tif_path}'",
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
            f"--matlib {self.matlib}",
        ]
        self._cmd_parts.update({key: value})
        return self


class NoSlopeFilter(BatchTiles2TifV4Base):
    def __init__(self, settings: Settings, utm_zone: UtmZone, tiles: Path):
        super().__init__(settings=settings, utm_zone=utm_zone, tiles=tiles)
        self.base_command()

    def base_command(self) -> Self:
        """No Slope Filter base command."""
        stage_dir = "20-no-slope-filter"

        key = "base"
        value = [
            f"python {self.script}",
            f"{self.utm_zone_working_dir / stage_dir} {self.tiles} earthdem 2",
            "--tif-format cog",
            "--output-set full",
            "--tile-buffer-meters 100",
            f"--ref-dem-path '{self.ref_dem_path}'",
            f"--cover-tif-path '{self.cover_tif_path}'",
            "--apply-ref-filter false",
            "--apply-water-fill false",
            "--fill-water-interp-method 2",
            "--register-to-ref none",
            "--add-sea-surface-height false",
            "--use-final-qc-mask false",
            "--tile-org pgc",
            "--process-by tile-file",
            f"--matlib {self.matlib}",
        ]
        self._cmd_parts.update({key: value})
        return self


class YesSlopeFilter(BatchTiles2TifV4Base):
    def __init__(self, settings: Settings, utm_zone: UtmZone, tiles: Path):
        super().__init__(settings=settings, utm_zone=utm_zone, tiles=tiles)
        self.base_command()

    def base_command(self) -> Self:
        """Yes Slope Filter base command."""
        stage_dir = "30-yes-slope-filter"

        key = "base"
        value = [
            f"python {self.script}",
            f"{self.utm_zone_working_dir / stage_dir} {self.tiles} earthdem 2",
            "--tif-format cog",
            "--output-set full",
            "--tile-buffer-meters 100",
            f"--ref-dem-path '{self.ref_dem_path}'",
            f"--cover-tif-path '{self.cover_tif_path}'",
            "--apply-ref-filter true",
            "--apply-water-fill false",
            "--fill-water-interp-method 2",
            "--register-to-ref none",
            "--add-sea-surface-height false",
            "--use-final-qc-mask false",
            "--tile-org pgc",
            "--process-by tile-file",
            f"--matlib {self.matlib}",
        ]
        self._cmd_parts.update({key: value})
        return self


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
    cmd = CoregDebug(settings=settings, utm_zone=utm_zone, tiles=tiles)

    if slurm:
        cmd = cmd.slurm(
            job_name_prefix="t2t_coreg_debug",
            job_ncpus=2,
            job_mem=99,
            job_walltime=4,
            job_queue="big_mem",
        )

    if dryrun:
        cmd = cmd.dryrun()

    if show_command:
        print(cmd.as_str(pretty=True))
        return

    subprocess.run(cmd.as_list(), check=True)


@click.command(short_help="Export final TIFs with or without slope filtering")
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
    if apply_slope_filter:
        cmd = YesSlopeFilter(settings=settings, utm_zone=utm_zone, tiles=tiles)
    else:
        cmd = NoSlopeFilter(settings=settings, utm_zone=utm_zone, tiles=tiles)

    cmd = cmd.delete_old_results() if rerun else cmd.keep_old_results()

    if slurm and apply_slope_filter:
        cmd = cmd.slurm(
            job_name_prefix="t2t_yes_slope_filter",
            job_ncpus=2,
            job_mem=99,
            job_walltime=4,
            job_queue="big_mem",
        )
    if slurm and not apply_slope_filter:
        cmd = cmd.slurm(
            job_name_prefix="t2t_no_slope_filter",
            job_ncpus=2,
            job_mem=16,
            job_walltime=4,
            job_queue="big_mem",
        )
    if dryrun:
        cmd = cmd.dryrun()
    if show_command:
        print(cmd.as_str(pretty=True))
        return

    subprocess.run(cmd.as_list(), check=True)
