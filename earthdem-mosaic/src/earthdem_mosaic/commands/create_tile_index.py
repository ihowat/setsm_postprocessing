import shlex
import subprocess
from pathlib import Path

import click

from earthdem_mosaic.commands._utils import EXISTING_DIR


@click.command(short_help="Create ")
@click.argument("stage-dir", nargs=1, type=EXISTING_DIR)
@click.argument("suffix", nargs=1, type=str)
@click.argument("tif-resolution", nargs=1, type=int)
@click.argument("src_nodata", nargs=1, type=int)
@click.option("--resampling-method", type=str, default="bilinear", show_default=True)
def mosaic_outputs(
        stage_dir: Path,
        suffix: str,
        tif_resolution: int,
        src_nodata: int,
        resampling_method: str,
) -> None:
    """Create a mosaics of all TIFs in a stage directory with the same suffix.

    The mosaic created will be written to the STAGE_DIR with the name STAGE_DIR + SUFFIX

        STAGE_DIR = '10-coregistration-debug'

        SUFFIX = '_2m_browse.tif'

        MOSAIC --> 10-coregistration-debug/10-coregistration-debug_2m_browse.vrt
    """
    if ".tif" not in suffix:
        click.echo("ERROR: This command only supports files with the extention '.tif'")
        exit(1)

    vrt = stage_dir / f"{stage_dir.name}{suffix.replace('.tif', '.vrt')}"

    tifs = [str(path) for path in stage_dir.glob(f"**/utm*{suffix}")]
    if not tifs:
        click.echo(f"ERROR: Did not find any TIFs matching the suffix '{suffix}'")
        exit(1)

    click.echo(f"Creating {vrt} from {len(tifs)} TIFs")

    cmd = "\n\t".join(
        [
            "gdalbuildvrt",
            f"-tr {tif_resolution} {tif_resolution}",
            "-tap",
            f"-srcnodata {src_nodata}",
            f"-r {resampling_method}",
            "-overwrite",
            f"{vrt}",
            f"{' '.join(tifs)}",
        ]
    )

    subprocess.run(shlex.split(cmd), check=True)
