import os
from pathlib import Path

import click

from earthdem_mosaic.config import Settings
from earthdem_mosaic.utm_zone import UtmZone


@click.command(short_help="Hardlink source matfiles to working directory")
@click.option("-v", "--verbose", is_flag=True)
@click.option("--dryrun", is_flag=True, help="Print actions without executing")
@click.argument("utm_zone", nargs=1, type=UtmZone)
@click.pass_obj
def link_source_matfiles(
    settings: Settings, verbose: bool, dryrun: bool, utm_zone: UtmZone
) -> None:
    """Hardlink the 2m .mat files from Settings.MATFILE_SOURCE_DIR to the 00-matfiles
    working directory of the UTM zone.

    This assumes that the directory structure below Settings.MATFILE_SOURE_DIR is the
    same as that below the 00-matfiles working directory. The assumed structure is:

        Structure: {src/dst dir} / {supertile} / {quartertile files}
        Example:    00-matfiles / utm20n_08_05 / utm20n_08_05_2_1_2m.mat
    """
    click.echo(
        f"Searching {settings.MATFILE_SOURCE_DIR} for files with suffix: _2m.mat",
    )
    matfiles = settings.MATFILE_SOURCE_DIR.glob(pattern=f"**/{utm_zone}*_2m.mat")

    src_supertiles_dir = settings.MATFILE_SOURCE_DIR
    dst_supertiles_dir = settings.WORKING_ZONES_DIR / f"{utm_zone}" / "00-matfiles"

    for src_matfile in matfiles:
        dst_matfile = Path(
            str(src_matfile).replace(str(src_supertiles_dir), str(dst_supertiles_dir)),
        )

        if not dst_matfile.exists():
            if verbose:
                click.echo(f"cp --link {src_matfile} {dst_matfile}")
            if not dryrun:
                dst_matfile.parent.mkdir(parents=True, exist_ok=True)
                os.link(src=src_matfile, dst=dst_matfile)

        # Most, but not necessarily all, matfiles will have an adjacent .fin file
        # Matfile processing scripts check for these .fin files and output warnings
        # if they are not found. Link any source .fin files that exist.
        src_finfile = src_matfile.with_suffix(".fin")
        if src_finfile.exists():
            dst_finfile = dst_matfile.with_suffix(".fin")
            if not dst_finfile.exists():
                if verbose:
                    click.echo(f"cp --link {src_finfile} {dst_finfile}")
                if not dryrun:
                    dst_finfile.parent.mkdir(parents=True, exist_ok=True)
                    os.link(src=src_finfile, dst=dst_finfile)
