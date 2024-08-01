import os
import shutil
from enum import StrEnum
from pathlib import Path
from typing import Optional

import click
import geopandas
from geopandas import GeoDataFrame

from earthdem_mosaic.commands._utils import EXISTING_FILE
from earthdem_mosaic.config import Settings
from earthdem_mosaic.utm_zone import UtmZone


class FinalProduct(StrEnum):
    BROWSE = "_browse.tif"
    COUNTMT = "_countmt.tif"
    COUNT = "_count.tif"
    DATAMASK = "_datamask.tif"
    DEM = "_dem.tif"
    FIN = ".fin"
    MAD = "_mad.tif"
    MAT = ".mat"
    MAXDATE = "_maxdate.tif"
    MINDATE = "_mindate.tif"
    META = "_meta.txt"


def get_correspoding_product_paths(product: Path, type: FinalProduct) -> list[Path]:
    """Given a path to a product, returns the paths of all products for that quartertile
    including the provided product.
    """
    return [product.parent / product.name.replace(type, p) for p in FinalProduct]


def get_selected_tiles(gdf: GeoDataFrame) -> list[Path]:
    yes_slope_filter = gdf.loc[gdf["skip_slope_filter"] == False][
        "path_yes_slope_filter"
    ].to_list()
    no_slope_filter = gdf.loc[gdf["skip_slope_filter"] == True][
        "path_no_slope_filter"
    ].to_list()
    return [Path(p) for p in [*yes_slope_filter, *no_slope_filter]]


def all_tiles_reviewed(gdf: GeoDataFrame) -> bool:
    return gdf["reviewed"].all()


def hardlink_file(src: Path, dst: Path) -> None:
    """Create a hardlink from src to dst, overwriting dst if it exists."""
    dst.unlink() if dst.exists() else None
    dst.parent.mkdir(exist_ok=True, parents=True)
    os.link(src=src, dst=dst)


def hardlink_tree(src: Path, dst: Path) -> None:
    """Create hardlinks for a directory (src) and its contents to dst. Does not
    overwrite dst files if they already exist.
    """
    if not src.exists():
        raise FileExistsError(f"{src} does not exist")
    if not src.is_dir():
        raise NotADirectoryError(f"{src} is not a directory")
    if dst.exists() and not dst.is_dir():
        raise NotADirectoryError(f"{dst} is not a directory")

    try:
        shutil.copytree(src, dst, copy_function=os.link)
    except FileExistsError:
        # Do not overwrite existing destination files
        pass


@click.command(short_help="Hardlink selected tiles to the final products directory.")
@click.option("-v", "--verbose", is_flag=True)
@click.option("--dryrun", is_flag=True, help="Print actions without executing")
@click.argument("utm_zone", nargs=1, type=UtmZone)
@click.argument("slope-filter-review", nargs=1, type=EXISTING_FILE)
@click.argument("skipreg-shapefile", nargs=1, type=EXISTING_FILE, required=False)
@click.pass_obj
def link_final_products(
        settings: Settings,
        verbose: bool,
        dryrun: bool,
        utm_zone: UtmZone,
        slope_filter_review: Path,
        skipreg_shapefile: Optional[Path],
) -> None:
    """Move the final matfiles and TIFs to the directory set in the
    EARTHDEM_MOSAIC_FINAL_PRODUCTS_DIR environment variable.

    WARNING: Overwrites the destination file if one already exists
    """
    click.echo("Determining source tiles from slope-filter-review")
    gdf = geopandas.read_file(slope_filter_review)
    if not all_tiles_reviewed(gdf):
        click.echo("Some tiles have not been reviewed (reviewed == False)")
        click.echo(
            "All tiles must indicate that they have been reviewed before proceeding."
        )
        click.echo("No files moved. Exiting...")
        exit(1)

    tiles = get_selected_tiles(gdf=gdf)

    click.echo(f"Linking products from {len(tiles)} to final products directory")
    for tile in tiles:
        products = get_correspoding_product_paths(
            product=tile, type=FinalProduct.BROWSE
        )

        for src in products:
            if not src.exists():
                click.echo(f"Source file not found: {src}")
                continue
            dst = (
                    settings.FINAL_PRODUCTS_DIR / f"{utm_zone}" / src.parent.name / src.name
            )
            if dryrun or verbose:
                click.echo(f"rm {dst}") if dst.exists() else None
                click.echo(f"cp --link {src} {dst}")
            if not dryrun:
                hardlink_file(src, dst)

    click.echo("Linking processing_logs to final products directory")
    src = settings.WORKING_ZONES_DIR / f"{utm_zone}" / "processing_logs"
    dst = settings.FINAL_PRODUCTS_DIR / f"{utm_zone}" / "processing_logs"
    if dryrun or verbose:
        click.echo(f"cp -R --link {src} {dst}")
    if not dryrun:
        hardlink_tree(src, dst)

    click.echo("Linking SLOPE_FILTER_REVIEW file to final products directory")
    src = slope_filter_review
    dst = settings.FINAL_PRODUCTS_DIR / f"{utm_zone}" / slope_filter_review.name
    if dryrun or verbose:
        click.echo(f"rm {dst}") if dst.exists() else None
        click.echo(f"cp --link {src} {dst}")
    if not dryrun:
        hardlink_file(src, dst)

    if skipreg_shapefile:
        click.echo("Linking SKIPREG_SHAPEFILE files to final products directory")
        shapefile_parts = skipreg_shapefile.parent.glob(f"{skipreg_shapefile.stem}.*")
        for src in shapefile_parts:
            dst = settings.FINAL_PRODUCTS_DIR / f"{utm_zone}" / src.name
            if dryrun or verbose:
                click.echo(f"rm {dst}") if dst.exists() else None
                click.echo(f"cp --link {src} {dst}")
            if not dryrun:
                hardlink_file(src, dst)
