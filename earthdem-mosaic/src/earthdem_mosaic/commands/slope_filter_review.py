from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import click
import geopandas
import numpy as np
import pandas as pd
import rasterio
import shapely

from earthdem_mosaic.config import Settings
from earthdem_mosaic.utm_zone import UtmZone


@dataclass
class RasterInfo:
    path: Path
    x_resolution: float
    y_resolution: float
    total_cell_count: int
    nodata_cell_count: int
    crs: rasterio.CRS
    footprint: shapely.Polygon

    @classmethod
    def from_geotiff(cls: RasterInfo, geotiff: Path) -> RasterInfo:
        with rasterio.open(geotiff) as raster:
            return cls(
                path=geotiff,
                x_resolution=raster.res[0],
                y_resolution=raster.res[1],
                total_cell_count=raster.height * raster.width,
                nodata_cell_count=cls._get_nodata_cell_count(raster=raster),
                crs=raster.crs,
                footprint=cls._get_raster_footprint(raster=raster),
            )

    @staticmethod
    def _get_raster_footprint(raster: rasterio.DatasetReader) -> shapely.Polygon:
        bounds = raster.bounds
        return shapely.box(
            xmin=bounds.left,
            xmax=bounds.right,
            ymin=bounds.bottom,
            ymax=bounds.top,
        )

    @staticmethod
    def _get_nodata_cell_count(raster: rasterio.DatasetReader) -> int:
        # Create an array where nodata cells are stored 1 and all others as 0
        nodata_cells: np.ndarray = np.where(raster.read_masks(1) == 0, 1, 0)
        # Since each nodata cell is equal to 1, sum and count are the same
        return nodata_cells.sum()

    @property
    def quartertile(self) -> str:
        regex = re.compile(r"utm\d{2}[ns]_\d{2}_\d{2}_\d_\d")
        return regex.search(self.path.stem).group()

    @property
    def nodata_area(self) -> float:
        return self.x_resolution * self.y_resolution * self.nodata_cell_count

    @property
    def nodata_ratio(self) -> float:
        return self.nodata_cell_count / self.total_cell_count

    @property
    def nodata_percent(self) -> float:
        return round(self.nodata_ratio * 100, 2)


def raster_info_gdf(geotifs: list[Path]) -> geopandas.GeoDataFrame:
    rasters = [RasterInfo.from_geotiff(g) for g in geotifs]
    gdf_content = {
        "data": {
            "quartertile": [r.quartertile for r in rasters],
            "path": [str(r.path) for r in rasters],
            "nodata_percent": [r.nodata_percent for r in rasters],
        },
        "crs": rasters[0].crs,
        "geometry": [r.footprint for r in rasters],
    }
    return geopandas.GeoDataFrame(**gdf_content)


@click.command(short_help="Create tile index with review fields")
@click.option("-v", "--verbose", is_flag=True)
@click.argument("utm_zone", nargs=1, type=UtmZone)
@click.pass_obj
def slope_filter_review(settings: Settings, verbose: bool, utm_zone: UtmZone) -> None:
    """Create a geopackage with features representing the footprint of the tiles to be
    reviewed with fields for indicating review status and outcome.
    """
    zone_dir = settings.WORKING_ZONES_DIR / f"{utm_zone}"
    no_slope_filter_dir = zone_dir / "20-no-slope-filter"
    yes_slope_filter_dir = zone_dir / "30-yes-slope-filter"
    review_features = zone_dir / f"{utm_zone}_slope_filter_review.gpkg"

    if review_features.exists():
        print(f"Slope filter review features already exist at: {review_features}")
        return

    print("Building index of no-slope-filter TIFs")
    no_slope_filter = raster_info_gdf(
        geotifs=list(no_slope_filter_dir.glob("**/*_browse.tif"))
    )
    print("Building index of yes-slope-filter TIFs")
    yes_slope_filter = raster_info_gdf(
        geotifs=list(yes_slope_filter_dir.glob("**/*_browse.tif"))
    )

    print("Merging indexes")
    gdf = pd.merge(
        left=no_slope_filter,
        right=yes_slope_filter.drop("geometry", axis="columns"),
        how="left",
        on="quartertile",
        suffixes=("_no_slope_filter", "_yes_slope_filter"),
    )

    # Compute change in nodata percentage
    gdf["nodata_percent_delta"] = (
        gdf["nodata_percent_yes_slope_filter"] - gdf["nodata_percent_no_slope_filter"]
    ).round(decimals=2)

    # Add review fields with default values
    gdf["reviewed"] = False
    gdf["skip_slope_filter"] = False

    EXCLUDE_TILE_THRESHOLD = 70.0  # percent
    gdf["exclude_tile"] = (
        gdf["nodata_percent_yes_slope_filter"] > EXCLUDE_TILE_THRESHOLD
    )

    # Reorder fields and sort rows by quartertile name
    gdf = gdf[
        [
            "quartertile",
            "path_no_slope_filter",
            "path_yes_slope_filter",
            "nodata_percent_no_slope_filter",
            "nodata_percent_yes_slope_filter",
            "nodata_percent_delta",
            "reviewed",
            "skip_slope_filter",
            "exclude_tile",
            "geometry",
        ]
    ]
    gdf = gdf.sort_values(by="quartertile", axis="index")

    print(f"Writing slope filter review features to: {review_features}")
    gdf.to_file(filename=review_features)
    # Set the permissions of the slope_filter_review geopackage to -rw-rw-r--
    review_features.chmod(0o664)
