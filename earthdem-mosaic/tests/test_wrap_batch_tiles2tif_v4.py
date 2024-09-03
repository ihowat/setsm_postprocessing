import os
import shlex
from pathlib import Path

from earthdem_mosaic.commands.wrap_batch_tiles2tif_v4 import NoSlopeFilter, YesSlopeFilter
from earthdem_mosaic.config import Settings
from earthdem_mosaic.utm_zone import UtmZone


def mock_settings() -> Settings:
    s = Settings()
    s.MATFILE_SOURCE_DIR = Path("/path/to/matfile_source_dir")
    s.LANDCOVER_DIR = Path("/path/to/landcover_dir")
    s.REFERENCE_DEM_DIR = Path("/path/to/reference_dem_dir")
    s.FINAL_PRODUCTS_DIR = Path("/path/to/final_products_dir")
    s.SETSM_POSTPROCESSING_PYTHON_DIR = (
        Path("/path/to/setsm_postprocessing_python_dir")
    )
    s.SETSM_POSTPROCESSING_MATLAB_DIR = (
        Path("/path/to/setsm_postprocessing_matlab_dir")
    )
    s.WORKING_ZONES_DIR = Path("/path/to/working_zones_dir")

    return s


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
            f"python {settings.SETSM_POSTPROCESSING_PYTHON_DIR / 'batch_tiles2tif_v4.py'}",
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


def test_no_slope_filter_base_command():
    settings = mock_settings()
    tiles = Path("/path/to/all_supertiles.txt")
    utm_zone = UtmZone.utm10n
    tiledir = Path(f"/path/to/working_zones_dir/{utm_zone}/20-no-slope-filter")
    apply_slope_filter = False

    expected = shlex.split(
        _export_tif_base_cmd(
            settings=settings,
            tiles=tiles,
            utm_zone=utm_zone,
            tiledir=tiledir,
            apply_slope_filter=apply_slope_filter,
        )
    )

    actual = NoSlopeFilter(settings=settings, utm_zone=utm_zone, tiles=tiles).as_list()

    assert actual == expected

def test_yes_slope_filter_base_command():
    settings = mock_settings()
    tiles = Path("/path/to/all_supertiles.txt")
    utm_zone = UtmZone.utm10n
    tiledir = Path(f"/path/to/working_zones_dir/{utm_zone}/30-yes-slope-filter")
    apply_slope_filter = True

    expected = shlex.split(
        _export_tif_base_cmd(
            settings=settings,
            tiles=tiles,
            utm_zone=utm_zone,
            tiledir=tiledir,
            apply_slope_filter=apply_slope_filter,
        )
    )

    actual = YesSlopeFilter(settings=settings, utm_zone=utm_zone, tiles=tiles).as_list()

    assert actual == expected
