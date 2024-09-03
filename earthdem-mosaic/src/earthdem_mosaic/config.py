import os
from pathlib import Path

from pydantic_settings import BaseSettings, SettingsConfigDict


def find_env_file() -> Path:
    from_env_variable = os.getenv("EARTHDEM_MOSAIC_ENV_FILE")
    if from_env_variable:
        return Path(from_env_variable)
    # Use the .env file at the root of the package
    return Path(__file__).resolve().parents[2] / ".env"


class Settings(BaseSettings):
    MATFILE_SOURCE_DIR: Path
    """Directory of matfiles output from batch_mosaicSubTiles.py"""

    LANDCOVER_DIR: Path
    """Directory of pre-tiled ESA WorldCover TIFs"""

    REFERENCE_DEM_DIR: Path
    """Directory of pre-tiled Copernicus 30m TIFs"""

    FINAL_PRODUCTS_DIR: Path
    """Directory to save final production artifacts"""

    SETSM_POSTPROCESSING_PYTHON_DIR: Path
    """Directory of Python batch processing scripts wrapping matlab algorithms"""

    SETSM_POSTPROCESSING_MATLAB_DIR: Path
    """Directory of matlab algorithms"""

    # Working directories
    WORKING_ZONES_DIR: Path
    """Directory of in process UTM Zones"""

    # Configure automatic discovery and loading of settings from env file
    # The prefix must be included in the entries in the env file. Example:
    #   .env entry: EARTHDEM_MOSAIC_MATLAB_SOURCE_DIR=/some/path
    #   will set:   Settings.MATLAB_SOURCE_DIR
    model_config = SettingsConfigDict(
        # The .env file is expected to be at the root of the package
        env_file=find_env_file(),
        env_prefix="EARTHDEM_MOSAIC_",
    )
