from pathlib import Path

from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    # Source and reference data
    MATLAB_SOURCE_DIR: Path
    LANDCOVER_DIR: Path
    REFERENCE_DEM_DIR: Path

    # Location of scripts called by this project
    SETSM_POSTPROCESSING_PYTHON_DIR: Path
    SETSM_POSTPROCESSING_MATLAB_DIR: Path

    # Working directories
    WORKING_ZONES_DIR: Path

    # Configure automatic discovery and loading of settings from env file
    # The prefix must be included in the entries in the env file. Example:
    #   .env entry: EARTHDEM_MOSAIC_MATLAB_SOURCE_DIR=/some/path
    #   will set:   Settings.MATLAB_SOURCE_DIR
    model_config = SettingsConfigDict(env_file=".env", env_prefix="EARTHDEM_MOSAIC_")
