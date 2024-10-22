import os
from pathlib import Path

from pydantic_settings import BaseSettings, SettingsConfigDict

ENV_FILE_PATH_VARIABLE = "REMA_MOSAIC_ENV_FILE"
ENV_PREFIX = "REMA_MOSAIC_"


def find_env_file() -> Path:
    from_env_variable = os.getenv(ENV_FILE_PATH_VARIABLE)
    if from_env_variable:
        return Path(from_env_variable)
    # Use the .env file at the root of the package
    return Path(__file__).resolve().parents[2] / ".env"


class Settings(BaseSettings):
    SETSM_POSTPROCESSING_PYTHON_DIR: Path
    """Directory of Python batch processing scripts wrapping matlab algorithms"""

    SETSM_POSTPROCESSING_MATLAB_DIR: Path
    """Directory of matlab algorithms"""

    # Configure automatic discovery and loading of settings from env file
    # The prefix must be included in the entries in the env file. Example:
    #   .env entry: REMA_MOSAIC_SETSM_POSTPROCESSING_PYTHON_DIR=/some/path
    #   will set:   Settings.SETSM_POSTPROCESSING_PYTHON_DIR
    model_config = SettingsConfigDict(
        # The .env file is expected to be at the root of the package
        env_file=find_env_file(),
        env_prefix=ENV_PREFIX,
    )
