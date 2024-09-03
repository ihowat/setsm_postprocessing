import re
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from earthdem_mosaic.config import Settings

UTM_ZONE_REGEX = re.compile(r"utm\d\d[n|s]")


def replace_in_path(p: Path, old: str, new: str) -> Path:
    return Path(str(p).replace(old, new))


class WorkingSubdirectory(StrEnum):
    Matfiles = "00-matfiles"
    CoregDebug = "10-coregistration-debug"
    NoSlopeFilter = "20-no-slope-filter"
    YesSlopeFilter = "30-yes-slope-filter"


class MatfileSuffix(StrEnum):
    Base = ".mat"
    Reg = "_reg.mat"
    RegFill = "_reg_fill.mat"
    RegFillMerge = "_reg_fill_merge.mat"


@dataclass
class Supertile:
    name: str
    _working_zones_dir: Path
    _landcover_dir: Path
    _reference_dem_dir: Path

    @classmethod
    def from_settings(cls, settings: Settings, name: str):
        return cls(
            name=name,
            _working_zones_dir=settings.WORKING_ZONES_DIR,
            _landcover_dir=settings.LANDCOVER_DIR,
            _reference_dem_dir=settings.REFERENCE_DEM_DIR,
        )

    @property
    def utm_zone(self) -> str:
        return UTM_ZONE_REGEX.search(self.name).group(0)

    @property
    def ref_landcover(self) -> Path:
        filename = f"{self.name}_10m_esa_worldcover_2021.tif"
        return self._landcover_dir / self.utm_zone / filename

    @property
    def ref_dem(self) -> Path:
        filename = f"{self.name}_10m_cop30_wgs84.tif"
        return self._reference_dem_dir / self.utm_zone / filename

    @property
    def skipreg_shp(self) -> Path | None:
        filename = f"{self.name}_skipreg.shp"
        return self._working_zones_dir / self.utm_zone / filename

    def supertile_dir(self, subdir: WorkingSubdirectory) -> Path:
        return self._working_zones_dir / self.utm_zone / f"{subdir}" / self.name

    @property
    def processing_logs_dir(self) -> Path:
        return self._working_zones_dir / self.utm_zone / "processing_logs"

    def find_base_matfiles(self) -> list[Path]:
        filenames = [
            f"{self.name}_1_1_2m.mat",
            f"{self.name}_1_2_2m.mat",
            f"{self.name}_2_1_2m.mat",
            f"{self.name}_2_2_2m.mat",
        ]
        paths = [
            self.supertile_dir(WorkingSubdirectory.Matfiles) / name
            for name in filenames
        ]
        return [p for p in paths if p.exists()]

    def derive_matfiles_from_base(
        self, subdir: WorkingSubdirectory, suffix: MatfileSuffix
    ) -> list[Path]:
        matfiles = self.find_base_matfiles()
        matfiles = [
            replace_in_path(p=mat, old=f"{WorkingSubdirectory.Matfiles}", new=subdir)
            for mat in matfiles
        ]
        return [
            replace_in_path(p=mat, old=f"{MatfileSuffix.Base}", new=suffix)
            for mat in matfiles
        ]
