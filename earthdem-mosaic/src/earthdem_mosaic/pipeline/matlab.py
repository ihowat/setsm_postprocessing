from dataclasses import dataclass
from pathlib import Path
from typing import Protocol


class MatlabFunction(Protocol):
    def as_str(self) -> str: ...


def build_matlab_command(funcs: list[MatlabFunction]) -> list[str]:
    joined = " ".join([f.as_str() for f in funcs])
    wrapped = wrap_in_try_catch(joined)
    return ["matlab", "-nojvm", "-nodisplay", "-nosplash", "-r", wrapped]


def wrap_in_try_catch(s: str) -> str:
    return f'"try; {s} catch e; disp(getReport(e)); exit(1); end; exit(0);"'


@dataclass
class Varargin:
    flags: list[str] | None = None
    kwargs: dict[[str], str] | None = None

    def as_str(self, prefix: str = ", ") -> str:
        params = []
        if self.flags:
            for flag in self.flags:
                params.append(f"'{flag}'")
        if self.kwargs:
            for k, v in self.kwargs.items():
                params.append(f"'{k}','{v}'")

        param_str = ", ".join(params)
        if prefix:
            return f"{prefix}{param_str}"
        return param_str

    def add_flag(self, flag: str) -> None:
        self.flags.append(f"'{flag}'")

    def add_kwarg(self, kwarg: dict[[str], str]) -> None:
        self.kwargs.update(kwarg)


@dataclass
class AddPath:
    path: Path

    def as_str(self) -> str:
        return f"addpath('{self.path}');"


@dataclass
class CopyFile:
    src: Path
    dst: Path

    def as_str(self) -> str:
        return f"copyfile('{self.src}', '{self.dst}');"


@dataclass
class BatchRegisterTilesToCOP30:
    supertile_dir: Path
    ref_dem: Path
    ref_landcover: Path
    skipreg_shp: Path | None = None

    def as_str(self) -> str:
        varargin = Varargin(
            flags=["registerBlobs"],
            kwargs={
                "resolution": "2m",
            },
        )
        if self.skipreg_shp:
            varargin.add_kwarg({"registerBlobsSkipregShp": f"{self.skipreg_shp}"})

        return f"batchRegisterTilesToCOP30('{self.supertile_dir}', '{self.ref_dem}', '{self.ref_landcover}'{varargin.as_str()});"
