from dataclasses import dataclass
from pathlib import Path
from typing import Optional


def parse_tiles_input(tiles: str) -> list[str]:
    # Text file with newline separated tile names
    if tiles.strip().endswith(".txt"):
        tiles_file = Path(tiles)
        if not tiles_file.exists():
            raise FileNotFoundError(f"{tiles} does not exist")
        with open(tiles_file, "r") as f:
            return sorted([line.strip() for line in f])

    # String with comma separated tile names
    if "," in tiles:
        return sorted(tiles.split(","))

    # Single tile name
    return [tiles]


def make_dirs_if_not_exist(paths: list[Path]) -> None:
    for p in paths:
        p.mkdir(exist_ok=True, parents=True)


@dataclass(kw_only=True)
class Command:
    program: Optional[str] = None
    script: Optional[str] = None
    args: Optional[list[str]] = None
    options: Optional[dict[str, str]] = None
    flags: Optional[list[str]] = None
    postfix: Optional[str] = None

    def to_pretty_str(self) -> str:
        cmd = []
        if self.program:
            cmd = [self.program]
        if self.script:
            cmd = [*cmd, self.script]
        if self.args:
            cmd = [*cmd, *self.args]
        if self.flags:
            cmd = [*cmd, *self.flags]
        if self.options:
            cmd = [*cmd, *[f"{k} {v}" for k, v in self.options.items()]]
        if self.postfix:
            cmd = [*cmd, self.postfix]

        return "\n\t".join(cmd)

    def to_list(self) -> list[str]:
        cmd = []
        if self.program:
            cmd = [self.program]
        if self.script:
            cmd = [*cmd, self.script]
        if self.args:
            cmd = [*cmd, *self.args]
        if self.flags:
            cmd = [*cmd, *self.flags]
        if self.options:
            cmd = [*cmd, *_unpack_options(self.options)]
        if self.postfix:
            cmd = [*cmd, self.postfix]

        return cmd


def _unpack_options(d: dict[str, str]) -> list[str]:
    options = []
    for key, value in d.items():
        if value:
            options.append(key)
            options.append(value)
    return options
