from dataclasses import dataclass
from pathlib import Path


@dataclass
class AddPath:
    path: str | Path

    def build(self) -> str:
        return f"addpath('{self.path}');"
