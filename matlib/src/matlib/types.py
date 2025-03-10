from typing import Protocol


class MatlabFunction(Protocol):
    def build(self) -> str: ...


class MatlabFunctionOptions(Protocol):
    def flag_keys(self) -> set[str]: ...
