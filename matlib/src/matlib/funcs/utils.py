from dataclasses import asdict
from typing import Any, Iterable

from matlib.types import MatlabFunctionOptions


def get_activated_flags(options: MatlabFunctionOptions) -> set[str]:
    all_flags = {k: v for k, v in asdict(options).items() if k in options.flag_keys()}
    activated = {str(k) for k, v in all_flags.items() if v}
    return activated


def get_activated_options(options: MatlabFunctionOptions) -> dict[str, Any]:
    all_options = {
        k: v for k, v in asdict(options).items() if k not in options.flag_keys()
    }
    activated = {k: v for k, v in all_options.items() if v is not None}
    return activated


def single_quote(s: str) -> str:
    return f"'{s}'"


def to_str(
    func_name: str, args: Iterable[str], options: dict, flags: Iterable[str]
) -> str:
    quoted_args = [single_quote(arg) for arg in args] if args else []

    quoted_options = []
    for k, v in options.items():
        if isinstance(v, int) or isinstance(v, float):
            # Do not wrap numeric values in single quotes
            quoted_options.append(f"{single_quote(k)},{v}")
        else:
            quoted_options.append(f"{single_quote(k)},{single_quote(v)}")

    quoted_flags = [single_quote(flag) for flag in flags]

    params = ", ".join([*quoted_args, *quoted_options, *quoted_flags])

    return f"{func_name}({params});"
