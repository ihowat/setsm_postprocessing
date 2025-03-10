from matlib.types import MatlabFunction


def build_matlab_script(funcs: list[MatlabFunction]) -> str:
    built = [func.build() for func in funcs]
    joined = " ".join(built)
    return f"try; {joined}; catch e; disp(getReport(e)); exit(1); end; exit(0);"


def build_bash_command(funcs: list[MatlabFunction]) -> list[str]:
    script = build_matlab_script(funcs)
    return ["matlab", "-nojvm", "-nodisplay", "-r", f'"{script}"']
