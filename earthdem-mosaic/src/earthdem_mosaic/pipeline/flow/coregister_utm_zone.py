from pathlib import Path

import parsl
from parsl import File
from parsl.config import Config

from earthdem_mosaic.config import Settings
from earthdem_mosaic.pipeline.stage.coregister import (
    coregister_bash_app,
    coregister_htex,
    coregister_inputs,
    coregister_outputs,
)
from earthdem_mosaic.pipeline.tile import Supertile

StageResult = tuple[Supertile, list[File]]


def main():
    settings = Settings()
    # Override WORKING_ZONES_DIR with testing location
    settings.WORKING_ZONES_DIR = Path("~/scratch/earthdem_parsl_testing")

    supertile_names = [
        "utm15s_66_07",
        "utm15s_67_07",
        "utm15s_66_06",
        "utm15s_67_06",
        "utm15s_67_05",
    ]
    supertiles = [Supertile.from_settings(settings, name) for name in supertile_names]

    config = Config(
        run_dir=str(settings.WORKING_ZONES_DIR / "utm15s" / "runinfo"),
        executors=[coregister_htex()],
    )

    with parsl.load(config):
        coreg_results: list[StageResult] = []
        for supertile in supertiles:
            future = coregister_bash_app(
                settings=settings,
                supertile=supertile,
                inputs=coregister_inputs(supertile),
                outputs=coregister_outputs(supertile),
                stdout=parsl.AUTO_LOGNAME,
                stderr=parsl.AUTO_LOGNAME,
            )
            outputs = [output.result() for output in future.outputs]
            coreg_results.append((supertile, outputs))

        for supertile, outputs in coreg_results:
            print(f"{supertile.name}:")
            for output in outputs:
                print(f"\t{output}")


if __name__ == "__main__":
    main()
