from datetime import datetime

import parsl
from parsl import File, bash_app
from parsl.channels import LocalChannel
from parsl.dataflow.futures import AppFuture
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SrunLauncher
from parsl.providers import SlurmProvider

from earthdem_mosaic.config import Settings
from earthdem_mosaic.pipeline.matlab import (
    AddPath,
    BatchRegisterTilesToCOP30,
    build_matlab_command,
)
from earthdem_mosaic.pipeline.tile import MatfileSuffix, Supertile, WorkingSubdirectory

STAGE_HTEX = "coregister_htex"


def coregister_htex(
    nodes_per_job: int = 1,
    cpus_per_task: int = 4,
    mem_per_task: int = 200,
    job_walltime: str = "48:00:00",
    partition: str = "big_mem",
) -> HighThroughputExecutor:
    return HighThroughputExecutor(
        label=STAGE_HTEX,
        cores_per_worker=cpus_per_task,
        mem_per_worker=mem_per_task,
        provider=SlurmProvider(
            partition=partition,
            channel=LocalChannel(),
            nodes_per_block=nodes_per_job,
            init_blocks=1,
            min_blocks=1,
            max_blocks=100,
            walltime=job_walltime,
            exclusive=False,
            launcher=SrunLauncher(),
        ),
    )


def coregister_inputs(supertile: Supertile) -> list[File]:
    return [File(mat) for mat in supertile.find_base_matfiles()]


def coregister_outputs(supertile: Supertile) -> list[File]:
    return [
        File(mat)
        for mat in supertile.derive_matfiles_from_base(
            subdir=WorkingSubdirectory.Matfiles, suffix=MatfileSuffix.Reg
        )
    ]


def coregister_log(supertile: Supertile) -> str:
    timestamp = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
    filename = f"reg2m_cop30_{supertile.name}.log.{timestamp}"
    return str(supertile.processing_logs_dir / filename)


@bash_app(executors=[STAGE_HTEX])
def coregister_bash_app(
    settings: Settings,
    supertile: Supertile,
    inputs: list[File],
    outputs: list[File],
    stdout: str,
    stderr: str,
) -> AppFuture:
    funcs = [
        AddPath(settings.SETSM_POSTPROCESSING_MATLAB_DIR),
        AddPath(settings.SETSM_POSTPROCESSING_PYTHON_DIR),
        BatchRegisterTilesToCOP30(
            supertile_dir=supertile.supertile_dir(subdir=WorkingSubdirectory.Matfiles),
            ref_dem=supertile.ref_dem,
            ref_landcover=supertile.ref_landcover,
            skipreg_shp=supertile.skipreg_shp,
        ),
    ]
    cmd = build_matlab_command(funcs)
    return " ".join(cmd)
