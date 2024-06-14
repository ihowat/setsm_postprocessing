import click

from earthdem_mosaic.config import Settings
from earthdem_mosaic.slurm.configurable_jobscript import ConfigurableJobscript
from earthdem_mosaic.utm_zone import UtmZone


def wrap_with_try_catch(cmd: str) -> str:
    """Wraps a matlab command in a try-catch block and prepares it for shell execution."""
    return f"""matlab -nojvm -nodisplay -nosplash -r "try; {cmd}; catch e; disp(getReport(e)); exit(1); end; exit(0);" """


@click.command(short_help="Create matfile tile neighbor index")
@click.option("--slurm", is_flag=True, help="Submit jobs to slurm cluster")
@click.option(
    "--show-command",
    is_flag=True,
    help="Print the generated command without executing",
)
@click.option(
    "--show-jobscript",
    is_flag=True,
    help="Print the generated jobscript without executing",
)
@click.argument("utm_zone", nargs=1, type=UtmZone)
@click.pass_obj
def create_neighbor_index(
    settings: Settings,
    slurm: bool,
    show_command: bool,
    show_jobscript: bool,
    utm_zone: UtmZone,
) -> None:
    tiledir = settings.WORKING_ZONES_DIR / str(utm_zone) / "00-matfiles"
    tile_index = (
        settings.WORKING_ZONES_DIR
        / str(utm_zone)
        / "tile_index_files/tileNeighborIndex_2m.mat"
    )
    priority_suffix = "_reg_fill_merge.mat"
    # To prevent the backfill with files other than those matching the priority suffix,
    # pass a secondary suffix that won't match any files.
    secondary_suffix = "_match_nothing.mat"

    cmd = " ".join(
        [
            f"addpath('{settings.SETSM_POSTPROCESSING_MATLAB_DIR}');",
            "tileNeighborIndex(",
            f"'{tiledir}',",
            "'org','pgc',",
            "'resolution','2m',",
            f"'outfile','{tile_index}',",
            f"'priority_suffix','{priority_suffix}',",
            f"'secondary_suffix','{secondary_suffix}'",
            ");",
        ]
    )

    cmd = wrap_with_try_catch(cmd)

    if show_command:
        print(cmd)
        return

    # Configure jobscript for slurm execution
    jobscript_dir = tiledir / "jobscripts"
    log_dir = tiledir / "logs"

    prefix = "create_neighbor_index"
    job_name = f"{prefix}_{utm_zone}"

    jobscript = ConfigurableJobscript(
        jobscript_name_prefix=prefix,
        jobscript_dir=jobscript_dir,
        sbatch_options=[
            "--nodes 1",
            "--ntasks 2",
            "--mem=10G",
            "--time 40:00:00",
            "--partition=batch",
            f"--job-name {job_name}",
            f"--output {log_dir}/%x.o%j",
            f"--error {log_dir}/%x.o%j",
        ],
        conda_env_name="earthdem-mosaic",
        command_to_run=cmd,
    )

    if show_jobscript:
        print(jobscript.content())
        return

    if slurm:
        print(f"Writing jobscript to: {jobscript.path}")
        jobscript.write()
        print(f"Submitting {jobscript.path.name} to queue")
        jobscript.submit()
