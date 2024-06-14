import click

from earthdem_mosaic.config import Settings
from earthdem_mosaic.utm_zone import UtmZone


@click.command()
@click.option("-v", "--verbose", is_flag=True)
@click.argument("utm_zone", nargs=1, type=UtmZone)
@click.pass_obj
def create_working_dirs(settings: Settings, verbose: bool, utm_zone: UtmZone) -> None:
    """Create working directories for a UTM zone"""
    zone_dir = settings.WORKING_ZONES_DIR / str(utm_zone)
    stage_dirs = [
        zone_dir / "00-matfiles",
        zone_dir / "10-coregistration-debug",
        zone_dir / "20-no-slope-filter",
        zone_dir / "30-yes-slope-filter",
    ]

    click.echo(f"Creating directory: {zone_dir}") if verbose else None
    zone_dir.mkdir(exist_ok=True)

    # create the stage directory itself and a 'logs' subdirectory
    for stage_dir in stage_dirs:
        click.echo(f"Creating directory: {stage_dir}") if verbose else None
        click.echo(f"Creating directory: {stage_dir / 'logs'}") if verbose else None

        stage_dir.mkdir(exist_ok=True)
        (stage_dir / "logs").mkdir(exist_ok=True)

    click.echo(f"Working directories created at: {zone_dir}")
