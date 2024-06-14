import click

from earthdem_mosaic.config import Settings


@click.command()
@click.pass_obj
def show_settings(settings: Settings) -> None:
    """Diplay settings values"""
    click.echo(settings.model_dump_json(indent=2))
