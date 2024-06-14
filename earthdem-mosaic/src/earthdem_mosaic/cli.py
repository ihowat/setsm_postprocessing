import click

from earthdem_mosaic.commands import COMMANDS
from earthdem_mosaic.config import Settings


@click.group
@click.pass_context
def cli(ctx: click.Context) -> None:
    """A CLI for executing EarthDEM mosaic production steps"""
    ctx.obj = Settings()


for cmd in COMMANDS:
    cli.add_command(cmd)
