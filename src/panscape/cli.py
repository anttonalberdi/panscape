from __future__ import annotations

import platform
import sys

import typer
from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from panscape import __version__
from panscape.commands import build, map, species, strain, update

console = Console()
SUBCOMMANDS = ["build", "map", "species", "strain", "update"]

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=(
        "PanScape command-line toolkit for species-level reference construction, mapping, "
        "and strain deconvolution workflows."
    ),
)

app.add_typer(build.app, name="build", help="Build species clusters and pangenome references.")
app.add_typer(map.app, name="map", help="Map reads to species references.")
app.add_typer(species.app, name="species", help="Create species-level mapping artifacts.")
app.add_typer(strain.app, name="strain", help="Run strain deconvolution entrypoint.")
app.add_typer(update.app, name="update", help="Update PanScape from GitHub.")


def _print_startup_intro(command_name: str) -> None:
    banner = Panel(
        f"[bold cyan]PanScape {__version__}[/bold cyan]\n"
        "[white]Species-resolved pangenome toolkit[/white]",
        title="[bold]CLI Start[/bold]",
        border_style="cyan",
        expand=False,
    )
    console.print(banner)

    stats = Table(
        title="[bold]Session Summary[/bold]",
        box=box.SIMPLE_HEAVY,
        show_header=False,
        expand=False,
    )
    stats.add_column("Key", style="bold cyan")
    stats.add_column("Value", style="white")
    stats.add_row("Command", command_name)
    stats.add_row("Subcommands", str(len(SUBCOMMANDS)))
    stats.add_row("Python", sys.version.split()[0])
    stats.add_row("Platform", f"{platform.system()} {platform.release()}")
    console.print(stats)


@app.callback()
def main(
    ctx: typer.Context,
    version: bool = typer.Option(False, "--version", help="Show PanScape version and exit."),
) -> None:
    if version:
        console.print(f"PanScape {__version__}")
        raise typer.Exit()

    if ctx.invoked_subcommand:
        _print_startup_intro(ctx.invoked_subcommand)
