from __future__ import annotations

import typer
from rich.console import Console

from panscape import __version__
from panscape.commands import build, map, species, strain

console = Console()

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


@app.callback()
def main(
    version: bool = typer.Option(False, "--version", help="Show PanScape version and exit."),
) -> None:
    if version:
        console.print(f"PanScape {__version__}")
        raise typer.Exit()
