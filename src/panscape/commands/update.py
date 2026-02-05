from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from panscape.exceptions import PanScapeError
from panscape.logging import configure_logging, get_logger
from panscape.utils.subprocess import CommandExecutionError, run_command

app = typer.Typer(help="Update PanScape by reinstalling from the GitHub repository.")
console = Console()


def run_update(
    *,
    log_file: Path | None,
    verbose: bool,
    quiet: bool,
    dry_run: bool,
) -> int:
    try:
        configure_logging(verbose=verbose, quiet=quiet, log_file=log_file)
        logger = get_logger("panscape.update")

        commands = [
            ["pip", "uninstall", "panscape", "-y"],
            ["pip", "install", "git+https://github.com/anttonalberdi/panscape.git"],
        ]

        for command in commands:
            logger.info("Running: %s", " ".join(command))
            result = run_command(command, dry_run=dry_run, logger=logger)
            if result.stdout.strip():
                logger.info(result.stdout.strip())
            if result.stderr.strip():
                logger.info(result.stderr.strip())

        if dry_run:
            logger.info("Dry-run requested; no changes were made.")
        else:
            logger.info("PanScape update completed.")

        return 0

    except CommandExecutionError as exc:
        console.print(f"[red]Update failed:[/red] {exc}")
        return 1
    except PanScapeError as exc:
        console.print(f"[red]Error:[/red] {exc}")
        return exc.exit_code
    except Exception as exc:  # pragma: no cover - defensive catch-all
        get_logger("panscape.update").exception("Unhandled update error")
        console.print(f"[red]Unexpected error:[/red] {exc}")
        return 1


@app.callback(invoke_without_command=True)
def update_callback(
    ctx: typer.Context,
    dry_run: bool = typer.Option(False, "--dry-run", help="Print commands without executing."),
    log_file: Path | None = typer.Option(None, "--log-file", help="Write JSON logs to this file."),
    verbose: bool = typer.Option(False, "--verbose", help="Enable verbose logging."),
    quiet: bool = typer.Option(False, "--quiet", help="Only show errors."),
) -> None:
    if ctx.invoked_subcommand is not None:
        return

    exit_code = run_update(
        log_file=log_file,
        verbose=verbose,
        quiet=quiet,
        dry_run=dry_run,
    )
    raise typer.Exit(exit_code)
