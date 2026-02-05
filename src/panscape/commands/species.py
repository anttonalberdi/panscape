from __future__ import annotations

import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import typer
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

from panscape.config import SpeciesConfig, merge_command_config
from panscape.exceptions import PanScapeError, PanScapeUsageError
from panscape.logging import configure_logging, get_logger
from panscape.manifest import create_run_manifest, finalize_manifest, write_manifest
from panscape.paths import OutputLayout, create_output_layout, species_output_dir
from panscape.utils.io import ensure_dir, write_json, write_tsv

app = typer.Typer(help="Derive species-level coverage and QC artifacts from map outputs.")
console = Console()


def _print_plan(step_plan: list[str]) -> None:
    console.print("[bold]Species step plan[/bold]")
    for idx, step in enumerate(step_plan, start=1):
        console.print(f"  {idx}. {step}")


def _discover_species_from_map(map_dir: Path) -> list[str]:
    species_ids: set[str] = set()
    if not map_dir.is_dir():
        return []

    for sample_dir in map_dir.iterdir():
        if not sample_dir.is_dir():
            continue
        for bam_path in sample_dir.glob("*.bam"):
            species_ids.add(bam_path.stem)

    return sorted(species_ids)


def _collect_species_rows(map_dir: Path, species_id: str) -> list[list[str]]:
    rows: list[list[str]] = []
    for sample_dir in sorted(path for path in map_dir.iterdir() if path.is_dir()):
        bam_path = sample_dir / f"{species_id}.bam"
        if bam_path.exists():
            rows.append([sample_dir.name, "0.0", "0.0"])

    if not rows:
        rows.append(["placeholder_sample", "0.0", "0.0"])

    return rows


def _write_species_placeholders(
    *,
    layout: OutputLayout,
    map_dir: Path,
    species_id: str,
    force: bool,
) -> list[Path]:
    species_dir = ensure_dir(species_output_dir(layout, species_id))

    coverage_rows = _collect_species_rows(map_dir, species_id)
    coverage_path = write_tsv(
        species_dir / "coverage_matrix.tsv",
        ["sample_id", "mean_coverage", "breadth"],
        coverage_rows,
        force=force,
    )
    qc_path = write_json(
        species_dir / "qc_summary.json",
        {
            "species_id": species_id,
            "sample_count": len(coverage_rows),
            "notes": "placeholder species-level QC",
        },
        force=force,
    )
    return [coverage_path, qc_path]


def run_species(
    *,
    config_path: Path | None,
    map_dir: Path | None,
    species_ids: list[str] | None,
    outdir: Path | None,
    threads: int | None,
    dry_run: bool | None,
    force: bool | None,
    log_file: Path | None,
    verbose: bool | None,
    quiet: bool | None,
) -> int:
    try:
        cfg = merge_command_config(
            config_path=config_path,
            section="species",
            model_cls=SpeciesConfig,
            cli_overrides={
                "map_dir": map_dir,
                "species_ids": species_ids,
                "outdir": outdir,
                "threads": threads,
                "dry_run": dry_run,
                "force": force,
                "log_file": log_file,
                "verbose": verbose,
                "quiet": quiet,
            },
        )
        configure_logging(verbose=cfg.verbose, quiet=cfg.quiet, log_file=cfg.log_file)
        logger = get_logger("panscape.species")

        resolved_map_dir = cfg.map_dir or (cfg.outdir / "map")
        if not resolved_map_dir.exists() or not resolved_map_dir.is_dir():
            raise PanScapeUsageError(f"Map directory does not exist: {resolved_map_dir}")

        selected_species = cfg.species_ids or _discover_species_from_map(resolved_map_dir)
        if not selected_species:
            raise PanScapeUsageError(
                "No species IDs discovered. Provide --species-id or ensure map outputs exist."
            )

        layout = create_output_layout(cfg.outdir)

        step_plan = [
            f"Inspect map outputs under {resolved_map_dir}",
            f"Select species IDs ({len(selected_species)} species)",
            f"Create species artifacts under {layout.species_dir}",
            "Write placeholder coverage matrices and QC summaries",
        ]

        manifest = create_run_manifest(
            command="species",
            argv=sys.argv,
            outdir=layout.root,
            dry_run=cfg.dry_run,
            threads=cfg.threads,
            config_path=config_path,
            input_paths=[resolved_map_dir],
            planned_steps=step_plan,
        )
        write_manifest(layout.root, manifest)

        _print_plan(step_plan)
        if cfg.dry_run:
            logger.info("Dry-run requested; stopping before writing species outputs.")
            finalize_manifest(manifest, status="dry-run", output_paths=[])
            write_manifest(layout.root, manifest)
            return 0

        outputs: list[Path] = []
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TimeElapsedColumn(),
            console=console,
        ) as progress:
            task_id = progress.add_task("Creating species placeholders", total=len(selected_species))
            with ThreadPoolExecutor(max_workers=cfg.threads) as pool:
                futures = [
                    pool.submit(
                        _write_species_placeholders,
                        layout=layout,
                        map_dir=resolved_map_dir,
                        species_id=species_id,
                        force=cfg.force,
                    )
                    for species_id in selected_species
                ]
                for future in as_completed(futures):
                    outputs.extend(future.result())
                    progress.advance(task_id)

        finalize_manifest(manifest, status="completed", output_paths=outputs)
        write_manifest(layout.root, manifest)
        logger.info("Species scaffold completed for %d species.", len(selected_species))
        return 0

    except PanScapeError as exc:
        console.print(f"[red]Error:[/red] {exc}")
        return exc.exit_code
    except Exception as exc:  # pragma: no cover - defensive catch-all
        get_logger("panscape.species").exception("Unhandled species error")
        console.print(f"[red]Unexpected error:[/red] {exc}")
        return 1


@app.callback(invoke_without_command=True)
def species_callback(
    ctx: typer.Context,
    config: Path | None = typer.Option(None, "--config", help="YAML config file."),
    map_dir: Path | None = typer.Option(None, "--map-dir", help="Input map output directory."),
    species_id: list[str] | None = typer.Option(
        None,
        "--species-id",
        help="Species ID to process. Repeat for multiple species.",
    ),
    outdir: Path | None = typer.Option(None, "--outdir", help="Output root directory."),
    threads: int | None = typer.Option(None, "--threads", min=1, help="Worker threads."),
    dry_run: bool | None = typer.Option(None, "--dry-run", help="Plan only, do not write outputs."),
    force: bool | None = typer.Option(
        None,
        "--force",
        help="Overwrite existing placeholder outputs.",
    ),
    log_file: Path | None = typer.Option(None, "--log-file", help="Write JSON logs to this file."),
    verbose: bool | None = typer.Option(None, "--verbose", help="Enable verbose logging."),
    quiet: bool | None = typer.Option(None, "--quiet", help="Only show errors."),
) -> None:
    if ctx.invoked_subcommand is not None:
        return

    exit_code = run_species(
        config_path=config,
        map_dir=map_dir,
        species_ids=species_id,
        outdir=outdir,
        threads=threads,
        dry_run=dry_run,
        force=force,
        log_file=log_file,
        verbose=verbose,
        quiet=quiet,
    )
    raise typer.Exit(exit_code)
