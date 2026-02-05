from __future__ import annotations

import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import typer
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

from panscape.config import BuildConfig, merge_command_config
from panscape.exceptions import PanScapeError, PanScapeUsageError
from panscape.logging import configure_logging, get_logger
from panscape.manifest import create_run_manifest, finalize_manifest, write_manifest
from panscape.paths import OutputLayout, create_output_layout, pangenome_dir
from panscape.utils.io import ensure_dir, write_text, write_tsv
from panscape.utils.validation import GenomeRecord, validate_genomes_manifest

app = typer.Typer(help="Create species clusters, pangenome refs, and mapper indices.")
console = Console()


def _assign_species_ids(genomes: list[GenomeRecord]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for idx, genome in enumerate(genomes, start=1):
        mapping[genome.genome_id] = f"sp{idx:04d}"
    return mapping


def _print_plan(step_plan: list[str]) -> None:
    console.print("[bold]Build step plan[/bold]")
    for idx, step in enumerate(step_plan, start=1):
        console.print(f"  {idx}. {step}")


def _create_species_placeholder_outputs(
    *, layout: OutputLayout, species_id: str, force: bool
) -> list[Path]:
    species_dir = pangenome_dir(layout, species_id)
    index_dir = ensure_dir(species_dir / "index")

    outputs: list[Path] = []
    outputs.append(
        write_text(
            species_dir / "genes.fna",
            f">{species_id}_gene_0001\nNNNNNNNNNN\n",
            force=force,
        )
    )
    outputs.append(
        write_text(
            species_dir / "genes.faa",
            f">{species_id}_gene_0001\nMXXXXXXXXX\n",
            force=force,
        )
    )
    outputs.append(
        write_tsv(
            species_dir / "gene_families.tsv",
            ["family_id", "gene_id"],
            [[f"{species_id}_fam_1", f"{species_id}_gene_0001"]],
            force=force,
        )
    )
    outputs.append(
        write_text(
            index_dir / "README.txt",
            "TODO: Populate mapper-specific index files here.\n",
            force=force,
        )
    )
    return outputs


def run_build(
    *,
    config_path: Path | None,
    genomes_tsv: Path | None,
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
            section="build",
            model_cls=BuildConfig,
            cli_overrides={
                "genomes_tsv": genomes_tsv,
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
        logger = get_logger("panscape.build")

        if cfg.genomes_tsv is None:
            raise PanScapeUsageError(
                "Missing genomes manifest. Provide --genomes-tsv or set build.genomes_tsv in config."
            )

        genomes = validate_genomes_manifest(cfg.genomes_tsv)
        layout = create_output_layout(cfg.outdir)

        species_by_genome = _assign_species_ids(genomes)
        species_ids = sorted(set(species_by_genome.values()))

        step_plan = [
            f"Validate genomes manifest ({cfg.genomes_tsv})",
            f"Create output layout under {layout.root}",
            "Assign placeholder species clusters and backbone genomes",
            "Create pangenome placeholder files per species",
        ]

        manifest = create_run_manifest(
            command="build",
            argv=sys.argv,
            outdir=layout.root,
            dry_run=cfg.dry_run,
            threads=cfg.threads,
            config_path=config_path,
            input_paths=[cfg.genomes_tsv],
            planned_steps=step_plan,
        )
        write_manifest(layout.root, manifest)

        _print_plan(step_plan)
        if cfg.dry_run:
            logger.info("Dry-run requested; stopping before writing build outputs.")
            finalize_manifest(manifest, status="dry-run", output_paths=[])
            write_manifest(layout.root, manifest)
            return 0

        outputs: list[Path] = []
        outputs.append(
            write_tsv(
                layout.build_dir / "species_clusters.tsv",
                ["genome_id", "species_id"],
                [[record.genome_id, species_by_genome[record.genome_id]] for record in genomes],
                force=cfg.force,
            )
        )

        backbones_dir = ensure_dir(layout.build_dir / "backbones")
        for species_id in species_ids:
            outputs.append(
                write_text(
                    backbones_dir / f"{species_id}.fna",
                    f">{species_id}_backbone\nNNNNNNNNNNNN\n",
                    force=cfg.force,
                )
            )

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TimeElapsedColumn(),
            console=console,
        ) as progress:
            task_id = progress.add_task("Creating pangenome placeholders", total=len(species_ids))
            with ThreadPoolExecutor(max_workers=cfg.threads) as pool:
                futures = [
                    pool.submit(
                        _create_species_placeholder_outputs,
                        layout=layout,
                        species_id=species_id,
                        force=cfg.force,
                    )
                    for species_id in species_ids
                ]
                for future in as_completed(futures):
                    outputs.extend(future.result())
                    progress.advance(task_id)

        finalize_manifest(manifest, status="completed", output_paths=outputs)
        write_manifest(layout.root, manifest)
        logger.info("Build scaffold completed with %d species.", len(species_ids))
        return 0

    except PanScapeError as exc:
        console.print(f"[red]Error:[/red] {exc}")
        return exc.exit_code
    except Exception as exc:  # pragma: no cover - defensive catch-all
        get_logger("panscape.build").exception("Unhandled build error")
        console.print(f"[red]Unexpected error:[/red] {exc}")
        return 1


@app.callback(invoke_without_command=True)
def build_callback(
    ctx: typer.Context,
    config: Path | None = typer.Option(None, "--config", help="YAML config file."),
    genomes_tsv: Path | None = typer.Option(
        None,
        "--genomes-tsv",
        help="Genome manifest TSV with columns: genome_id, fasta_path, ...",
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

    exit_code = run_build(
        config_path=config,
        genomes_tsv=genomes_tsv,
        outdir=outdir,
        threads=threads,
        dry_run=dry_run,
        force=force,
        log_file=log_file,
        verbose=verbose,
        quiet=quiet,
    )
    raise typer.Exit(exit_code)
