from __future__ import annotations

import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import typer
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

from panscape.config import StrainConfig, merge_command_config
from panscape.exceptions import PanScapeError, PanScapeUsageError
from panscape.logging import configure_logging, get_logger
from panscape.manifest import create_run_manifest, finalize_manifest, write_manifest
from panscape.paths import OutputLayout, create_output_layout, strain_k_dir
from panscape.utils.io import ensure_dir, write_json, write_tsv
from panscape.utils.validation import validate_optional_file

app = typer.Typer(help="Run strain deconvolution (placeholder) from matrices + BAM/VCF inputs.")
console = Console()


def _print_plan(step_plan: list[str]) -> None:
    console.print("[bold]Strain step plan[/bold]")
    for idx, step in enumerate(step_plan, start=1):
        console.print(f"  {idx}. {step}")


def _discover_species_ids(outdir: Path) -> list[str]:
    species_root = outdir / "species"
    if not species_root.is_dir():
        return []
    return sorted(path.name for path in species_root.iterdir() if path.is_dir())


def _write_k_placeholders(
    *,
    layout: OutputLayout,
    species_id: str,
    k_value: int,
    force: bool,
) -> list[Path]:
    k_dir = ensure_dir(strain_k_dir(layout, species_id, k_value))
    abundance = write_tsv(
        k_dir / "abundance.tsv",
        ["sample_id", "strain_id", "relative_abundance"],
        [["placeholder_sample", f"{species_id}_strain_1", "1.0"]],
        force=force,
    )
    gene_content = write_tsv(
        k_dir / "gene_content.tsv",
        ["strain_id", "gene_id", "present"],
        [[f"{species_id}_strain_1", f"{species_id}_gene_0001", "1"]],
        force=force,
    )
    model_json = write_json(
        k_dir / "model.json",
        {
            "species_id": species_id,
            "k": k_value,
            "status": "placeholder",
            "notes": "TODO: Implement probabilistic strain deconvolution model.",
        },
        force=force,
    )
    return [abundance, gene_content, model_json]


def run_strain(
    *,
    config_path: Path | None,
    species_ids: list[str] | None,
    coverage_matrix: Path | None,
    bam_dir: Path | None,
    vcf_path: Path | None,
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
            section="strain",
            model_cls=StrainConfig,
            cli_overrides={
                "species_ids": species_ids,
                "coverage_matrix": coverage_matrix,
                "bam_dir": bam_dir,
                "vcf_path": vcf_path,
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
        logger = get_logger("panscape.strain")

        validate_optional_file(cfg.coverage_matrix, "coverage matrix")
        validate_optional_file(cfg.vcf_path, "VCF file")
        if cfg.bam_dir is not None and (not cfg.bam_dir.exists() or not cfg.bam_dir.is_dir()):
            raise PanScapeUsageError(f"BAM directory does not exist: {cfg.bam_dir}")

        selected_species = cfg.species_ids or _discover_species_ids(cfg.outdir)
        if not selected_species:
            selected_species = ["sp0001"]
            logger.warning(
                "No species IDs supplied or discovered; using fallback placeholder species_id=%s",
                selected_species[0],
            )

        layout = create_output_layout(cfg.outdir)
        k_values = [2, 3]

        step_plan = [
            "Validate optional matrix/BAM/VCF inputs",
            f"Select species IDs ({len(selected_species)} species)",
            f"Create strain output layout under {layout.strain_dir}",
            f"Write placeholder outputs for K={','.join(str(k) for k in k_values)}",
        ]

        input_paths = [path for path in [cfg.coverage_matrix, cfg.bam_dir, cfg.vcf_path] if path is not None]
        manifest = create_run_manifest(
            command="strain",
            argv=sys.argv,
            outdir=layout.root,
            dry_run=cfg.dry_run,
            threads=cfg.threads,
            config_path=config_path,
            input_paths=input_paths,
            planned_steps=step_plan,
        )
        write_manifest(layout.root, manifest)

        _print_plan(step_plan)
        if cfg.dry_run:
            logger.info("Dry-run requested; stopping before writing strain outputs.")
            finalize_manifest(manifest, status="dry-run", output_paths=[])
            write_manifest(layout.root, manifest)
            return 0

        outputs: list[Path] = []
        tasks = [(species_id, k) for species_id in selected_species for k in k_values]
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TimeElapsedColumn(),
            console=console,
        ) as progress:
            task_id = progress.add_task("Creating strain placeholders", total=len(tasks))
            with ThreadPoolExecutor(max_workers=cfg.threads) as pool:
                futures = [
                    pool.submit(
                        _write_k_placeholders,
                        layout=layout,
                        species_id=species_id,
                        k_value=k_value,
                        force=cfg.force,
                    )
                    for species_id, k_value in tasks
                ]
                for future in as_completed(futures):
                    outputs.extend(future.result())
                    progress.advance(task_id)

        finalize_manifest(manifest, status="completed", output_paths=outputs)
        write_manifest(layout.root, manifest)
        logger.info("Strain scaffold completed for %d species.", len(selected_species))
        return 0

    except PanScapeError as exc:
        console.print(f"[red]Error:[/red] {exc}")
        return exc.exit_code
    except Exception as exc:  # pragma: no cover - defensive catch-all
        get_logger("panscape.strain").exception("Unhandled strain error")
        console.print(f"[red]Unexpected error:[/red] {exc}")
        return 1


@app.callback(invoke_without_command=True)
def strain_callback(
    ctx: typer.Context,
    config: Path | None = typer.Option(None, "--config", help="YAML config file."),
    species_id: list[str] | None = typer.Option(
        None,
        "--species-id",
        help="Species ID to process. Repeat for multiple species.",
    ),
    coverage_matrix: Path | None = typer.Option(
        None,
        "--coverage-matrix",
        help="Coverage matrix TSV (optional).",
    ),
    bam_dir: Path | None = typer.Option(None, "--bam-dir", help="Directory with BAMs (optional)."),
    vcf_path: Path | None = typer.Option(None, "--vcf", help="VCF path (optional)."),
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

    exit_code = run_strain(
        config_path=config,
        species_ids=species_id,
        coverage_matrix=coverage_matrix,
        bam_dir=bam_dir,
        vcf_path=vcf_path,
        outdir=outdir,
        threads=threads,
        dry_run=dry_run,
        force=force,
        log_file=log_file,
        verbose=verbose,
        quiet=quiet,
    )
    raise typer.Exit(exit_code)
