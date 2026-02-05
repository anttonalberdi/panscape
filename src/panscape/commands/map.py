from __future__ import annotations

import csv
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import typer
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

from panscape.config import MapConfig, merge_command_config
from panscape.exceptions import PanScapeError, PanScapeUsageError
from panscape.logging import configure_logging, get_logger
from panscape.manifest import create_run_manifest, finalize_manifest, write_manifest
from panscape.paths import OutputLayout, create_output_layout, map_sample_dir
from panscape.runners.mapper import MapperRunner
from panscape.utils.io import ensure_dir, write_json, write_text
from panscape.utils.subprocess import shell_join
from panscape.utils.validation import SampleRecord, validate_samples_manifest

app = typer.Typer(help="Map reads to species references and produce sample/species BAM outputs.")
console = Console()


def _print_plan(step_plan: list[str]) -> None:
    console.print("[bold]Map step plan[/bold]")
    for idx, step in enumerate(step_plan, start=1):
        console.print(f"  {idx}. {step}")


def _read_species_from_clusters(species_clusters_tsv: Path) -> list[str]:
    species_ids: set[str] = set()
    if not species_clusters_tsv.exists():
        return []

    with species_clusters_tsv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            species_id = (row.get("species_id") or "").strip()
            if species_id:
                species_ids.add(species_id)

    return sorted(species_ids)


def _discover_species_ids(cfg: MapConfig) -> list[str]:
    candidates: list[Path] = []

    if cfg.references_dir is not None:
        if (cfg.references_dir / "pangenomes").is_dir():
            candidates.append(cfg.references_dir / "pangenomes")
        elif cfg.references_dir.is_dir():
            candidates.append(cfg.references_dir)

    candidates.append(cfg.outdir / "build" / "pangenomes")

    for candidate in candidates:
        if candidate.is_dir():
            species_ids = sorted(path.name for path in candidate.iterdir() if path.is_dir())
            if species_ids:
                return species_ids

    return _read_species_from_clusters(cfg.outdir / "build" / "species_clusters.tsv")


def _write_bam_placeholder(
    *,
    layout: OutputLayout,
    sample: SampleRecord,
    species_id: str,
    mapper: MapperRunner,
    threads: int,
    force: bool,
) -> Path:
    sample_dir = ensure_dir(map_sample_dir(layout, sample.sample_id))
    bam_path = sample_dir / f"{species_id}.bam"

    reference = layout.build_dir / "pangenomes" / species_id / "genes.fna"
    command = mapper.command(
        mapper.plan_map_args(
            reference_fasta=reference,
            r1=sample.r1,
            r2=sample.r2,
            output_bam=bam_path,
            threads=threads,
        )
    )

    content = "\n".join(
        [
            "# PanScape placeholder BAM (not a real BAM binary)",
            f"# sample_id: {sample.sample_id}",
            f"# species_id: {species_id}",
            f"# mapper: {mapper.executable}",
            f"# planned_command: {shell_join(command)}",
            "@HD\tVN:1.6\tSO:unknown",
            f"@RG\tID:{sample.sample_id}\tSM:{sample.sample_id}",
            "",
        ]
    )
    return write_text(bam_path, content, force=force)


def run_map(
    *,
    config_path: Path | None,
    samples_tsv: Path | None,
    references_dir: Path | None,
    mapper: str | None,
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
            section="map",
            model_cls=MapConfig,
            cli_overrides={
                "samples_tsv": samples_tsv,
                "references_dir": references_dir,
                "mapper": mapper,
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
        logger = get_logger("panscape.map")

        if cfg.samples_tsv is None:
            raise PanScapeUsageError(
                "Missing samples manifest. Provide --samples-tsv or set map.samples_tsv in config."
            )

        samples = validate_samples_manifest(cfg.samples_tsv)
        layout = create_output_layout(cfg.outdir)

        species_ids = _discover_species_ids(cfg)
        if not species_ids:
            logger.warning(
                "No species references detected. Falling back to a default placeholder species_id."
            )
            species_ids = ["sp0001"]

        step_plan = [
            f"Validate samples manifest ({cfg.samples_tsv})",
            f"Discover species references ({len(species_ids)} species)",
            f"Create map outputs under {layout.map_dir}",
            f"Generate placeholder BAMs with mapper adapter ({cfg.mapper})",
        ]

        manifest = create_run_manifest(
            command="map",
            argv=sys.argv,
            outdir=layout.root,
            dry_run=cfg.dry_run,
            threads=cfg.threads,
            config_path=config_path,
            input_paths=[cfg.samples_tsv],
            planned_steps=step_plan,
        )
        write_manifest(layout.root, manifest)

        _print_plan(step_plan)
        if cfg.dry_run:
            logger.info("Dry-run requested; stopping before writing map outputs.")
            finalize_manifest(manifest, status="dry-run", output_paths=[])
            write_manifest(layout.root, manifest)
            return 0

        mapper_runner = MapperRunner(cfg.mapper)
        outputs: list[Path] = []

        tasks = [(sample, species_id) for sample in samples for species_id in species_ids]
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TimeElapsedColumn(),
            console=console,
        ) as progress:
            task_id = progress.add_task("Creating BAM placeholders", total=len(tasks))
            with ThreadPoolExecutor(max_workers=cfg.threads) as pool:
                futures = [
                    pool.submit(
                        _write_bam_placeholder,
                        layout=layout,
                        sample=sample,
                        species_id=species_id,
                        mapper=mapper_runner,
                        threads=cfg.threads,
                        force=cfg.force,
                    )
                    for sample, species_id in tasks
                ]
                for future in as_completed(futures):
                    outputs.append(future.result())
                    progress.advance(task_id)

        for sample in samples:
            sample_dir = ensure_dir(map_sample_dir(layout, sample.sample_id))
            qc_payload = {
                "sample_id": sample.sample_id,
                "mapper": cfg.mapper,
                "species_count": len(species_ids),
                "notes": "placeholder metrics",
            }
            outputs.append(write_json(sample_dir / "qc.json", qc_payload, force=cfg.force))

        finalize_manifest(manifest, status="completed", output_paths=outputs)
        write_manifest(layout.root, manifest)
        logger.info("Map scaffold completed for %d samples.", len(samples))
        return 0

    except PanScapeError as exc:
        console.print(f"[red]Error:[/red] {exc}")
        return exc.exit_code
    except Exception as exc:  # pragma: no cover - defensive catch-all
        get_logger("panscape.map").exception("Unhandled map error")
        console.print(f"[red]Unexpected error:[/red] {exc}")
        return 1


@app.callback(invoke_without_command=True)
def map_callback(
    ctx: typer.Context,
    config: Path | None = typer.Option(None, "--config", help="YAML config file."),
    samples_tsv: Path | None = typer.Option(
        None,
        "--samples-tsv",
        help="Sample manifest TSV with columns: sample_id, r1, r2(optional).",
    ),
    references_dir: Path | None = typer.Option(
        None,
        "--references-dir",
        help="Directory containing pangenome references (optional).",
    ),
    mapper: str | None = typer.Option(
        None,
        "--mapper",
        help="Mapper executable name or path (e.g., minimap2, bwa).",
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

    exit_code = run_map(
        config_path=config,
        samples_tsv=samples_tsv,
        references_dir=references_dir,
        mapper=mapper,
        outdir=outdir,
        threads=threads,
        dry_run=dry_run,
        force=force,
        log_file=log_file,
        verbose=verbose,
        quiet=quiet,
    )
    raise typer.Exit(exit_code)
