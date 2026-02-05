from __future__ import annotations

import csv
import math
import re
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable

import typer
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

from panscape.build.derep import DerepGene, DerepGroup, dereplicate_family
from panscape.build.fasta import AssemblyStats, FastaRecord, normalize_fasta_file, write_fasta_records
from panscape.build.gene_calling import CalledGene, call_genes, write_gene_outputs
from panscape.build.graph import connected_components
from panscape.build.kmer import mock_mash_distance, unique_kmer_fractions
from panscape.build.scoring import backbone_score
from panscape.config import BuildConfig, merge_command_config
from panscape.exceptions import PanScapeError, PanScapeUsageError
from panscape.logging import configure_logging, get_logger
from panscape.manifest import create_run_manifest, finalize_manifest, write_manifest
from panscape.paths import create_output_layout
from panscape.runners.checkm2 import CheckM2Runner
from panscape.runners.mash import MashRunner
from panscape.runners.mmseqs import MMseqsRunner
from panscape.runners.skani import SkaniRunner
from panscape.utils.io import ensure_dir, write_json, write_text, write_tsv
from panscape.utils.validation import FASTA_SUFFIXES, validate_existing_file

app = typer.Typer(help="Build species clusters, backbones, and species pangenome catalogs.")
console = Console()


@dataclass(frozen=True, slots=True)
class GenomeInput:
    genome_id: str
    fasta_path: Path
    completeness: float | None
    contamination: float | None


@dataclass(frozen=True, slots=True)
class GenomeContext:
    genome: GenomeInput
    normalized_fasta: Path
    normalized_records: list[FastaRecord]
    assembly_stats: AssemblyStats


@dataclass(frozen=True, slots=True)
class GenomeQC:
    completeness: float | None
    contamination: float | None
    qc_source: str
    passes_thresholds: bool | None
    keep_for_clustering: bool


@dataclass(frozen=True, slots=True)
class FamilyCluster:
    family_id: str
    representative_gene_id: str
    member_gene_ids: tuple[str, ...]


def _parse_optional_float(raw_value: str, *, field_name: str, row_number: int) -> float | None:
    value = raw_value.strip()
    if value == "":
        return None
    try:
        return float(value)
    except ValueError as exc:
        raise PanScapeUsageError(
            f"Invalid float value in {field_name} at row {row_number}: {raw_value!r}"
        ) from exc


def _resolve_manifest_path(base_dir: Path, raw_value: str) -> Path:
    path = Path(raw_value).expanduser()
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


def parse_genomes_manifest(genomes_tsv: Path) -> list[GenomeInput]:
    """Read genomes.tsv and validate required columns + referenced FASTAs."""

    validate_existing_file(genomes_tsv, "genomes.tsv")

    with genomes_tsv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise PanScapeUsageError(f"genomes.tsv has no header row: {genomes_tsv}")

        header = [name.strip() for name in reader.fieldnames]
        required = {"genome_id", "fasta_path"}
        missing = required.difference(header)
        if missing:
            raise PanScapeUsageError(
                f"Missing required columns in genomes.tsv ({genomes_tsv}): {', '.join(sorted(missing))}"
            )

        rows = [dict(row) for row in reader]

    if not rows:
        raise PanScapeUsageError(f"genomes.tsv has no data rows: {genomes_tsv}")

    observed: set[str] = set()
    genomes: list[GenomeInput] = []

    for idx, row in enumerate(rows, start=2):
        genome_id = (row.get("genome_id") or "").strip()
        fasta_raw = (row.get("fasta_path") or "").strip()
        completeness_raw = row.get("completeness") or ""
        contamination_raw = row.get("contamination") or ""

        if genome_id == "":
            raise PanScapeUsageError(f"Missing genome_id at row {idx} in {genomes_tsv}")
        if genome_id in observed:
            raise PanScapeUsageError(f"Duplicated genome_id in genomes.tsv: {genome_id}")
        observed.add(genome_id)

        if fasta_raw == "":
            raise PanScapeUsageError(f"Missing fasta_path at row {idx} in {genomes_tsv}")

        fasta_path = _resolve_manifest_path(genomes_tsv.parent, fasta_raw)
        validate_existing_file(fasta_path, f"FASTA for genome {genome_id}")
        lowered = fasta_path.name.lower()
        if not any(lowered.endswith(suffix) for suffix in FASTA_SUFFIXES):
            raise PanScapeUsageError(
                f"Unsupported FASTA extension for genome {genome_id}: {fasta_path}"
            )

        genomes.append(
            GenomeInput(
                genome_id=genome_id,
                fasta_path=fasta_path,
                completeness=_parse_optional_float(
                    completeness_raw,
                    field_name="completeness",
                    row_number=idx,
                ),
                contamination=_parse_optional_float(
                    contamination_raw,
                    field_name="contamination",
                    row_number=idx,
                ),
            )
        )

    return sorted(genomes, key=lambda genome: genome.genome_id)


def _format_optional_float(value: float | None, digits: int = 4) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def _pair_key(left: str, right: str) -> tuple[str, str]:
    return (left, right) if left <= right else (right, left)


def _all_pairs(items: Iterable[str]) -> list[tuple[str, str]]:
    sorted_items = sorted(items)
    return [(left, right) for left, right in combinations(sorted_items, 2)]


def _parse_mash_dist_output(stdout: str) -> dict[tuple[str, str], float]:
    distances: dict[tuple[str, str], float] = {}

    for raw_line in stdout.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        fields = line.split("\t")
        if len(fields) < 3:
            continue

        query = Path(fields[0]).stem
        reference = Path(fields[1]).stem
        if query == reference:
            continue

        try:
            distance = float(fields[2])
        except ValueError:
            continue

        distances[_pair_key(query, reference)] = distance

    return distances


def _parse_skani_ani_output(stdout: str) -> float:
    """Parse ANI from skani output with tolerant handling across versions."""

    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines:
        raise PanScapeUsageError("skani returned no output; unable to parse ANI")

    header_tokens = re.split(r"\s+", lines[0].lower())
    if "ani" in header_tokens and len(lines) >= 2:
        ani_idx = header_tokens.index("ani")
        value_tokens = re.split(r"\s+", lines[-1])
        if ani_idx < len(value_tokens):
            try:
                ani = float(value_tokens[ani_idx])
                return ani / 100.0 if ani > 1.0 else ani
            except ValueError:
                pass

    float_candidates: list[float] = []
    for line in lines:
        for token in re.split(r"\s+", line):
            try:
                float_candidates.append(float(token))
            except ValueError:
                continue

    if not float_candidates:
        raise PanScapeUsageError("Unable to parse ANI value from skani output")

    candidate = max(value for value in float_candidates if value <= 100.0)
    return candidate / 100.0 if candidate > 1.0 else candidate


def _copy_file(src: Path, dst: Path, *, force: bool) -> Path:
    ensure_dir(dst.parent)
    if dst.exists() and not force:
        raise PanScapeUsageError(f"Refusing to overwrite existing file without --force: {dst}")
    shutil.copyfile(src, dst)
    return dst


def _select_backbone(
    genome_ids: Iterable[str],
    *,
    stats_by_genome: dict[str, AssemblyStats],
    qc_by_genome: dict[str, GenomeQC],
    logger_name: str,
) -> tuple[str, dict[str, float]]:
    logger = get_logger(logger_name)
    scores: dict[str, float] = {}

    for genome_id in sorted(genome_ids):
        qc = qc_by_genome[genome_id]
        stats = stats_by_genome[genome_id]

        if qc.completeness is None or qc.contamination is None:
            logger.warning(
                "Genome %s is missing QC metrics; using neutral values for backbone score.",
                genome_id,
            )

        scores[genome_id] = backbone_score(
            completeness=qc.completeness,
            contamination=qc.contamination,
            n50=stats.n50,
            genome_size=stats.genome_size,
            contig_count=stats.contig_count,
        )

    backbone_genome_id = sorted(scores, key=lambda gid: (-scores[gid], gid))[0]
    return backbone_genome_id, scores


def _mock_checkm2_predictions(contexts: dict[str, GenomeContext]) -> dict[str, tuple[float, float]]:
    # TODO: Replace heuristic QC simulation with a benchmarked surrogate model if
    # mock-mode QC should emulate CheckM2 score distributions more realistically.
    predictions: dict[str, tuple[float, float]] = {}
    for genome_id, context in contexts.items():
        stats = context.assembly_stats
        completeness = min(99.9, 60.0 + (math.log10(max(stats.genome_size, 1)) * 12.0))
        contamination = min(15.0, max(0.1, stats.contig_count * 0.2))
        predictions[genome_id] = (completeness, contamination)
    return predictions


def _parse_checkm2_report(report_path: Path) -> dict[str, tuple[float, float]]:
    with report_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise PanScapeUsageError(f"CheckM2 report has no header row: {report_path}")

        header_by_lower = {name.lower(): name for name in reader.fieldnames}

        name_col = header_by_lower.get("name") or header_by_lower.get("bin")
        comp_col = header_by_lower.get("completeness")
        contam_col = header_by_lower.get("contamination")

        if name_col is None or comp_col is None or contam_col is None:
            raise PanScapeUsageError(
                f"CheckM2 report missing expected columns (name/completeness/contamination): {report_path}"
            )

        predictions: dict[str, tuple[float, float]] = {}
        for row in reader:
            name = (row.get(name_col) or "").strip()
            if name == "":
                continue

            genome_id = Path(name).stem
            try:
                completeness = float((row.get(comp_col) or "").strip())
                contamination = float((row.get(contam_col) or "").strip())
            except ValueError:
                continue

            predictions[genome_id] = (completeness, contamination)

    return predictions


def _find_checkm2_report(checkm2_dir: Path) -> Path:
    candidates = [
        checkm2_dir / "quality_report.tsv",
        checkm2_dir / "predictions.tsv",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate

    tsv_files = sorted(checkm2_dir.rglob("*.tsv"))
    if not tsv_files:
        raise PanScapeUsageError(f"Unable to find CheckM2 report in {checkm2_dir}")

    return tsv_files[0]


def _cluster_proteins_mock(genes: list[CalledGene]) -> list[FamilyCluster]:
    # TODO: Replace mock clustering with a lightweight in-Python protein similarity graph
    # when mmseqs2 is unavailable but higher fidelity than exact-sequence grouping is needed.
    groups: dict[str, list[str]] = {}
    for gene in sorted(genes, key=lambda item: item.gene_id):
        key = gene.aa_sequence or "X"
        groups.setdefault(key, []).append(gene.gene_id)

    clusters: list[FamilyCluster] = []
    for idx, member_gene_ids in enumerate(
        sorted((sorted(members) for members in groups.values()), key=lambda members: members[0]),
        start=1,
    ):
        family_id = f"family_{idx:06d}"
        clusters.append(
            FamilyCluster(
                family_id=family_id,
                representative_gene_id=member_gene_ids[0],
                member_gene_ids=tuple(member_gene_ids),
            )
        )

    return clusters


def _cluster_proteins_mmseqs(
    *,
    species_dir: Path,
    genes: list[CalledGene],
    runner: MMseqsRunner,
    cfg: BuildConfig,
) -> list[FamilyCluster]:
    if len(genes) <= 1:
        return _cluster_proteins_mock(genes)

    all_faa = species_dir / "all_proteins.faa"
    protein_records = [FastaRecord(header=gene.gene_id, sequence=(gene.aa_sequence or "X")) for gene in genes]
    write_fasta_records(all_faa, protein_records, force=cfg.force)

    output_prefix = species_dir / "mmseqs_cluster"
    tmp_dir = ensure_dir(species_dir / "mmseqs_tmp")

    runner.easy_cluster(
        input_faa=all_faa,
        output_prefix=output_prefix,
        tmp_dir=tmp_dir,
        min_seq_id=cfg.mmseqs_min_aa_id,
        coverage=cfg.mmseqs_cov,
        threads=cfg.threads,
        dry_run=cfg.dry_run,
    )

    candidate_cluster_paths = [
        Path(f"{output_prefix}_cluster.tsv"),
        output_prefix.with_suffix(".tsv"),
        species_dir / "mmseqs_cluster.tsv",
    ]
    cluster_path = next((path for path in candidate_cluster_paths if path.exists()), None)
    if cluster_path is None:
        raise PanScapeUsageError(
            f"Unable to find mmseqs cluster TSV output for species directory: {species_dir}"
        )

    clusters_by_rep: dict[str, set[str]] = {}
    with cluster_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) == 1:
                rep = fields[0]
                member = fields[0]
            else:
                rep = fields[0]
                member = fields[1]
            clusters_by_rep.setdefault(rep, set()).add(member)

    all_gene_ids = {gene.gene_id for gene in genes}
    observed_members = {member for members in clusters_by_rep.values() for member in members}
    for missing_gene in sorted(all_gene_ids.difference(observed_members)):
        clusters_by_rep[missing_gene] = {missing_gene}

    clusters: list[FamilyCluster] = []
    for idx, rep_gene_id in enumerate(
        sorted(clusters_by_rep, key=lambda rep: (sorted(clusters_by_rep[rep])[0], rep)),
        start=1,
    ):
        members = tuple(sorted(clusters_by_rep[rep_gene_id]))
        representative = rep_gene_id if rep_gene_id in members else members[0]
        clusters.append(
            FamilyCluster(
                family_id=f"family_{idx:06d}",
                representative_gene_id=representative,
                member_gene_ids=members,
            )
        )

    return clusters


def _family_tier(
    *,
    genomes_present: int,
    denominator: int,
    prevalence_fraction: float,
    rare_min_prevalence: int,
) -> str:
    if genomes_present < rare_min_prevalence:
        return "rare"
    if denominator >= 20 and prevalence_fraction < 0.05:
        return "rare"
    if prevalence_fraction >= 0.95:
        return "core"
    return "accessory"


def _print_plan(step_plan: list[str]) -> None:
    console.print("[bold]Build step plan[/bold]")
    for idx, step in enumerate(step_plan, start=1):
        console.print(f"  {idx}. {step}")


def run_build(
    *,
    config_path: Path | None,
    genomes_tsv: Path | None,
    outdir: Path | None,
    run_checkm2: bool | None,
    checkm2_db: Path | None,
    min_completeness: float | None,
    max_contamination: float | None,
    mash_sketch_size: int | None,
    mash_k: int | None,
    mash_threshold: float | None,
    species_ani: float | None,
    strain_ani: float | None,
    mmseqs_min_aa_id: float | None,
    mmseqs_cov: float | None,
    within_family_derep_nt: float | None,
    mappability_k: int | None,
    min_gene_len: int | None,
    rare_min_prevalence: int | None,
    use_high_quality_for_prevalence: bool | None,
    hq_min_completeness: float | None,
    hq_max_contamination: float | None,
    keep_singletons: bool | None,
    mock: bool | None,
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
                "run_checkm2": run_checkm2,
                "checkm2_db": checkm2_db,
                "min_completeness": min_completeness,
                "max_contamination": max_contamination,
                "mash_sketch_size": mash_sketch_size,
                "mash_k": mash_k,
                "mash_threshold": mash_threshold,
                "species_ani": species_ani,
                "strain_ani": strain_ani,
                "mmseqs_min_aa_id": mmseqs_min_aa_id,
                "mmseqs_cov": mmseqs_cov,
                "within_family_derep_nt": within_family_derep_nt,
                "mappability_k": mappability_k,
                "min_gene_len": min_gene_len,
                "rare_min_prevalence": rare_min_prevalence,
                "use_high_quality_for_prevalence": use_high_quality_for_prevalence,
                "hq_min_completeness": hq_min_completeness,
                "hq_max_contamination": hq_max_contamination,
                "keep_singletons": keep_singletons,
                "mock": mock,
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
        if cfg.checkm2_db is not None and not cfg.checkm2_db.exists():
            raise PanScapeUsageError(f"CheckM2 DB path does not exist: {cfg.checkm2_db}")

        if outdir is None and cfg.outdir == Path("panscape_out") and config_path is None:
            raise PanScapeUsageError("Missing output directory. Provide --outdir or set build.outdir in config.")

        genomes = parse_genomes_manifest(cfg.genomes_tsv)
        if len(genomes) == 0:
            raise PanScapeUsageError("No genomes available for build.")

        layout = create_output_layout(cfg.outdir)
        build_dir = ensure_dir(layout.build_dir)

        if build_dir.exists() and any(build_dir.iterdir()) and not cfg.force and not cfg.dry_run:
            raise PanScapeUsageError(
                f"Build output directory already contains files: {build_dir}. Use --force to overwrite."
            )

        qc_dir = ensure_dir(build_dir / "qc")
        mash_dir = ensure_dir(build_dir / "mash")
        mash_sketch_dir = ensure_dir(mash_dir / "mash_sketches")
        ani_dir = ensure_dir(build_dir / "ani")
        clusters_dir = ensure_dir(build_dir / "clusters")
        backbones_dir = ensure_dir(build_dir / "backbones")
        genes_dir = ensure_dir(build_dir / "genes")
        normalized_dir = ensure_dir(build_dir / "normalized_fastas")
        pangenomes_dir = ensure_dir(build_dir / "pangenomes")

        mash_runner = MashRunner()
        skani_runner = SkaniRunner()
        mmseqs_runner = MMseqsRunner()
        checkm2_runner = CheckM2Runner()

        missing_qc_ids = [
            genome.genome_id
            for genome in genomes
            if genome.completeness is None or genome.contamination is None
        ]

        tool_versions: dict[str, str] = {}
        if cfg.mock:
            tool_versions = {
                "mash": "mock",
                "skani": "mock",
                "mmseqs": "mock",
                "checkm2": "mock" if cfg.run_checkm2 and missing_qc_ids else "not-used",
            }
        else:
            required_tools = [
                ("mash", mash_runner),
                ("skani", skani_runner),
                ("mmseqs", mmseqs_runner),
            ]
            for tool_name, runner in required_tools:
                if not runner.is_available():
                    raise PanScapeUsageError(
                        f"Required external tool not found in PATH: {tool_name}. "
                        f"Install {tool_name} or run with --mock."
                    )
                tool_versions[tool_name] = runner.version(dry_run=cfg.dry_run)

            if cfg.run_checkm2 and missing_qc_ids:
                if not checkm2_runner.is_available():
                    raise PanScapeUsageError(
                        "--run-checkm2 requested but `checkm2` binary is not available in PATH."
                    )
                tool_versions["checkm2"] = checkm2_runner.version(dry_run=cfg.dry_run)
            else:
                tool_versions["checkm2"] = "not-used"

        step_plan = [
            "Parse and validate genomes.tsv",
            "Normalize FASTA headers and compute assembly stats",
            "Compute genome QC (input values and optional CheckM2)",
            "Run Mash preclustering",
            "Run skani ANI edges for species/strain thresholds",
            "Build deterministic species and strain clusters",
            "Select species backbone genomes and write backbone FASTAs",
            "Call genes per genome using pyrodigal (or mock caller)",
            "Cluster proteins into gene families via mmseqs2 (or mock clustering)",
            "Compute prevalence tiers and within-family dereplication",
            "Compute gene mappability and final tiers",
            "Write per-species catalog summaries and mapper index plans",
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
            parameters=cfg.model_dump(mode="json"),
            tool_versions=tool_versions,
        )
        write_manifest(layout.root, manifest)

        _print_plan(step_plan)

        if cfg.dry_run:
            logger.info("Dry-run requested; stopping before build outputs are written.")
            finalize_manifest(manifest, status="dry-run", output_paths=[])
            write_manifest(layout.root, manifest)
            return 0

        output_paths: list[Path] = []

        # Step 2: FASTA normalization + stats
        contexts: dict[str, GenomeContext] = {}
        assembly_rows: list[list[str]] = []

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TimeElapsedColumn(),
            console=console,
        ) as progress:
            task = progress.add_task("Normalizing FASTAs", total=len(genomes))
            for genome in genomes:
                normalized_fasta = normalized_dir / f"{genome.genome_id}.fna"
                normalized_records, stats = normalize_fasta_file(
                    genome.genome_id,
                    genome.fasta_path,
                    normalized_fasta,
                    force=cfg.force,
                )
                contexts[genome.genome_id] = GenomeContext(
                    genome=genome,
                    normalized_fasta=normalized_fasta,
                    normalized_records=normalized_records,
                    assembly_stats=stats,
                )
                assembly_rows.append(
                    [
                        genome.genome_id,
                        str(stats.genome_size),
                        str(stats.n50),
                        str(stats.contig_count),
                    ]
                )
                progress.advance(task)

        assembly_stats_path = write_tsv(
            qc_dir / "assembly_stats.tsv",
            ["genome_id", "genome_size", "n50", "contig_count"],
            sorted(assembly_rows, key=lambda row: row[0]),
            force=cfg.force,
        )
        output_paths.append(assembly_stats_path)

        # Step 4: QC table + filtering
        checkm2_predictions: dict[str, tuple[float, float]] = {}
        if missing_qc_ids and cfg.run_checkm2:
            if cfg.mock:
                checkm2_predictions = _mock_checkm2_predictions(contexts)
            else:
                checkm2_dir = ensure_dir(qc_dir / "checkm2")
                checkm2_runner.predict(
                    input_dir=normalized_dir,
                    output_dir=checkm2_dir,
                    threads=cfg.threads,
                    database_path=cfg.checkm2_db,
                    dry_run=cfg.dry_run,
                )
                report_path = _find_checkm2_report(checkm2_dir)
                checkm2_predictions = _parse_checkm2_report(report_path)

        qc_by_genome: dict[str, GenomeQC] = {}
        qc_rows: list[list[str]] = []
        kept_genomes: list[str] = []

        for genome in genomes:
            completeness = genome.completeness
            contamination = genome.contamination
            qc_source = "input"

            if completeness is None or contamination is None:
                prediction = checkm2_predictions.get(genome.genome_id)
                if prediction is not None:
                    completeness, contamination = prediction
                    qc_source = "checkm2"
                else:
                    qc_source = "missing"
                    if cfg.run_checkm2:
                        logger.warning(
                            "CheckM2 did not provide QC values for genome %s; keeping with missing QC.",
                            genome.genome_id,
                        )

            passes_thresholds: bool | None
            keep_for_clustering: bool
            if completeness is None or contamination is None:
                passes_thresholds = None
                keep_for_clustering = True
                logger.warning(
                    "Genome %s has missing QC values and is kept for downstream clustering.",
                    genome.genome_id,
                )
            else:
                passes_thresholds = (
                    completeness >= cfg.min_completeness
                    and contamination <= cfg.max_contamination
                )
                keep_for_clustering = passes_thresholds

            if keep_for_clustering:
                kept_genomes.append(genome.genome_id)

            qc_entry = GenomeQC(
                completeness=completeness,
                contamination=contamination,
                qc_source=qc_source,
                passes_thresholds=passes_thresholds,
                keep_for_clustering=keep_for_clustering,
            )
            qc_by_genome[genome.genome_id] = qc_entry

            qc_rows.append(
                [
                    genome.genome_id,
                    _format_optional_float(completeness, digits=3),
                    _format_optional_float(contamination, digits=3),
                    qc_source,
                    "" if passes_thresholds is None else ("1" if passes_thresholds else "0"),
                    "1" if keep_for_clustering else "0",
                ]
            )

        genome_qc_path = write_tsv(
            qc_dir / "genome_qc.tsv",
            [
                "genome_id",
                "completeness",
                "contamination",
                "qc_source",
                "passes_thresholds",
                "keep_for_clustering",
            ],
            sorted(qc_rows, key=lambda row: row[0]),
            force=cfg.force,
        )
        output_paths.append(genome_qc_path)

        kept_genomes = sorted(set(kept_genomes))
        if len(kept_genomes) == 0:
            raise PanScapeUsageError(
                "No genomes passed QC filters. Relax thresholds or provide higher-quality assemblies."
            )

        # Step 5: Mash preclustering
        mash_distances: dict[tuple[str, str], float] = {}

        if len(kept_genomes) >= 2:
            if cfg.mock:
                concatenated_sequences = {
                    genome_id: "".join(record.sequence for record in contexts[genome_id].normalized_records)
                    for genome_id in kept_genomes
                }
                for left, right in _all_pairs(kept_genomes):
                    mash_distances[(left, right)] = mock_mash_distance(
                        concatenated_sequences[left],
                        concatenated_sequences[right],
                        k=cfg.mash_k,
                        sketch_size=cfg.mash_sketch_size,
                    )
                output_paths.append(
                    write_text(
                        mash_sketch_dir / "MOCK_SKETCHES.txt",
                        "Mock Mash mode: sketches approximated with deterministic k-mer minhash.\n",
                        force=cfg.force,
                    )
                )
            else:
                sketch_prefix = mash_sketch_dir / "genomes"
                mash_runner.sketch(
                    input_fastas=[contexts[genome_id].normalized_fasta for genome_id in kept_genomes],
                    out_prefix=sketch_prefix,
                    kmer_size=cfg.mash_k,
                    sketch_size=cfg.mash_sketch_size,
                    dry_run=cfg.dry_run,
                )
                dist_result = mash_runner.dist(
                    sketch_path=Path(f"{sketch_prefix}.msh"),
                    dry_run=cfg.dry_run,
                )
                mash_distances = _parse_mash_dist_output(dist_result.stdout)

        mash_edges = [pair for pair, dist in mash_distances.items() if dist <= cfg.mash_threshold]
        precluster_components = connected_components(kept_genomes, mash_edges)

        precluster_id_by_members: dict[tuple[str, ...], str] = {}
        for idx, component in enumerate(precluster_components, start=1):
            precluster_id_by_members[tuple(component)] = f"precluster_{idx:06d}"

        genome_to_precluster: dict[str, str] = {}
        precluster_rows: list[list[str]] = []
        for component in precluster_components:
            precluster_id = precluster_id_by_members[tuple(component)]
            for genome_id in component:
                genome_to_precluster[genome_id] = precluster_id
                precluster_rows.append([precluster_id, genome_id])

        preclusters_path = write_tsv(
            mash_dir / "preclusters.tsv",
            ["precluster_id", "genome_id"],
            sorted(precluster_rows, key=lambda row: (row[0], row[1])),
            force=cfg.force,
        )
        output_paths.append(preclusters_path)

        # Step 6: skani ANI edges
        ani_jobs: list[tuple[str, str, str]] = []
        for component in precluster_components:
            precluster_id = precluster_id_by_members[tuple(component)]
            for left, right in _all_pairs(component):
                ani_jobs.append((precluster_id, left, right))

        ani_scores: dict[tuple[str, str], float] = {}

        def _compute_ani(job: tuple[str, str, str]) -> tuple[tuple[str, str], float]:
            _, left, right = job
            pair = _pair_key(left, right)
            if cfg.mock:
                # TODO: Swap this approximation for a calibrated regression against real skani ANI
                # if mock-mode behavior needs to mimic production edge densities more closely.
                distance = mash_distances.get(pair)
                if distance is None:
                    distance = mock_mash_distance(
                        "".join(record.sequence for record in contexts[left].normalized_records),
                        "".join(record.sequence for record in contexts[right].normalized_records),
                        k=cfg.mash_k,
                        sketch_size=cfg.mash_sketch_size,
                    )
                ani = max(0.0, min(1.0, 1.0 - distance))
                return pair, ani

            result = skani_runner.dist_pair(
                query=contexts[left].normalized_fasta,
                reference=contexts[right].normalized_fasta,
                dry_run=cfg.dry_run,
            )
            ani_value = _parse_skani_ani_output(result.stdout)
            return pair, ani_value

        if ani_jobs:
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TimeElapsedColumn(),
                console=console,
            ) as progress:
                task = progress.add_task("Computing ANI edges", total=len(ani_jobs))
                with ThreadPoolExecutor(max_workers=cfg.threads) as pool:
                    futures = [pool.submit(_compute_ani, job) for job in ani_jobs]
                    for future in as_completed(futures):
                        pair, ani_value = future.result()
                        ani_scores[pair] = ani_value
                        progress.advance(task)

        species_edge_rows: list[list[str]] = []
        strain_edge_rows: list[list[str]] = []

        for precluster_id, left, right in sorted(ani_jobs, key=lambda item: (item[0], item[1], item[2])):
            ani_value = ani_scores[_pair_key(left, right)]
            row = [precluster_id, left, right, f"{ani_value:.6f}"]
            if ani_value >= cfg.species_ani:
                species_edge_rows.append(row)
            if ani_value >= cfg.strain_ani:
                strain_edge_rows.append(row)

        species_edges_path = write_tsv(
            ani_dir / "skani_edges_species.tsv",
            ["precluster_id", "genome_a", "genome_b", "ani"],
            species_edge_rows,
            force=cfg.force,
        )
        strain_edges_path = write_tsv(
            ani_dir / "skani_edges_strain.tsv",
            ["precluster_id", "genome_a", "genome_b", "ani"],
            strain_edge_rows,
            force=cfg.force,
        )
        output_paths.extend([species_edges_path, strain_edges_path])

        # Step 7 + 8: species/strain components and backbones
        species_components = connected_components(
            kept_genomes,
            [(row[1], row[2]) for row in species_edge_rows],
        )

        stats_by_genome = {genome_id: contexts[genome_id].assembly_stats for genome_id in kept_genomes}

        component_backbone_info: list[tuple[tuple[str, ...], str, dict[str, float]]] = []
        for component in species_components:
            backbone_genome_id, scores = _select_backbone(
                component,
                stats_by_genome=stats_by_genome,
                qc_by_genome=qc_by_genome,
                logger_name="panscape.build.backbone",
            )
            component_backbone_info.append((tuple(component), backbone_genome_id, scores))

        component_backbone_info.sort(key=lambda item: (item[1], item[0]))

        species_by_component: dict[tuple[str, ...], str] = {}
        component_scores: dict[tuple[str, ...], dict[str, float]] = {}
        backbone_by_species: dict[str, str] = {}

        for idx, (component, backbone_genome_id, scores) in enumerate(component_backbone_info, start=1):
            species_id = f"species_{idx:06d}"
            species_by_component[component] = species_id
            component_scores[component] = scores
            backbone_by_species[species_id] = backbone_genome_id

        species_members: dict[str, tuple[str, ...]] = {
            species_by_component[component]: component
            for component in species_by_component
        }

        species_cluster_rows: list[list[str]] = []
        for species_id, genomes_in_species in sorted(species_members.items()):
            backbone_id = backbone_by_species[species_id]
            for genome_id in genomes_in_species:
                species_cluster_rows.append(
                    [
                        species_id,
                        genome_id,
                        "1" if genome_id == backbone_id else "0",
                    ]
                )

        species_clusters_path = write_tsv(
            clusters_dir / "species_clusters.tsv",
            ["species_id", "genome_id", "is_backbone"],
            species_cluster_rows,
            force=cfg.force,
        )
        output_paths.append(species_clusters_path)

        strain_bin_rows: list[list[str]] = []
        strain_edge_pairs = {(row[1], row[2]) for row in strain_edge_rows}
        reverse_pairs = {(right, left) for left, right in strain_edge_pairs}
        strain_edge_pairs |= reverse_pairs

        for species_id, genomes_in_species in sorted(species_members.items()):
            species_nodes = list(genomes_in_species)
            species_strain_edges = [
                (left, right)
                for left, right in combinations(sorted(species_nodes), 2)
                if (left, right) in strain_edge_pairs
            ]
            strain_components = connected_components(species_nodes, species_strain_edges)
            for strain_idx, strain_component in enumerate(strain_components, start=1):
                strain_bin_id = f"strain_{strain_idx:06d}"
                for genome_id in strain_component:
                    strain_bin_rows.append([species_id, strain_bin_id, genome_id])

        strain_bins_path = write_tsv(
            clusters_dir / "strain_bins.tsv",
            ["species_id", "strain_bin_id", "genome_id"],
            sorted(strain_bin_rows, key=lambda row: (row[0], row[1], row[2])),
            force=cfg.force,
        )
        output_paths.append(strain_bins_path)

        backbone_rows: list[list[str]] = []
        for species_id, genomes_in_species in sorted(species_members.items()):
            backbone_id = backbone_by_species[species_id]
            component = tuple(genomes_in_species)
            score = component_scores[component][backbone_id]
            qc = qc_by_genome[backbone_id]
            stats = stats_by_genome[backbone_id]

            backbone_rows.append(
                [
                    species_id,
                    backbone_id,
                    f"{score:.6f}",
                    _format_optional_float(qc.completeness, digits=3),
                    _format_optional_float(qc.contamination, digits=3),
                    str(stats.n50),
                    str(stats.genome_size),
                    str(stats.contig_count),
                ]
            )

            backbone_fasta_path = backbones_dir / f"{species_id}.fna"
            output_paths.append(
                _copy_file(
                    contexts[backbone_id].normalized_fasta,
                    backbone_fasta_path,
                    force=cfg.force,
                )
            )

        backbone_table_path = write_tsv(
            backbones_dir / "backbone_table.tsv",
            [
                "species_id",
                "backbone_genome_id",
                "backbone_score",
                "completeness",
                "contamination",
                "n50",
                "genome_size",
                "contig_count",
            ],
            sorted(backbone_rows, key=lambda row: row[0]),
            force=cfg.force,
        )
        output_paths.append(backbone_table_path)

        # Step 9: gene calling
        genes_by_genome: dict[str, list[CalledGene]] = {}

        def _call_genes_for_genome(genome_id: str) -> tuple[str, list[CalledGene], tuple[Path, Path, Path]]:
            context = contexts[genome_id]
            genes = call_genes(
                genome_id=genome_id,
                records=context.normalized_records,
                min_gene_len=cfg.min_gene_len,
                mock=cfg.mock,
            )
            paths = write_gene_outputs(
                genome_id=genome_id,
                genes=genes,
                genes_dir=genes_dir,
                force=cfg.force,
            )
            return genome_id, genes, paths

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TimeElapsedColumn(),
            console=console,
        ) as progress:
            task = progress.add_task("Calling genes", total=len(kept_genomes))
            with ThreadPoolExecutor(max_workers=cfg.threads) as pool:
                futures = [pool.submit(_call_genes_for_genome, genome_id) for genome_id in kept_genomes]
                for future in as_completed(futures):
                    genome_id, genes, paths = future.result()
                    genes_by_genome[genome_id] = genes
                    output_paths.extend(list(paths))
                    progress.advance(task)

        # Step 10-14: per-species pangenome catalogs
        for species_id, genomes_in_species in sorted(species_members.items()):
            species_dir = ensure_dir(pangenomes_dir / species_id)
            species_genes = [
                gene
                for genome_id in genomes_in_species
                for gene in genes_by_genome.get(genome_id, [])
            ]
            genes_by_id = {gene.gene_id: gene for gene in species_genes}

            if cfg.mock:
                family_clusters = _cluster_proteins_mock(species_genes)
            else:
                family_clusters = _cluster_proteins_mmseqs(
                    species_dir=species_dir,
                    genes=species_genes,
                    runner=mmseqs_runner,
                    cfg=cfg,
                )

            family_by_gene: dict[str, str] = {}
            for cluster in family_clusters:
                for gene_id in cluster.member_gene_ids:
                    family_by_gene[gene_id] = cluster.family_id

            gene_families_rows = [
                [gene_id, family_id, genes_by_id[gene_id].genome_id]
                for gene_id, family_id in sorted(
                    family_by_gene.items(),
                    key=lambda item: (item[1], item[0]),
                )
            ]
            gene_families_path = write_tsv(
                species_dir / "gene_families.tsv",
                ["gene_id", "family_id", "genome_id"],
                gene_families_rows,
                force=cfg.force,
            )
            output_paths.append(gene_families_path)

            representative_faa_records: list[FastaRecord] = []
            for cluster in sorted(family_clusters, key=lambda item: item.family_id):
                representative = genes_by_id.get(cluster.representative_gene_id)
                if representative is None:
                    continue
                representative_faa_records.append(
                    FastaRecord(
                        header=cluster.representative_gene_id,
                        sequence=representative.aa_sequence or "X",
                    )
                )

            representatives_faa_path = species_dir / "representatives.faa"
            write_fasta_records(representatives_faa_path, representative_faa_records, force=cfg.force)
            output_paths.append(representatives_faa_path)

            # Family prevalence
            species_genome_set = set(genomes_in_species)
            hq_genomes = {
                genome_id
                for genome_id in genomes_in_species
                if (
                    qc_by_genome[genome_id].completeness is not None
                    and qc_by_genome[genome_id].contamination is not None
                    and qc_by_genome[genome_id].completeness >= cfg.hq_min_completeness
                    and qc_by_genome[genome_id].contamination <= cfg.hq_max_contamination
                )
            }

            if cfg.use_high_quality_for_prevalence and hq_genomes:
                prevalence_denominator = hq_genomes
                prevalence_mode = "high_quality"
            else:
                prevalence_denominator = species_genome_set
                prevalence_mode = "all_genomes"

            denominator_size = len(prevalence_denominator)

            family_to_member_genomes: dict[str, set[str]] = {cluster.family_id: set() for cluster in family_clusters}
            for gene_id, family_id in family_by_gene.items():
                family_to_member_genomes[family_id].add(genes_by_id[gene_id].genome_id)

            family_prevalence_rows: list[list[str]] = []
            family_prevalence_map: dict[str, tuple[int, int, float]] = {}

            for family_id in sorted(family_to_member_genomes):
                present_all = family_to_member_genomes[family_id]
                present_in_denominator = len(present_all & prevalence_denominator)
                fraction = (
                    float(present_in_denominator / denominator_size)
                    if denominator_size > 0
                    else 0.0
                )
                family_prevalence_map[family_id] = (present_in_denominator, denominator_size, fraction)
                family_prevalence_rows.append(
                    [
                        family_id,
                        str(present_in_denominator),
                        str(denominator_size),
                        f"{fraction:.6f}",
                        prevalence_mode,
                    ]
                )

            family_prevalence_path = write_tsv(
                species_dir / "family_prevalence.tsv",
                [
                    "family_id",
                    "genomes_present",
                    "denominator_genomes",
                    "prevalence_fraction",
                    "prevalence_mode",
                ],
                family_prevalence_rows,
                force=cfg.force,
            )
            output_paths.append(family_prevalence_path)

            # Within-family dereplication
            backbone_genome_id = backbone_by_species[species_id]

            def _preference_key(member: DerepGene) -> tuple:
                qc = qc_by_genome[member.genome_id]
                qc_score = (
                    (qc.completeness if qc.completeness is not None else 85.0)
                    - (5.0 * (qc.contamination if qc.contamination is not None else 5.0))
                )
                return (
                    0 if member.genome_id == backbone_genome_id else 1,
                    -qc_score,
                    -len(member.nt_sequence),
                    member.gene_id,
                )

            derep_groups_by_family: dict[str, list[DerepGroup]] = {}
            for cluster in sorted(family_clusters, key=lambda item: item.family_id):
                members: list[DerepGene] = []
                for gene_id in cluster.member_gene_ids:
                    gene = genes_by_id[gene_id]
                    members.append(
                        DerepGene(
                            gene_id=gene.gene_id,
                            genome_id=gene.genome_id,
                            family_id=cluster.family_id,
                            nt_sequence=gene.nt_sequence,
                            aa_sequence=gene.aa_sequence,
                        )
                    )

                derep_groups_by_family[cluster.family_id] = dereplicate_family(
                    members,
                    identity_threshold=cfg.within_family_derep_nt,
                    k=cfg.mappability_k,
                    preference_key=_preference_key,
                )

            representative_rows: list[list[str]] = []
            selected_representative_gene_ids: list[str] = []
            representative_family_by_gene: dict[str, str] = {}

            for family_id in sorted(derep_groups_by_family):
                groups = derep_groups_by_family[family_id]
                family_member_genome_count = len(family_to_member_genomes.get(family_id, set()))
                include_family = cfg.keep_singletons or family_member_genome_count > 1

                for group in groups:
                    representative_gene_id = group.representative_gene_id
                    representative_gene = genes_by_id[representative_gene_id]
                    is_backbone_rep = representative_gene.genome_id == backbone_genome_id
                    representative_rows.append(
                        [
                            family_id,
                            group.group_id,
                            representative_gene_id,
                            str(len(group.member_gene_ids)),
                            ",".join(group.member_gene_ids),
                            "1" if is_backbone_rep else "0",
                            "1" if include_family else "0",
                        ]
                    )

                    if include_family:
                        selected_representative_gene_ids.append(representative_gene_id)
                        representative_family_by_gene[representative_gene_id] = family_id

            selected_representative_gene_ids = sorted(set(selected_representative_gene_ids))

            gene_representatives_path = write_tsv(
                species_dir / "gene_representatives.tsv",
                [
                    "family_id",
                    "group_id",
                    "rep_gene_id",
                    "member_count",
                    "member_gene_ids",
                    "is_backbone_rep",
                    "included_in_catalog",
                ],
                sorted(representative_rows, key=lambda row: (row[0], row[1], row[2])),
                force=cfg.force,
            )
            output_paths.append(gene_representatives_path)

            genes_fna_records = [
                FastaRecord(header=gene_id, sequence=genes_by_id[gene_id].nt_sequence)
                for gene_id in selected_representative_gene_ids
            ]
            genes_faa_records = [
                FastaRecord(header=gene_id, sequence=(genes_by_id[gene_id].aa_sequence or "X"))
                for gene_id in selected_representative_gene_ids
            ]

            genes_fna_path = species_dir / "genes.fna"
            genes_faa_path = species_dir / "genes.faa"
            write_fasta_records(genes_fna_path, genes_fna_records, force=cfg.force)
            write_fasta_records(genes_faa_path, genes_faa_records, force=cfg.force)
            output_paths.extend([genes_fna_path, genes_faa_path])

            mappability_results = unique_kmer_fractions(
                {gene_id: genes_by_id[gene_id].nt_sequence for gene_id in selected_representative_gene_ids},
                cfg.mappability_k,
            )

            gene_mappability_rows = [
                [
                    gene_id,
                    str(result.total_kmers),
                    str(result.unique_kmers),
                    f"{result.unique_fraction:.6f}",
                ]
                for gene_id, result in sorted(mappability_results.items())
            ]

            gene_mappability_path = write_tsv(
                species_dir / "gene_mappability.tsv",
                ["rep_gene_id", "total_kmers", "unique_kmers", "unique_fraction"],
                gene_mappability_rows,
                force=cfg.force,
            )
            output_paths.append(gene_mappability_path)

            family_tier_by_family: dict[str, str] = {}
            for family_id, (present_count, denominator, prevalence_fraction) in family_prevalence_map.items():
                family_tier_by_family[family_id] = _family_tier(
                    genomes_present=present_count,
                    denominator=denominator,
                    prevalence_fraction=prevalence_fraction,
                    rare_min_prevalence=cfg.rare_min_prevalence,
                )

            gene_tier_rows: list[list[str]] = []
            for gene_id in selected_representative_gene_ids:
                family_id = representative_family_by_gene[gene_id]
                present_count, _, prevalence_fraction = family_prevalence_map[family_id]
                family_tier = family_tier_by_family[family_id]
                mappability_fraction = mappability_results[gene_id].unique_fraction
                strain_ok = (
                    mappability_fraction >= cfg.strain_ok_min_mappability
                    and family_tier != "rare"
                )
                gene_tier_rows.append(
                    [
                        gene_id,
                        family_id,
                        family_tier,
                        str(present_count),
                        f"{prevalence_fraction:.6f}",
                        f"{mappability_fraction:.6f}",
                        "1" if strain_ok else "0",
                    ]
                )

            gene_tiers_path = write_tsv(
                species_dir / "gene_tiers.tsv",
                [
                    "rep_gene_id",
                    "family_id",
                    "family_tier",
                    "genomes_present",
                    "prevalence_fraction",
                    "unique_kmer_fraction",
                    "strain_ok",
                ],
                gene_tier_rows,
                force=cfg.force,
            )
            output_paths.append(gene_tiers_path)

            summary_rows = [
                ["species_id", species_id],
                ["genomes_in_species", str(len(genomes_in_species))],
                ["total_called_genes", str(len(species_genes))],
                ["total_families", str(len(family_clusters))],
                ["core_families", str(sum(1 for tier in family_tier_by_family.values() if tier == "core"))],
                [
                    "accessory_families",
                    str(sum(1 for tier in family_tier_by_family.values() if tier == "accessory")),
                ],
                ["rare_families", str(sum(1 for tier in family_tier_by_family.values() if tier == "rare"))],
                ["catalog_representatives", str(len(selected_representative_gene_ids))],
                [
                    "strain_ok_representatives",
                    str(sum(1 for row in gene_tier_rows if row[-1] == "1")),
                ],
            ]

            catalog_summary_path = write_tsv(
                species_dir / "catalog_summary.tsv",
                ["metric", "value"],
                summary_rows,
                force=cfg.force,
            )
            output_paths.append(catalog_summary_path)

            index_plan = {
                "species_id": species_id,
                "backbone_fasta": str(backbones_dir / f"{species_id}.fna"),
                "genes_fna": str(genes_fna_path),
                "genes_faa": str(genes_faa_path),
                "recommended_mappers": ["minimap2", "bwa"],
                "notes": [
                    "Build mapper indices in `panscape map` based on this plan.",
                    "Backbone and gene-catalog references should be indexed separately.",
                ],
            }
            index_plan_path = write_json(
                species_dir / "index" / "INDEX_PLAN.json",
                index_plan,
                force=cfg.force,
            )
            output_paths.append(index_plan_path)

        finalize_manifest(manifest, status="completed", output_paths=output_paths)
        write_manifest(layout.root, manifest)

        logger.info("Build pipeline completed successfully.")
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
        help="Genome manifest TSV with columns: genome_id, fasta_path, completeness, contamination.",
    ),
    outdir: Path | None = typer.Option(None, "--outdir", help="Output root directory."),
    run_checkm2: bool | None = typer.Option(
        None,
        "--run-checkm2/--no-run-checkm2",
        help="Run CheckM2 for genomes missing completeness/contamination.",
    ),
    checkm2_db: Path | None = typer.Option(None, "--checkm2-db", help="Optional CheckM2 DB path."),
    min_completeness: float | None = typer.Option(None, "--min-completeness", help="Minimum completeness threshold."),
    max_contamination: float | None = typer.Option(None, "--max-contamination", help="Maximum contamination threshold."),
    mash_sketch_size: int | None = typer.Option(None, "--mash-sketch-size", min=1, help="Mash sketch size."),
    mash_k: int | None = typer.Option(None, "--mash-k", min=1, help="Mash k-mer size."),
    mash_threshold: float | None = typer.Option(None, "--mash-threshold", min=0.0, max=1.0, help="Mash distance threshold for preclusters."),
    species_ani: float | None = typer.Option(None, "--species-ani", min=0.0, max=1.0, help="Species ANI threshold."),
    strain_ani: float | None = typer.Option(None, "--strain-ani", min=0.0, max=1.0, help="Strain ANI threshold."),
    mmseqs_min_aa_id: float | None = typer.Option(None, "--mmseqs-min-aa-id", min=0.0, max=1.0, help="mmseqs minimum amino-acid identity."),
    mmseqs_cov: float | None = typer.Option(None, "--mmseqs-cov", min=0.0, max=1.0, help="mmseqs coverage threshold."),
    within_family_derep_nt: float | None = typer.Option(None, "--within-family-derep-nt", min=0.0, max=1.0, help="Within-family nucleotide derep threshold."),
    mappability_k: int | None = typer.Option(None, "--mappability-k", min=1, help="k-mer length for mappability scoring."),
    min_gene_len: int | None = typer.Option(None, "--min-gene-len", min=30, help="Minimum nucleotide CDS length."),
    rare_min_prevalence: int | None = typer.Option(None, "--rare-min-prevalence", min=1, help="Families with prevalence below this count are rare."),
    use_high_quality_for_prevalence: bool | None = typer.Option(
        None,
        "--use-high-quality-for-prevalence/--no-use-high-quality-for-prevalence",
        help="Use HQ genomes as prevalence denominator when available.",
    ),
    hq_min_completeness: float | None = typer.Option(None, "--hq-min-completeness", min=0.0, max=100.0, help="HQ completeness threshold."),
    hq_max_contamination: float | None = typer.Option(None, "--hq-max-contamination", min=0.0, max=100.0, help="HQ contamination threshold."),
    keep_singletons: bool | None = typer.Option(
        None,
        "--keep-singletons/--drop-singletons",
        help="Keep or drop singleton-family representatives.",
    ),
    mock: bool | None = typer.Option(None, "--mock/--no-mock", help="Mock mode without external binaries."),
    threads: int | None = typer.Option(None, "--threads", min=1, help="Worker threads."),
    dry_run: bool | None = typer.Option(None, "--dry-run", help="Plan only, do not write outputs."),
    force: bool | None = typer.Option(None, "--force", help="Overwrite existing output files."),
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
        run_checkm2=run_checkm2,
        checkm2_db=checkm2_db,
        min_completeness=min_completeness,
        max_contamination=max_contamination,
        mash_sketch_size=mash_sketch_size,
        mash_k=mash_k,
        mash_threshold=mash_threshold,
        species_ani=species_ani,
        strain_ani=strain_ani,
        mmseqs_min_aa_id=mmseqs_min_aa_id,
        mmseqs_cov=mmseqs_cov,
        within_family_derep_nt=within_family_derep_nt,
        mappability_k=mappability_k,
        min_gene_len=min_gene_len,
        rare_min_prevalence=rare_min_prevalence,
        use_high_quality_for_prevalence=use_high_quality_for_prevalence,
        hq_min_completeness=hq_min_completeness,
        hq_max_contamination=hq_max_contamination,
        keep_singletons=keep_singletons,
        mock=mock,
        threads=threads,
        dry_run=dry_run,
        force=force,
        log_file=log_file,
        verbose=verbose,
        quiet=quiet,
    )
    raise typer.Exit(exit_code)
