from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from panscape.exceptions import PanScapeUsageError

FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")
FASTQ_SUFFIXES = (".fq", ".fastq", ".fq.gz", ".fastq.gz")


@dataclass(frozen=True, slots=True)
class GenomeRecord:
    genome_id: str
    fasta_path: Path
    completeness: float | None = None
    contamination: float | None = None
    sample_id: str | None = None


@dataclass(frozen=True, slots=True)
class SampleRecord:
    sample_id: str
    r1: Path
    r2: Path | None = None


def _matches_suffix(path: Path, suffixes: Iterable[str]) -> bool:
    lowered = path.name.lower()
    return any(lowered.endswith(suffix) for suffix in suffixes)


def _resolve_manifest_path(base_dir: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


def validate_existing_file(path: Path, label: str) -> None:
    if not path.exists():
        raise PanScapeUsageError(f"{label} does not exist: {path}")
    if not path.is_file():
        raise PanScapeUsageError(f"{label} is not a file: {path}")


def validate_optional_file(path: Path | None, label: str) -> None:
    if path is None:
        return
    validate_existing_file(path, label)


def _read_manifest(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise PanScapeUsageError(f"Manifest TSV has no header row: {path}")

        header = [name.strip() for name in reader.fieldnames]
        rows: list[dict[str, str]] = []
        for row in reader:
            normalized = {str(key): (value or "").strip() for key, value in row.items()}
            rows.append(normalized)

    return header, rows


def _require_columns(path: Path, header: list[str], required: set[str]) -> None:
    missing = required.difference(header)
    if missing:
        missing_cols = ", ".join(sorted(missing))
        raise PanScapeUsageError(f"Missing required columns in {path}: {missing_cols}")


def _parse_optional_float(raw: str, field_name: str, row_idx: int) -> float | None:
    if raw == "":
        return None
    try:
        return float(raw)
    except ValueError as exc:
        raise PanScapeUsageError(
            f"Invalid float in column `{field_name}` at row {row_idx + 2}: {raw!r}"
        ) from exc


def validate_genomes_manifest(manifest_path: Path) -> list[GenomeRecord]:
    validate_existing_file(manifest_path, "genomes.tsv")
    header, rows = _read_manifest(manifest_path)
    _require_columns(manifest_path, header, {"genome_id", "fasta_path"})

    if not rows:
        raise PanScapeUsageError(f"genomes.tsv has no data rows: {manifest_path}")

    records: list[GenomeRecord] = []
    for idx, row in enumerate(rows):
        genome_id = row.get("genome_id", "")
        fasta_raw = row.get("fasta_path", "")

        if genome_id == "":
            raise PanScapeUsageError(f"Missing genome_id at row {idx + 2} in {manifest_path}")
        if fasta_raw == "":
            raise PanScapeUsageError(f"Missing fasta_path at row {idx + 2} in {manifest_path}")

        fasta_path = _resolve_manifest_path(manifest_path.parent, fasta_raw)
        validate_existing_file(fasta_path, f"FASTA for genome `{genome_id}`")
        if not _matches_suffix(fasta_path, FASTA_SUFFIXES):
            raise PanScapeUsageError(
                f"Unsupported FASTA extension for genome `{genome_id}`: {fasta_path}"
            )

        records.append(
            GenomeRecord(
                genome_id=genome_id,
                fasta_path=fasta_path,
                completeness=_parse_optional_float(row.get("completeness", ""), "completeness", idx),
                contamination=_parse_optional_float(
                    row.get("contamination", ""), "contamination", idx
                ),
                sample_id=row.get("sample_id", "") or None,
            )
        )

    return records


def validate_samples_manifest(manifest_path: Path) -> list[SampleRecord]:
    validate_existing_file(manifest_path, "samples.tsv")
    header, rows = _read_manifest(manifest_path)
    _require_columns(manifest_path, header, {"sample_id", "r1"})

    if not rows:
        raise PanScapeUsageError(f"samples.tsv has no data rows: {manifest_path}")

    records: list[SampleRecord] = []
    for idx, row in enumerate(rows):
        sample_id = row.get("sample_id", "")
        r1_raw = row.get("r1", "")
        r2_raw = row.get("r2", "")

        if sample_id == "":
            raise PanScapeUsageError(f"Missing sample_id at row {idx + 2} in {manifest_path}")
        if r1_raw == "":
            raise PanScapeUsageError(f"Missing r1 at row {idx + 2} in {manifest_path}")

        r1 = _resolve_manifest_path(manifest_path.parent, r1_raw)
        validate_existing_file(r1, f"R1 for sample `{sample_id}`")
        if not _matches_suffix(r1, FASTQ_SUFFIXES):
            raise PanScapeUsageError(f"Unsupported FASTQ extension for sample `{sample_id}`: {r1}")

        r2 = None
        if r2_raw:
            r2 = _resolve_manifest_path(manifest_path.parent, r2_raw)
            validate_existing_file(r2, f"R2 for sample `{sample_id}`")
            if not _matches_suffix(r2, FASTQ_SUFFIXES):
                raise PanScapeUsageError(
                    f"Unsupported FASTQ extension for sample `{sample_id}`: {r2}"
                )

        records.append(SampleRecord(sample_id=sample_id, r1=r1, r2=r2))

    return records
