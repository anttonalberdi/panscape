from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from panscape.exceptions import PanScapeUsageError
from panscape.utils.io import ensure_dir


@dataclass(frozen=True, slots=True)
class FastaRecord:
    """Simple FASTA record."""

    header: str
    sequence: str


@dataclass(frozen=True, slots=True)
class AssemblyStats:
    """Assembly summary statistics for one genome FASTA."""

    genome_size: int
    contig_count: int
    n50: int


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _write_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "wt", encoding="utf-8")
    return path.open("w", encoding="utf-8")


def read_fasta_records(path: Path) -> list[FastaRecord]:
    """Read FASTA records, preserving order."""

    records: list[FastaRecord] = []
    header: str | None = None
    seq_chunks: list[str] = []

    with _open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(seq_chunks).upper()
                    records.append(FastaRecord(header=header, sequence=sequence))
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)

    if header is not None:
        sequence = "".join(seq_chunks).upper()
        records.append(FastaRecord(header=header, sequence=sequence))

    return records


def write_fasta_records(
    path: Path,
    records: Iterable[FastaRecord],
    line_width: int = 80,
    *,
    force: bool = False,
) -> Path:
    """Write FASTA records using deterministic line wrapping."""

    ensure_dir(path.parent)
    if path.exists() and not force:
        raise PanScapeUsageError(f"Refusing to overwrite existing file without --force: {path}")
    with _write_text(path) as handle:
        for record in records:
            handle.write(f">{record.header}\n")
            sequence = record.sequence
            for start in range(0, len(sequence), line_width):
                handle.write(f"{sequence[start:start + line_width]}\n")

    return path


def compute_assembly_stats(records: Iterable[FastaRecord]) -> AssemblyStats:
    """Compute genome size, contig count, and N50 from FASTA records."""

    lengths = sorted((len(record.sequence) for record in records), reverse=True)
    contig_count = len(lengths)
    genome_size = sum(lengths)

    if not lengths:
        return AssemblyStats(genome_size=0, contig_count=0, n50=0)

    half = genome_size / 2
    cumulative = 0
    n50 = 0
    for length in lengths:
        cumulative += length
        if cumulative >= half:
            n50 = length
            break

    return AssemblyStats(genome_size=genome_size, contig_count=contig_count, n50=n50)


def normalize_fasta_records(genome_id: str, records: Iterable[FastaRecord]) -> list[FastaRecord]:
    """Rename contigs deterministically to stable IDs."""

    normalized: list[FastaRecord] = []
    for idx, record in enumerate(records, start=1):
        normalized.append(
            FastaRecord(
                header=f"{genome_id}_contig_{idx:06d}",
                sequence=record.sequence.upper(),
            )
        )
    return normalized


def normalize_fasta_file(
    genome_id: str,
    input_fasta: Path,
    output_fasta: Path,
    *,
    force: bool = False,
) -> tuple[list[FastaRecord], AssemblyStats]:
    """Normalize contig names and write a canonical FASTA file."""

    records = read_fasta_records(input_fasta)
    normalized = normalize_fasta_records(genome_id, records)
    write_fasta_records(output_fasta, normalized, force=force)
    return normalized, compute_assembly_stats(normalized)
