from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from panscape.build.fasta import FastaRecord
from panscape.exceptions import PanScapeUsageError
from panscape.utils.io import write_text


@dataclass(frozen=True, slots=True)
class CalledGene:
    gene_id: str
    genome_id: str
    contig_id: str
    start: int
    end: int
    strand: str
    nt_sequence: str
    aa_sequence: str


CODON_TABLE: dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1]


def translate_dna(sequence: str) -> str:
    seq = sequence.upper()
    amino_acids: list[str] = []
    for idx in range(0, len(seq) - 2, 3):
        codon = seq[idx : idx + 3]
        amino_acids.append(CODON_TABLE.get(codon, "X"))
    return "".join(amino_acids).rstrip("*")


def _call_genes_mock(
    *,
    genome_id: str,
    records: Iterable[FastaRecord],
    min_gene_len: int,
) -> list[CalledGene]:
    genes: list[CalledGene] = []
    counter = 0

    for record in records:
        sequence = record.sequence.upper()
        if len(sequence) < min_gene_len:
            continue

        for start0 in range(0, len(sequence) - min_gene_len + 1, min_gene_len):
            end0 = start0 + min_gene_len
            nt_sequence = sequence[start0:end0]
            if len(nt_sequence) < min_gene_len:
                continue

            counter += 1
            gene_id = f"{genome_id}_gene_{counter:06d}"
            genes.append(
                CalledGene(
                    gene_id=gene_id,
                    genome_id=genome_id,
                    contig_id=record.header,
                    start=start0 + 1,
                    end=end0,
                    strand="+",
                    nt_sequence=nt_sequence,
                    aa_sequence=translate_dna(nt_sequence),
                )
            )

    return genes


def _call_genes_pyrodigal(
    *,
    genome_id: str,
    records: Iterable[FastaRecord],
    min_gene_len: int | None,
) -> list[CalledGene]:
    try:
        import pyrodigal
    except ImportError as exc:  # pragma: no cover - pyrodigal not required for mock tests
        raise PanScapeUsageError(
            "pyrodigal is required for non-mock gene calling. Install the `pyrodigal` package or use --mock."
        ) from exc

    finder = pyrodigal.GeneFinder(meta=True)
    genes: list[CalledGene] = []
    counter = 0

    for record in records:
        try:
            predictions = finder.find_genes(record.sequence)
        except TypeError:
            predictions = finder.find_genes(record.sequence.encode("ascii"))

        for prediction in predictions:
            begin = int(getattr(prediction, "begin", getattr(prediction, "start", 0)))
            end = int(getattr(prediction, "end", getattr(prediction, "stop", 0)))
            strand_value = int(getattr(prediction, "strand", 1))
            strand = "+" if strand_value >= 0 else "-"

            if begin <= 0 or end <= 0 or end < begin:
                continue

            nt_sequence = record.sequence[begin - 1 : end]
            if strand == "-":
                nt_sequence = reverse_complement(nt_sequence)
            if min_gene_len is not None and len(nt_sequence) < min_gene_len:
                continue

            if hasattr(prediction, "translate"):
                aa_sequence = str(prediction.translate())
            else:
                aa_sequence = translate_dna(nt_sequence)

            counter += 1
            gene_id = f"{genome_id}_gene_{counter:06d}"
            genes.append(
                CalledGene(
                    gene_id=gene_id,
                    genome_id=genome_id,
                    contig_id=record.header,
                    start=begin,
                    end=end,
                    strand=strand,
                    nt_sequence=nt_sequence,
                    aa_sequence=aa_sequence,
                )
            )

    return genes


def call_genes(
    *,
    genome_id: str,
    records: Iterable[FastaRecord],
    min_gene_len: int,
    mock: bool,
    apply_min_gene_len: bool = True,
) -> list[CalledGene]:
    """Call genes from normalized FASTA records."""

    if mock:
        return _call_genes_mock(genome_id=genome_id, records=records, min_gene_len=min_gene_len)

    resolved_min_len: int | None = min_gene_len if apply_min_gene_len else None
    return _call_genes_pyrodigal(genome_id=genome_id, records=records, min_gene_len=resolved_min_len)


def write_gene_outputs(
    *,
    genome_id: str,
    genes: Iterable[CalledGene],
    genes_dir: Path,
    force: bool,
) -> tuple[Path, Path, Path]:
    """Write per-genome FFN/FAA/GFF files."""

    genes_sorted = sorted(genes, key=lambda gene: (gene.contig_id, gene.start, gene.gene_id))

    ffn_lines: list[str] = []
    faa_lines: list[str] = []
    gff_lines: list[str] = ["##gff-version 3"]

    for gene in genes_sorted:
        ffn_lines.append(f">{gene.gene_id}")
        ffn_lines.append(gene.nt_sequence)

        faa_lines.append(f">{gene.gene_id}")
        faa_lines.append(gene.aa_sequence)

        attributes = f"ID={gene.gene_id};genome_id={gene.genome_id}"
        gff_lines.append(
            "\t".join(
                [
                    gene.contig_id,
                    "PanScape",
                    "CDS",
                    str(gene.start),
                    str(gene.end),
                    ".",
                    gene.strand,
                    "0",
                    attributes,
                ]
            )
        )

    if len(ffn_lines) == 0:
        ffn_lines.append("# No genes were called for this genome.")
    if len(faa_lines) == 0:
        faa_lines.append("# No genes were called for this genome.")

    ffn_path = write_text(genes_dir / f"{genome_id}.ffn", "\n".join(ffn_lines) + "\n", force=force)
    faa_path = write_text(genes_dir / f"{genome_id}.faa", "\n".join(faa_lines) + "\n", force=force)
    gff_path = write_text(genes_dir / f"{genome_id}.gff", "\n".join(gff_lines) + "\n", force=force)

    return ffn_path, faa_path, gff_path
