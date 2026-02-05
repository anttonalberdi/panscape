from __future__ import annotations

from panscape.build.fasta import FastaRecord, compute_assembly_stats, normalize_fasta_records


def test_compute_assembly_stats_n50() -> None:
    records = [
        FastaRecord(header="a", sequence="A" * 100),
        FastaRecord(header="b", sequence="C" * 80),
        FastaRecord(header="c", sequence="G" * 20),
    ]

    stats = compute_assembly_stats(records)

    assert stats.genome_size == 200
    assert stats.contig_count == 3
    assert stats.n50 == 100


def test_normalize_fasta_records_headers() -> None:
    records = [
        FastaRecord(header="x", sequence="AAAA"),
        FastaRecord(header="y", sequence="TTTT"),
    ]

    normalized = normalize_fasta_records("genomeA", records)

    assert [record.header for record in normalized] == [
        "genomeA_contig_000001",
        "genomeA_contig_000002",
    ]
