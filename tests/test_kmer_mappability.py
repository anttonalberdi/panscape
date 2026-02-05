from __future__ import annotations

from panscape.build.kmer import unique_kmer_fractions


def test_unique_kmer_fraction() -> None:
    sequences = {
        "gene1": "AAATTT",
        "gene2": "CCGCGG",
        "gene3": "AAACCC",
    }

    scores = unique_kmer_fractions(sequences, k=3)

    assert scores["gene2"].unique_fraction == 1.0
    assert scores["gene1"].unique_fraction < 1.0
    assert scores["gene3"].unique_fraction < 1.0
