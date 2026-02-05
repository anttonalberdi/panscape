from __future__ import annotations

from panscape.build.scoring import backbone_score


def test_backbone_score_formula() -> None:
    score = backbone_score(
        completeness=95.0,
        contamination=2.0,
        n50=50000,
        genome_size=5000000,
        contig_count=100,
    )

    expected = 95.0 - (5.0 * 2.0) + (0.001 * 50000) - (0.000001 * 5000000) - (0.01 * 100)
    assert abs(score - expected) < 1e-9


def test_backbone_score_uses_neutral_values_when_qc_missing() -> None:
    score_missing = backbone_score(
        completeness=None,
        contamination=None,
        n50=10000,
        genome_size=2000000,
        contig_count=20,
    )
    score_explicit = backbone_score(
        completeness=85.0,
        contamination=5.0,
        n50=10000,
        genome_size=2000000,
        contig_count=20,
    )

    assert abs(score_missing - score_explicit) < 1e-9
