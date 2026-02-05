from __future__ import annotations

from dataclasses import dataclass

NEUTRAL_COMPLETENESS = 85.0
NEUTRAL_CONTAMINATION = 5.0


@dataclass(frozen=True, slots=True)
class QCMetrics:
    completeness: float | None
    contamination: float | None


def backbone_score(
    *,
    completeness: float | None,
    contamination: float | None,
    n50: int,
    genome_size: int,
    contig_count: int,
    neutral_completeness: float = NEUTRAL_COMPLETENESS,
    neutral_contamination: float = NEUTRAL_CONTAMINATION,
) -> float:
    """Deterministic backbone score used for representative genome selection."""

    resolved_completeness = completeness if completeness is not None else neutral_completeness
    resolved_contamination = contamination if contamination is not None else neutral_contamination

    return (
        resolved_completeness
        - (5.0 * resolved_contamination)
        + (0.001 * float(n50))
        - (0.000001 * float(genome_size))
        - (0.01 * float(contig_count))
    )
