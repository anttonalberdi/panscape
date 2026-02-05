from __future__ import annotations

from panscape.build.fasta import FastaRecord
from panscape.build.representation import (
    build_consensus_dereplicated_records,
    select_medoid_genome,
)


def test_select_medoid_genome_uses_mean_ani_then_quality() -> None:
    genome_ids = ["g1", "g2", "g3"]
    ani_by_pair = {
        ("g1", "g2"): 0.99,
        ("g1", "g3"): 0.96,
        ("g2", "g3"): 0.95,
    }
    quality_scores = {"g1": 80.0, "g2": 90.0, "g3": 85.0}

    medoid = select_medoid_genome(
        genome_ids,
        ani_by_pair=ani_by_pair,
        quality_scores=quality_scores,
    )

    # g1 mean ANI = (0.99 + 0.96)/2 = 0.975, highest among the three.
    assert medoid == "g1"


def test_build_consensus_dereplicates_identical_contigs() -> None:
    records_by_genome = {
        "g1": [
            FastaRecord(header="g1_c1", sequence="AAAA"),
            FastaRecord(header="g1_c2", sequence="CCCC"),
        ],
        "g2": [FastaRecord(header="g2_c1", sequence="AAAA")],
    }

    consensus_records = build_consensus_dereplicated_records(
        records_by_genome=records_by_genome,
        genome_ids=["g1", "g2"],
        quality_scores={"g1": 100.0, "g2": 90.0},
        anchor_genome_id="g1",
        derep_identity=1.0,
        k=2,
    )

    assert len(consensus_records) == 2
    assert [record.sequence for record in consensus_records] == ["AAAA", "CCCC"]


def test_build_consensus_majority_vote_for_cluster() -> None:
    records_by_genome = {
        "g1": [FastaRecord(header="g1_c1", sequence="ACGT")],
        "g2": [FastaRecord(header="g2_c1", sequence="ACGT")],
        "g3": [FastaRecord(header="g3_c1", sequence="ACGA")],
    }

    consensus_records = build_consensus_dereplicated_records(
        records_by_genome=records_by_genome,
        genome_ids=["g1", "g2", "g3"],
        quality_scores={"g1": 100.0, "g2": 90.0, "g3": 80.0},
        anchor_genome_id="g1",
        derep_identity=0.40,
        k=2,
    )

    assert len(consensus_records) == 1
    assert consensus_records[0].sequence == "ACGT"
