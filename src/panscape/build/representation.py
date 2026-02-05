from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

from panscape.build.fasta import FastaRecord
from panscape.build.kmer import jaccard_similarity, kmer_set


@dataclass(frozen=True, slots=True)
class ContigEntry:
    genome_id: str
    contig_id: str
    sequence: str


def _pair_key(left: str, right: str) -> tuple[str, str]:
    return (left, right) if left <= right else (right, left)


def select_medoid_genome(
    genome_ids: Iterable[str],
    *,
    ani_by_pair: Mapping[tuple[str, str], float],
    quality_scores: Mapping[str, float],
) -> str:
    """Select the genome with highest mean ANI within the component."""

    nodes = sorted(set(genome_ids))
    if len(nodes) == 0:
        raise ValueError("Cannot select medoid from an empty genome set.")
    if len(nodes) == 1:
        return nodes[0]

    mean_ani_by_genome: dict[str, float] = {}
    for genome_id in nodes:
        totals = 0.0
        count = 0
        for other_id in nodes:
            if other_id == genome_id:
                continue
            pair = _pair_key(genome_id, other_id)
            ani = ani_by_pair.get(pair, 0.0)
            totals += ani
            count += 1
        mean_ani_by_genome[genome_id] = (totals / float(count)) if count else 0.0

    # Deterministic tie-breaking by quality score then genome_id.
    return sorted(
        nodes,
        key=lambda genome_id: (
            -mean_ani_by_genome[genome_id],
            -quality_scores.get(genome_id, 0.0),
            genome_id,
        ),
    )[0]


def _short_sequence_identity(left: str, right: str) -> float:
    if not left and not right:
        return 1.0
    if not left or not right:
        return 0.0

    overlap = min(len(left), len(right))
    if overlap == 0:
        return 0.0

    matches = sum(1 for idx in range(overlap) if left[idx] == right[idx])
    return matches / float(max(len(left), len(right)))


def _sequence_similarity(
    left: str,
    right: str,
    *,
    k: int,
    left_kmers: set[str] | None = None,
    right_kmers: set[str] | None = None,
) -> float:
    if left == right:
        return 1.0
    if min(len(left), len(right)) < k:
        return _short_sequence_identity(left, right)

    left_resolved = left_kmers if left_kmers is not None else kmer_set(left, k)
    right_resolved = right_kmers if right_kmers is not None else kmer_set(right, k)
    return jaccard_similarity(left_resolved, right_resolved)


def _consensus_sequence(sequences: Sequence[str]) -> str:
    if len(sequences) == 0:
        return "N"
    if len(sequences) == 1:
        return sequences[0]

    max_len = max(len(sequence) for sequence in sequences)
    consensus: list[str] = []

    for idx in range(max_len):
        counts: dict[str, int] = {}
        for sequence in sequences:
            if idx >= len(sequence):
                continue
            base = sequence[idx].upper()
            if base not in {"A", "C", "G", "T", "N"}:
                base = "N"
            counts[base] = counts.get(base, 0) + 1

        if not counts:
            continue

        consensus_base = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]
        consensus.append(consensus_base)

    return "".join(consensus) if consensus else "N"


@dataclass(slots=True)
class _ClusterState:
    representative: ContigEntry
    representative_kmers: set[str]
    members: list[ContigEntry]


def build_consensus_dereplicated_records(
    *,
    records_by_genome: Mapping[str, Sequence[FastaRecord]],
    genome_ids: Iterable[str],
    quality_scores: Mapping[str, float],
    anchor_genome_id: str,
    derep_identity: float = 0.995,
    k: int = 21,
) -> list[FastaRecord]:
    """Build a species consensus by dereplicating near-identical contigs first.

    The implementation is intentionally deterministic and lightweight:
    it clusters contigs by sequence similarity proxy (k-mer Jaccard for long
    contigs and positional identity for short contigs), then emits one consensus
    sequence per cluster.
    """

    entries: list[ContigEntry] = []
    for genome_id in sorted(set(genome_ids)):
        for record in records_by_genome.get(genome_id, []):
            entries.append(
                ContigEntry(
                    genome_id=genome_id,
                    contig_id=record.header,
                    sequence=record.sequence.upper(),
                )
            )

    def _entry_sort_key(entry: ContigEntry) -> tuple[int, float, int, str, str]:
        return (
            0 if entry.genome_id == anchor_genome_id else 1,
            -quality_scores.get(entry.genome_id, 0.0),
            -len(entry.sequence),
            entry.genome_id,
            entry.contig_id,
        )

    clusters: list[_ClusterState] = []
    for entry in sorted(entries, key=_entry_sort_key):
        entry_kmers = kmer_set(entry.sequence, k)
        assigned = False
        for cluster in clusters:
            similarity = _sequence_similarity(
                entry.sequence,
                cluster.representative.sequence,
                k=k,
                left_kmers=entry_kmers,
                right_kmers=cluster.representative_kmers,
            )
            if similarity < derep_identity:
                continue

            cluster.members.append(entry)
            if _entry_sort_key(entry) < _entry_sort_key(cluster.representative):
                cluster.representative = entry
                cluster.representative_kmers = entry_kmers
            assigned = True
            break

        if not assigned:
            clusters.append(
                _ClusterState(
                    representative=entry,
                    representative_kmers=entry_kmers,
                    members=[entry],
                )
            )

    output_records: list[FastaRecord] = []
    for idx, cluster in enumerate(sorted(clusters, key=lambda item: _entry_sort_key(item.representative)), start=1):
        consensus = _consensus_sequence([member.sequence for member in cluster.members])
        output_records.append(
            FastaRecord(
                header=f"consensus_contig_{idx:06d}",
                sequence=consensus,
            )
        )

    return output_records
