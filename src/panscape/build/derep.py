from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

from panscape.build.kmer import jaccard_similarity, kmer_set


@dataclass(frozen=True, slots=True)
class DerepGene:
    gene_id: str
    genome_id: str
    family_id: str
    nt_sequence: str
    aa_sequence: str


@dataclass(frozen=True, slots=True)
class DerepGroup:
    family_id: str
    group_id: str
    representative_gene_id: str
    member_gene_ids: tuple[str, ...]


def approximate_nt_identity(seq_left: str, seq_right: str, *, k: int) -> float:
    """Approximate nucleotide identity using direct or k-mer comparisons."""

    left = seq_left.upper()
    right = seq_right.upper()

    if left == right:
        return 1.0

    if len(left) == len(right) and len(left) < k:
        if len(left) == 0:
            return 0.0
        matches = sum(1 for a, b in zip(left, right, strict=True) if a == b)
        return matches / len(left)

    left_kmers = kmer_set(left, k)
    right_kmers = kmer_set(right, k)
    return jaccard_similarity(left_kmers, right_kmers)


def dereplicate_family(
    members: Sequence[DerepGene],
    *,
    identity_threshold: float,
    k: int,
    preference_key: Callable[[DerepGene], tuple],
) -> list[DerepGroup]:
    """Collapse near-identical genes in one family using deterministic greedy grouping."""

    if not members:
        return []

    sorted_members = sorted(members, key=lambda member: (preference_key(member), member.gene_id))

    group_representatives: list[DerepGene] = []
    group_members: list[list[DerepGene]] = []

    for member in sorted_members:
        assigned = False
        for idx, representative in enumerate(group_representatives):
            identity = approximate_nt_identity(member.nt_sequence, representative.nt_sequence, k=k)
            if identity >= identity_threshold:
                group_members[idx].append(member)
                if preference_key(member) < preference_key(representative):
                    group_representatives[idx] = member
                assigned = True
                break

        if not assigned:
            group_representatives.append(member)
            group_members.append([member])

    ordered_indices = sorted(
        range(len(group_representatives)),
        key=lambda idx: (group_representatives[idx].gene_id, tuple(g.gene_id for g in group_members[idx])),
    )

    groups: list[DerepGroup] = []
    for order_idx, group_idx in enumerate(ordered_indices, start=1):
        representative = group_representatives[group_idx]
        members_sorted = tuple(sorted(gene.gene_id for gene in group_members[group_idx]))
        groups.append(
            DerepGroup(
                family_id=representative.family_id,
                group_id=f"group_{order_idx:03d}",
                representative_gene_id=representative.gene_id,
                member_gene_ids=members_sorted,
            )
        )

    groups.sort(key=lambda group: group.group_id)
    return groups
