from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Mapping


@dataclass(frozen=True, slots=True)
class KmerUniqueness:
    total_kmers: int
    unique_kmers: int
    unique_fraction: float


def _stable_hash64(value: str) -> int:
    digest = hashlib.blake2b(value.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, byteorder="big", signed=False)


def iter_kmers(sequence: str, k: int) -> list[str]:
    normalized = sequence.upper()
    if k <= 0 or len(normalized) < k:
        return []
    return [normalized[idx : idx + k] for idx in range(0, len(normalized) - k + 1)]


def kmer_set(sequence: str, k: int) -> set[str]:
    return set(iter_kmers(sequence, k))


def jaccard_similarity(left: set[str], right: set[str]) -> float:
    if not left and not right:
        return 1.0
    union = left | right
    if not union:
        return 1.0
    return len(left & right) / len(union)


def minhash_sketch(kmers: set[str], sketch_size: int) -> set[str]:
    if sketch_size <= 0:
        return set()
    hashed = sorted((( _stable_hash64(kmer), kmer) for kmer in kmers), key=lambda item: item[0])
    return {kmer for _, kmer in hashed[:sketch_size]}


def mock_mash_distance(
    seq_left: str,
    seq_right: str,
    *,
    k: int,
    sketch_size: int,
) -> float:
    """Compute a deterministic Mash-distance approximation from minhash sketches."""

    left_sketch = minhash_sketch(kmer_set(seq_left, k), sketch_size)
    right_sketch = minhash_sketch(kmer_set(seq_right, k), sketch_size)
    similarity = jaccard_similarity(left_sketch, right_sketch)
    if similarity <= 0.0:
        return 1.0
    mash_distance = (-1.0 / float(k)) * math.log((2.0 * similarity) / (1.0 + similarity))
    return max(0.0, min(1.0, mash_distance))


def unique_kmer_fractions(sequences: Mapping[str, str], k: int) -> dict[str, KmerUniqueness]:
    """Calculate per-sequence fraction of k-mers unique within a catalog."""

    per_gene_kmers: dict[str, list[str]] = {}
    global_counts: dict[str, int] = {}

    for gene_id in sorted(sequences):
        kmers = iter_kmers(sequences[gene_id], k)
        per_gene_kmers[gene_id] = kmers
        for kmer in kmers:
            global_counts[kmer] = global_counts.get(kmer, 0) + 1

    results: dict[str, KmerUniqueness] = {}
    for gene_id in sorted(per_gene_kmers):
        kmers = per_gene_kmers[gene_id]
        total = len(kmers)
        unique = sum(1 for kmer in kmers if global_counts[kmer] == 1)
        fraction = float(unique / total) if total > 0 else 0.0
        results[gene_id] = KmerUniqueness(
            total_kmers=total,
            unique_kmers=unique,
            unique_fraction=fraction,
        )

    return results
