from __future__ import annotations

from typing import Iterable

from panscape.build.gene_calling import CalledGene


def filter_genes_by_length(
    genes: Iterable[CalledGene],
    min_gene_len: int,
) -> list[CalledGene]:
    """Filter called genes by nucleotide sequence length."""

    if min_gene_len <= 0:
        return list(genes)

    return [gene for gene in genes if len(gene.nt_sequence) >= min_gene_len]
