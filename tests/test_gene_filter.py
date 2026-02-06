from __future__ import annotations

from panscape.build.gene_calling import CalledGene
from panscape.build.gene_filter import filter_genes_by_length


def test_filter_genes_by_length() -> None:
    genes = [
        CalledGene(
            gene_id="g1",
            genome_id="genome1",
            contig_id="contig1",
            start=1,
            end=100,
            strand="+",
            nt_sequence="A" * 50,
            aa_sequence="M" * 16,
        ),
        CalledGene(
            gene_id="g2",
            genome_id="genome1",
            contig_id="contig1",
            start=200,
            end=600,
            strand="+",
            nt_sequence="A" * 300,
            aa_sequence="M" * 100,
        ),
    ]

    filtered = filter_genes_by_length(genes, min_gene_len=100)
    assert [gene.gene_id for gene in filtered] == ["g2"]
