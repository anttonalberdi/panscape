from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from panscape.cli import app

runner = CliRunner()


def _assert_nonempty_tsv(path: Path) -> None:
    content = path.read_text(encoding="utf-8").strip().splitlines()
    assert len(content) >= 2


def test_build_mock_integration(tmp_path: Path) -> None:
    data_dir = Path(__file__).parent / "data"
    genomes_tsv = data_dir / "genomes.tsv"
    outdir = tmp_path / "build_run"

    result = runner.invoke(
        app,
        [
            "build",
            "--mock",
            "--genomes-tsv",
            str(genomes_tsv),
            "--outdir",
            str(outdir),
            "--threads",
            "2",
        ],
    )

    assert result.exit_code == 0, result.stdout

    manifest_path = outdir / "panscape_manifest.json"
    assert manifest_path.exists()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["command"] == "build"
    assert manifest["status"] == "completed"
    assert manifest["tool_versions"]["mash"] == "mock"

    _assert_nonempty_tsv(outdir / "build" / "qc" / "assembly_stats.tsv")
    _assert_nonempty_tsv(outdir / "build" / "qc" / "genome_qc.tsv")
    _assert_nonempty_tsv(outdir / "build" / "mash" / "preclusters.tsv")
    _assert_nonempty_tsv(outdir / "build" / "ani" / "skani_edges_species.tsv")
    _assert_nonempty_tsv(outdir / "build" / "clusters" / "species_clusters.tsv")
    _assert_nonempty_tsv(outdir / "build" / "clusters" / "strain_bins.tsv")
    _assert_nonempty_tsv(outdir / "build" / "backbones" / "backbone_table.tsv")

    for genome_id in ["g1", "g2", "g3"]:
        assert (outdir / "build" / "genes" / f"{genome_id}.ffn").exists()
        assert (outdir / "build" / "genes" / f"{genome_id}.faa").exists()
        assert (outdir / "build" / "genes" / f"{genome_id}.gff").exists()

    species_dirs = sorted(path for path in (outdir / "build" / "pangenomes").iterdir() if path.is_dir())
    assert species_dirs

    for species_dir in species_dirs:
        assert (species_dir / "representatives.faa").exists()
        _assert_nonempty_tsv(species_dir / "gene_families.tsv")
        _assert_nonempty_tsv(species_dir / "family_prevalence.tsv")
        assert (species_dir / "genes.fna").exists()
        assert (species_dir / "genes.faa").exists()
        _assert_nonempty_tsv(species_dir / "gene_representatives.tsv")
        _assert_nonempty_tsv(species_dir / "gene_mappability.tsv")
        _assert_nonempty_tsv(species_dir / "gene_tiers.tsv")
        _assert_nonempty_tsv(species_dir / "catalog_summary.tsv")
        assert (species_dir / "index" / "INDEX_PLAN.json").exists()
