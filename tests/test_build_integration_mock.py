from __future__ import annotations

import csv
import json
from pathlib import Path

from typer.testing import CliRunner

from panscape.cli import app
from panscape.commands.build import (
    _looks_like_checkm2_multiprocessing_bug,
    detect_fasta_extension,
    resolve_checkm2_db_path,
)

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


def test_build_mock_defaults_to_checkm2_when_qc_missing(tmp_path: Path) -> None:
    g1 = tmp_path / "g1.fna"
    g2 = tmp_path / "g2.fna"
    g1.write_text(">g1\nACGTACGTACGTACGTACGTACGT\n", encoding="utf-8")
    g2.write_text(">g2\nACGTACGTACGTACGTACGTACGA\n", encoding="utf-8")

    outdir = tmp_path / "build_run_no_qc"
    result = runner.invoke(
        app,
        [
            "build",
            "--mock",
            "--genomes-files",
            f"{g1},{g2}",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.stdout

    qc_path = outdir / "build" / "qc" / "genome_qc.tsv"
    assert qc_path.exists()
    with qc_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert rows
    assert all(row["qc_source"] == "checkm2" for row in rows)


def test_build_mock_accepts_checkm2_db_directory_or_dmnd_file(tmp_path: Path) -> None:
    g1 = tmp_path / "g1.fna"
    g2 = tmp_path / "g2.fna"
    g1.write_text(">g1\nACGTACGTACGTACGTACGTACGT\n", encoding="utf-8")
    g2.write_text(">g2\nACGTACGTACGTACGTACGTACGA\n", encoding="utf-8")

    db_dir = tmp_path / "checkm2_db"
    db_dir.mkdir()
    dmnd_file = db_dir / "uniref100.KO.1.dmnd"
    dmnd_file.write_text("placeholder", encoding="utf-8")

    outdir_dir_mode = tmp_path / "build_run_db_dir"
    result_dir_mode = runner.invoke(
        app,
        [
            "build",
            "--mock",
            "--genomes-files",
            f"{g1},{g2}",
            "--checkm2-db",
            str(db_dir),
            "--outdir",
            str(outdir_dir_mode),
        ],
    )
    assert result_dir_mode.exit_code == 0, result_dir_mode.stdout

    outdir_file_mode = tmp_path / "build_run_db_file"
    result_file_mode = runner.invoke(
        app,
        [
            "build",
            "--mock",
            "--genomes-files",
            f"{g1},{g2}",
            "--checkm2-db",
            str(dmnd_file),
            "--outdir",
            str(outdir_file_mode),
        ],
    )
    assert result_file_mode.exit_code == 0, result_file_mode.stdout


def test_resolve_checkm2_db_path_returns_dmnd_file_for_both_inputs(tmp_path: Path) -> None:
    db_dir = tmp_path / "checkm2_db"
    db_dir.mkdir()
    dmnd_file = db_dir / "uniref100.KO.1.dmnd"
    dmnd_file.write_text("placeholder", encoding="utf-8")

    assert resolve_checkm2_db_path(db_dir) == dmnd_file
    assert resolve_checkm2_db_path(dmnd_file) == dmnd_file


def test_detect_checkm2_multiprocessing_bug_patterns() -> None:
    assert _looks_like_checkm2_multiprocessing_bug(
        "AttributeError: 'Predictor' object has no attribute '__set_up_prodigal_thread'"
    )
    assert _looks_like_checkm2_multiprocessing_bug(
        "AttributeError: 'Predictor' object has no attribute '__reportProgress'"
    )
    assert _looks_like_checkm2_multiprocessing_bug(
        "ERROR: No protein files were generated in build/qc/checkm2/protein_files"
    )
    assert not _looks_like_checkm2_multiprocessing_bug("some other command failure")


def test_detect_fasta_extension_for_checkm2() -> None:
    assert detect_fasta_extension([Path("a.fna"), Path("b.fna.gz")]) == "fna"
    assert detect_fasta_extension([Path("a.fa.gz"), Path("b.fa.gz"), Path("c.fna")]) == "fa.gz"


def test_build_mock_consensus_representation_outputs_consensus_backbones(tmp_path: Path) -> None:
    data_dir = Path(__file__).parent / "data"
    genomes_tsv = data_dir / "genomes.tsv"
    outdir = tmp_path / "build_run_consensus"

    result = runner.invoke(
        app,
        [
            "build",
            "--mock",
            "--species-representation",
            "consensus",
            "--genomes-tsv",
            str(genomes_tsv),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.stdout

    backbone_table_path = outdir / "build" / "backbones" / "backbone_table.tsv"
    with backbone_table_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert rows
    for row in rows:
        assert row["representation_mode"] == "consensus"
        assert row["source_genome_id"] == ""
        species_id = row["species_id"]
        backbone_fasta = outdir / "build" / "backbones" / f"{species_id}.fna"
        assert backbone_fasta.exists()
        first_header = backbone_fasta.read_text(encoding="utf-8").splitlines()[0]
        assert first_header.startswith(f">{species_id}_contig_")

        index_plan = json.loads(
            (outdir / "build" / "pangenomes" / species_id / "index" / "INDEX_PLAN.json").read_text(
                encoding="utf-8"
            )
        )
        assert index_plan["representation_mode"] == "consensus"
        assert index_plan["anchor_genome_id"] == row["backbone_genome_id"]
        assert index_plan["source_genome_id"] is None


def test_build_mock_medoid_representation_sets_source_genome(tmp_path: Path) -> None:
    data_dir = Path(__file__).parent / "data"
    genomes_tsv = data_dir / "genomes.tsv"
    outdir = tmp_path / "build_run_medoid"

    result = runner.invoke(
        app,
        [
            "build",
            "--mock",
            "--species-representation",
            "medoid",
            "--genomes-tsv",
            str(genomes_tsv),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.stdout

    with (outdir / "build" / "backbones" / "backbone_table.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    assert all(row["representation_mode"] == "medoid" for row in rows)
    assert all(row["source_genome_id"] != "" for row in rows)
