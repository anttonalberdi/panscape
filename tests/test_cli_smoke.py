from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from panscape.cli import app

runner = CliRunner()


def test_root_help_smoke() -> None:
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "build" in result.stdout
    assert "map" in result.stdout
    assert "species" in result.stdout
    assert "strain" in result.stdout


def test_build_dry_run_writes_manifest(tmp_path: Path) -> None:
    fasta_path = tmp_path / "g1.fna"
    fasta_path.write_text(">g1\nACGTACGT\n", encoding="utf-8")

    genomes_tsv = tmp_path / "genomes.tsv"
    genomes_tsv.write_text("genome_id\tfasta_path\ng1\tg1.fna\n", encoding="utf-8")

    outdir = tmp_path / "run_build"
    result = runner.invoke(
        app,
        [
            "build",
            "--mock",
            "--genomes-tsv",
            str(genomes_tsv),
            "--outdir",
            str(outdir),
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    manifest_path = outdir / "panscape_manifest.json"
    assert manifest_path.exists()

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["command"] == "build"
    assert manifest["status"] == "dry-run"
    assert not (outdir / "build" / "species_clusters.tsv").exists()


def test_map_dry_run_writes_manifest(tmp_path: Path) -> None:
    r1 = tmp_path / "sample1_R1.fastq"
    r1.write_text("@r1\nACGT\n+\n####\n", encoding="utf-8")

    samples_tsv = tmp_path / "samples.tsv"
    samples_tsv.write_text("sample_id\tr1\ns1\tsample1_R1.fastq\n", encoding="utf-8")

    outdir = tmp_path / "run_map"
    result = runner.invoke(
        app,
        [
            "map",
            "--samples-tsv",
            str(samples_tsv),
            "--outdir",
            str(outdir),
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    manifest_path = outdir / "panscape_manifest.json"
    assert manifest_path.exists()

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["command"] == "map"
    assert manifest["status"] == "dry-run"
    assert not list((outdir / "map").glob("**/*.bam"))
