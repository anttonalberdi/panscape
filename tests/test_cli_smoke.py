from __future__ import annotations

import gzip
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
    assert "update" in result.stdout


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


def test_build_dry_run_with_genomes_files_relative_paths(tmp_path: Path, monkeypatch) -> None:
    (tmp_path / "gA.fna").write_text(">gA\nACGTACGTACGT\n", encoding="utf-8")
    (tmp_path / "gB.fna").write_text(">gB\nACGTACGTACGA\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    outdir = tmp_path / "run_build_files"
    result = runner.invoke(
        app,
        [
            "build",
            "--mock",
            "--genomes-files",
            "gA.fna,gB.fna",
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
    assert len(manifest["input_paths"]) == 2


def test_build_fails_when_genome_input_modes_are_combined(tmp_path: Path) -> None:
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()
    (genomes_dir / "gA.fna").write_text(">gA\nACGTACGTACGT\n", encoding="utf-8")
    genomes_tsv = tmp_path / "genomes.tsv"
    genomes_tsv.write_text("genome_id\tfasta_path\ngA\tgenomes/gA.fna\n", encoding="utf-8")

    outdir = tmp_path / "run_build_invalid"
    result = runner.invoke(
        app,
        [
            "build",
            "--mock",
            "--genomes-tsv",
            str(genomes_tsv),
            "--genomes-files",
            str(genomes_dir / "gA.fna"),
            "--outdir",
            str(outdir),
            "--dry-run",
        ],
    )

    assert result.exit_code != 0
    assert "mutually" in result.stdout
    assert "exclusive" in result.stdout


def test_build_dry_run_with_genomes_files_absolute_paths(tmp_path: Path) -> None:
    fasta_a = tmp_path / "g_abs_1.fna"
    fasta_b = tmp_path / "g_abs_2.fasta"
    fasta_a.write_text(">g_abs_1\nACGTACGTACGT\n", encoding="utf-8")
    fasta_b.write_text(">g_abs_2\nACGTACGTACGA\n", encoding="utf-8")

    outdir = tmp_path / "run_build_abs_files"
    result = runner.invoke(
        app,
        [
            "build",
            "--mock",
            "--genomes-files",
            f"{fasta_a},{fasta_b}",
            "--outdir",
            str(outdir),
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    manifest = json.loads((outdir / "panscape_manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "dry-run"
    assert len(manifest["input_paths"]) == 2


def test_build_dry_run_with_genomes_path_autodetects_fasta_and_gz(tmp_path: Path) -> None:
    genomes_dir = tmp_path / "autodetect"
    genomes_dir.mkdir()
    (genomes_dir / "auto1.fna").write_text(">auto1\nACGTACGTACGT\n", encoding="utf-8")
    with gzip.open(genomes_dir / "auto2.fa.gz", "wt", encoding="utf-8") as handle:
        handle.write(">auto2\nACGTACGTACGA\n")

    outdir = tmp_path / "run_build_autodetect"
    result = runner.invoke(
        app,
        [
            "build",
            "--mock",
            "--genomes-path",
            str(genomes_dir),
            "--outdir",
            str(outdir),
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    manifest = json.loads((outdir / "panscape_manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "dry-run"
    assert len(manifest["input_paths"]) == 2
    assert any(path.endswith("auto2.fa.gz") for path in manifest["input_paths"])


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


def test_update_dry_run_smoke() -> None:
    result = runner.invoke(app, ["update", "--dry-run"])
    assert result.exit_code == 0
