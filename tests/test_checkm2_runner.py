from __future__ import annotations

from pathlib import Path

from panscape.runners.checkm2 import CheckM2Runner
from panscape.utils.subprocess import CommandResult


def test_checkm2_predict_uses_directory_input_and_extension(monkeypatch) -> None:
    runner = CheckM2Runner()
    captured: dict[str, object] = {}

    def _fake_run(self, args, **kwargs):  # type: ignore[no-untyped-def]
        captured["args"] = [str(arg) for arg in args]
        captured["kwargs"] = kwargs
        return CommandResult(command=["checkm2"], returncode=0, stdout="", stderr="", dry_run=False)

    monkeypatch.setattr(CheckM2Runner, "run", _fake_run)

    runner.predict(
        input_dir=Path("in_dir"),
        output_dir=Path("out_dir"),
        threads=8,
        database_path=Path("db/uniref100.KO.1.dmnd"),
        extension="fna",
        force=True,
        dry_run=False,
    )

    args = captured["args"]
    assert isinstance(args, list)
    assert "--input" in args
    assert "in_dir" in args
    assert "--extension" in args
    assert "fna" in args
    assert "--genes" not in args


def test_checkm2_predict_supports_genes_mode_with_multiple_inputs(monkeypatch) -> None:
    runner = CheckM2Runner()
    captured: dict[str, object] = {}

    def _fake_run(self, args, **kwargs):  # type: ignore[no-untyped-def]
        captured["args"] = [str(arg) for arg in args]
        captured["kwargs"] = kwargs
        return CommandResult(command=["checkm2"], returncode=0, stdout="", stderr="", dry_run=False)

    monkeypatch.setattr(CheckM2Runner, "run", _fake_run)

    runner.predict(
        input_paths=[Path("g1.faa"), Path("g2.faa")],
        output_dir=Path("out_dir"),
        threads=4,
        database_path=Path("db/uniref100.KO.1.dmnd"),
        genes=True,
        force=True,
        dry_run=False,
    )

    args = captured["args"]
    assert isinstance(args, list)
    assert "--input" in args
    assert "g1.faa" in args
    assert "g2.faa" in args
    assert "--genes" in args
    assert "--extension" not in args
