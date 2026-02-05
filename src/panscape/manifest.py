from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Sequence

from panscape import __version__


@dataclass(slots=True)
class RunManifest:
    command: str
    argv: list[str]
    started_at: str
    ended_at: str | None
    status: str
    cwd: str
    outdir: str
    dry_run: bool
    threads: int
    config_path: str | None
    git_commit: str | None
    versions: dict[str, str]
    input_paths: list[str]
    output_paths: list[str] = field(default_factory=list)
    planned_steps: list[str] = field(default_factory=list)


def _utcnow_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _package_version(package_name: str) -> str:
    try:
        return version(package_name)
    except PackageNotFoundError:
        return "unknown"


def _detect_git_commit(cwd: Path) -> str | None:
    try:
        process = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=cwd,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return process.stdout.strip() or None


def create_run_manifest(
    *,
    command: str,
    argv: Sequence[str],
    outdir: Path,
    dry_run: bool,
    threads: int,
    config_path: Path | None,
    input_paths: Sequence[Path],
    planned_steps: Sequence[str],
) -> RunManifest:
    return RunManifest(
        command=command,
        argv=list(argv),
        started_at=_utcnow_iso(),
        ended_at=None,
        status="running",
        cwd=str(Path.cwd()),
        outdir=str(outdir),
        dry_run=dry_run,
        threads=threads,
        config_path=str(config_path) if config_path is not None else None,
        git_commit=_detect_git_commit(Path.cwd()),
        versions={
            "panscape": __version__,
            "python": sys.version.split()[0],
            "typer": _package_version("typer"),
            "pydantic": _package_version("pydantic"),
            "rich": _package_version("rich"),
        },
        input_paths=[str(path) for path in input_paths],
        planned_steps=list(planned_steps),
    )


def finalize_manifest(
    manifest: RunManifest,
    *,
    status: str,
    output_paths: Sequence[Path | str],
) -> RunManifest:
    manifest.status = status
    manifest.ended_at = _utcnow_iso()
    manifest.output_paths = [str(path) for path in output_paths]
    return manifest


def write_manifest(outdir: Path, manifest: RunManifest) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / "panscape_manifest.json"
    manifest_path.write_text(json.dumps(asdict(manifest), indent=2, ensure_ascii=True), encoding="utf-8")
    return manifest_path
