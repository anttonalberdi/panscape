from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from panscape.exceptions import PanScapeUsageError


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _check_overwrite(path: Path, force: bool) -> None:
    if path.exists() and not force:
        raise PanScapeUsageError(f"Refusing to overwrite existing file without --force: {path}")


def write_text(path: Path, content: str, *, force: bool = False) -> Path:
    ensure_dir(path.parent)
    _check_overwrite(path, force)
    path.write_text(content, encoding="utf-8")
    return path


def write_json(path: Path, payload: Mapping[str, Any], *, force: bool = False) -> Path:
    ensure_dir(path.parent)
    _check_overwrite(path, force)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True), encoding="utf-8")
    return path


def write_tsv(
    path: Path,
    header: Sequence[str],
    rows: Iterable[Sequence[Any]],
    *,
    force: bool = False,
) -> Path:
    ensure_dir(path.parent)
    _check_overwrite(path, force)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        for row in rows:
            writer.writerow(list(row))

    return path


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise PanScapeUsageError(f"TSV file has no header row: {path}")
        return [dict(row) for row in reader]
