from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class OutputLayout:
    root: Path
    build_dir: Path
    map_dir: Path
    species_dir: Path
    strain_dir: Path


def sanitize_identifier(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned.strip("_") or "unknown"


def create_output_layout(outdir: Path) -> OutputLayout:
    root = outdir
    build_dir = root / "build"
    map_dir = root / "map"
    species_dir = root / "species"
    strain_dir = root / "strain"

    for path in (root, build_dir, map_dir, species_dir, strain_dir):
        path.mkdir(parents=True, exist_ok=True)

    return OutputLayout(
        root=root,
        build_dir=build_dir,
        map_dir=map_dir,
        species_dir=species_dir,
        strain_dir=strain_dir,
    )


def pangenome_dir(layout: OutputLayout, species_id: str) -> Path:
    return layout.build_dir / "pangenomes" / sanitize_identifier(species_id)


def map_sample_dir(layout: OutputLayout, sample_id: str) -> Path:
    return layout.map_dir / sanitize_identifier(sample_id)


def species_output_dir(layout: OutputLayout, species_id: str) -> Path:
    return layout.species_dir / sanitize_identifier(species_id)


def strain_k_dir(layout: OutputLayout, species_id: str, k_value: int) -> Path:
    return layout.strain_dir / sanitize_identifier(species_id) / f"K{k_value}"
