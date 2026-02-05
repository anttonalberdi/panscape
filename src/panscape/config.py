from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, TypeVar

import yaml
from pydantic import BaseModel, ConfigDict, Field, PositiveInt, ValidationError, model_validator

from panscape.exceptions import PanScapeUsageError


class CommonConfig(BaseModel):
    """Shared command options across PanScape subcommands."""

    model_config = ConfigDict(extra="forbid")

    outdir: Path = Field(default_factory=Path.cwd)
    threads: PositiveInt = 1
    dry_run: bool = False
    force: bool = False
    log_file: Path | None = None
    verbose: bool = False
    quiet: bool = False

    @model_validator(mode="after")
    def _validate_verbosity(self) -> "CommonConfig":
        if self.verbose and self.quiet:
            raise ValueError("`verbose` and `quiet` cannot both be true.")
        return self


class BuildConfig(CommonConfig):
    genomes_tsv: Path | None = None
    genomes_path: Path | None = None
    genomes_files: list[str] = Field(default_factory=list)
    run_checkm2: bool = False
    checkm2_db: Path | None = None

    min_completeness: float = Field(default=70.0, ge=0.0, le=100.0)
    max_contamination: float = Field(default=10.0, ge=0.0, le=100.0)

    mash_sketch_size: PositiveInt = 10000
    mash_k: PositiveInt = 21
    mash_threshold: float = Field(default=0.10, ge=0.0, le=1.0)

    species_ani: float = Field(default=0.95, ge=0.0, le=1.0)
    strain_ani: float = Field(default=0.995, ge=0.0, le=1.0)

    mmseqs_min_aa_id: float = Field(default=0.9, ge=0.0, le=1.0)
    mmseqs_cov: float = Field(default=0.8, ge=0.0, le=1.0)

    within_family_derep_nt: float = Field(default=0.995, ge=0.0, le=1.0)
    mappability_k: PositiveInt = 31
    min_gene_len: PositiveInt = 300

    rare_min_prevalence: PositiveInt = 2
    use_high_quality_for_prevalence: bool = True
    hq_min_completeness: float = Field(default=90.0, ge=0.0, le=100.0)
    hq_max_contamination: float = Field(default=5.0, ge=0.0, le=100.0)

    keep_singletons: bool = True
    mock: bool = False

    strain_ok_min_mappability: float = Field(default=0.2, ge=0.0, le=1.0)

    @model_validator(mode="after")
    def _validate_threshold_relationships(self) -> "BuildConfig":
        if self.strain_ani < self.species_ani:
            raise ValueError("`strain_ani` must be >= `species_ani`.")
        return self


class MapConfig(CommonConfig):
    samples_tsv: Path | None = None
    mapper: str = "minimap2"
    references_dir: Path | None = None


class SpeciesConfig(CommonConfig):
    map_dir: Path | None = None
    species_ids: list[str] = Field(default_factory=list)


class StrainConfig(CommonConfig):
    species_ids: list[str] = Field(default_factory=list)
    coverage_matrix: Path | None = None
    bam_dir: Path | None = None
    vcf_path: Path | None = None


class PanScapeConfig(BaseModel):
    """Top-level YAML config model."""

    model_config = ConfigDict(extra="forbid")

    build: BuildConfig | None = None
    map: MapConfig | None = None
    species: SpeciesConfig | None = None
    strain: StrainConfig | None = None


def load_config(config_path: Path | None) -> PanScapeConfig:
    """Load and validate a YAML config file."""

    if config_path is None:
        return PanScapeConfig()

    if not config_path.exists():
        raise PanScapeUsageError(f"Config file does not exist: {config_path}")

    if not config_path.is_file():
        raise PanScapeUsageError(f"Config path is not a file: {config_path}")

    payload_raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if payload_raw is None:
        payload_raw = {}

    if not isinstance(payload_raw, dict):
        raise PanScapeUsageError("Config YAML must be a key/value mapping at the top level.")

    try:
        return PanScapeConfig.model_validate(payload_raw)
    except ValidationError as exc:
        raise PanScapeUsageError(f"Invalid config file: {config_path}\n{exc}") from exc


T = TypeVar("T", bound=CommonConfig)


def merge_command_config(
    *,
    config_path: Path | None,
    section: str,
    model_cls: type[T],
    cli_overrides: Mapping[str, Any],
) -> T:
    """Merge YAML config values with explicit CLI overrides and validate."""

    root = load_config(config_path)
    section_model = getattr(root, section)

    merged: dict[str, Any] = {}
    if section_model is not None:
        merged.update(section_model.model_dump(exclude_none=True))

    for key, value in cli_overrides.items():
        if value is not None:
            merged[key] = value

    try:
        return model_cls.model_validate(merged)
    except ValidationError as exc:
        raise PanScapeUsageError(f"Invalid merged config for `{section}`:\n{exc}") from exc
