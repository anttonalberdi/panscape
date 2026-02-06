# PanScape

PanScape is a Python CLI scaffold for species-level pangenome-style reference building,
read mapping, species artifact generation, and strain deconvolution workflows.

> Status: initial implementation scaffold (placeholder bioinformatics logic).

## Install

### Conda (recommended)

```bash
conda env create -f environment.yml
conda activate panscape
```

To update an existing env:

```bash
conda env update -f environment.yml --prune
```

### pip + venv

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .[dev]
```

## Quickstart

```bash
# Build species-level references (plan-only)
panscape build --genomes-tsv genomes.tsv --outdir runs/demo --dry-run

# Map reads (plan-only)
panscape map --samples-tsv samples.tsv --outdir runs/demo --dry-run

# Derive species-level artifacts
panscape species --outdir runs/demo --dry-run

# Run strain deconvolution entrypoint
panscape strain --outdir runs/demo --dry-run
```

All commands support:

- `--config PATH` (YAML)
- `--outdir PATH`
- `--threads INT`
- `--dry-run`
- `--force`
- `--log-file PATH`
- `--verbose` / `--quiet`

## Input manifests

### `genomes.tsv`

Required columns:

- `genome_id`
- `fasta_path`

Optional columns:

- `completeness`
- `contamination`
- `sample_id`

### `samples.tsv`

Required columns:

- `sample_id`
- `r1`

Optional columns:

- `r2` (for paired-end)
- additional metadata columns

## Example config (`panscape.yaml`)

```yaml
build:
  genomes_tsv: genomes.tsv
  outdir: runs/demo
  threads: 8

map:
  samples_tsv: samples.tsv
  mapper: minimap2
  outdir: runs/demo
```

CLI options override config values.

## Output layout

```text
outdir/
  panscape_manifest.json
  build/
    checkm2/quality_report.tsv
    species_clusters.tsv
    backbones/
    pangenomes/<species_id>/{genes.fna,genes.faa,gene_families.tsv,index/}
  map/<sample_id>/<species_id>.bam
  map/<sample_id>/qc.json
  species/<species_id>/{coverage_matrix.tsv,qc_summary.json}
  strain/<species_id>/K*/{abundance.tsv,gene_content.tsv,model.json}
```

## Development

```bash
pytest
ruff check .
black --check .
mypy src
```

## Notes

- External bioinformatics tools are abstracted under `src/panscape/runners/`.
- Current command implementations validate inputs, create directories, write manifests,
  print execution plans, and create placeholder outputs.
- TODO comments mark where production bioinformatics logic should be implemented.
