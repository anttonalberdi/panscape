from __future__ import annotations

from pathlib import Path

from panscape.runners.base import ToolRunner
from panscape.utils.subprocess import CommandResult


class MapperRunner(ToolRunner):
    """Generic mapper adapter so workflow code is mapper-agnostic."""

    def __init__(self, executable: str = "minimap2") -> None:
        super().__init__(executable)

    def plan_map_args(
        self,
        *,
        reference_fasta: Path,
        r1: Path,
        r2: Path | None,
        output_bam: Path,
        threads: int,
    ) -> list[str]:
        reads = [str(r1)]
        if r2 is not None:
            reads.append(str(r2))

        # TODO: Implement tool-specific argument builders for minimap2/bwa/etc.
        return [
            "-t",
            str(threads),
            "--reference",
            str(reference_fasta),
            "--reads",
            *reads,
            "--output-bam",
            str(output_bam),
        ]

    def map_reads(
        self,
        *,
        reference_fasta: Path,
        r1: Path,
        r2: Path | None,
        output_bam: Path,
        threads: int,
        dry_run: bool = False,
    ) -> CommandResult:
        return self.run(
            self.plan_map_args(
                reference_fasta=reference_fasta,
                r1=r1,
                r2=r2,
                output_bam=output_bam,
                threads=threads,
            ),
            dry_run=dry_run,
        )
