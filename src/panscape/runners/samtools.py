from __future__ import annotations

from pathlib import Path

from panscape.runners.base import ToolRunner
from panscape.utils.subprocess import CommandResult


class SamtoolsRunner(ToolRunner):
    def __init__(self, executable: str = "samtools") -> None:
        super().__init__(executable)

    def sort(self, *, input_bam: Path, output_bam: Path, threads: int, dry_run: bool = False) -> CommandResult:
        return self.run(
            ["sort", "-@", str(threads), "-o", str(output_bam), str(input_bam)],
            dry_run=dry_run,
        )

    def index(self, *, input_bam: Path, dry_run: bool = False) -> CommandResult:
        return self.run(["index", str(input_bam)], dry_run=dry_run)
