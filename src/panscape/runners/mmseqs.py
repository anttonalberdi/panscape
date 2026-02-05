from __future__ import annotations

from pathlib import Path

from panscape.runners.base import ToolRunner
from panscape.utils.subprocess import CommandResult


class MMseqsRunner(ToolRunner):
    def __init__(self, executable: str = "mmseqs") -> None:
        super().__init__(executable)

    def cluster(
        self,
        *,
        input_db: Path,
        output_db: Path,
        tmp_dir: Path,
        min_seq_id: float,
        dry_run: bool = False,
    ) -> CommandResult:
        # TODO: Replace with production-grade MMseqs2 workflow and parameters.
        return self.run(
            [
                "cluster",
                str(input_db),
                str(output_db),
                str(tmp_dir),
                "--min-seq-id",
                str(min_seq_id),
            ],
            dry_run=dry_run,
        )
