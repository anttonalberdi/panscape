from __future__ import annotations

from pathlib import Path

from panscape.runners.base import ToolRunner
from panscape.utils.subprocess import CommandResult


class ANIRunner(ToolRunner):
    """Generic ANI adapter (e.g., skani, fastANI)."""

    def __init__(self, executable: str = "skani") -> None:
        super().__init__(executable)

    def pairwise(
        self,
        *,
        query_list: Path,
        reference_list: Path,
        output_tsv: Path,
        dry_run: bool = False,
    ) -> CommandResult:
        # TODO: Support fastANI/skani-specific arguments and output parsing.
        return self.run(
            [
                "search",
                "--ql",
                str(query_list),
                "--rl",
                str(reference_list),
                "--output",
                str(output_tsv),
            ],
            dry_run=dry_run,
        )
