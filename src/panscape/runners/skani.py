from __future__ import annotations

import shutil
from pathlib import Path

from panscape.runners.base import ToolRunner
from panscape.utils.subprocess import CommandResult


class SkaniRunner(ToolRunner):
    """Wrapper around skani commands for ANI calculations."""

    def __init__(self, executable: str = "skani") -> None:
        super().__init__(executable)

    def is_available(self) -> bool:
        return shutil.which(self.executable) is not None

    def version(self, *, dry_run: bool = False) -> str:
        result = self.run(["--version"], dry_run=dry_run, check=False)
        version_line = result.stdout.strip() or result.stderr.strip()
        return version_line or "unknown"

    def dist_pair(self, *, query: Path, reference: Path, dry_run: bool = False) -> CommandResult:
        return self.run(
            [
                "dist",
                "--query",
                str(query),
                "--ref",
                str(reference),
                "--min-af",
                "0",
            ],
            dry_run=dry_run,
        )
