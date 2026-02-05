from __future__ import annotations

import shutil
from pathlib import Path

from panscape.runners.base import ToolRunner
from panscape.utils.subprocess import CommandExecutionError, CommandResult


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
        # Newer skani releases expect positional QUERY and REFERENCE inputs.
        positional = self.run(
            ["dist", str(query), str(reference), "--min-af", "0"],
            dry_run=dry_run,
            check=False,
        )
        if dry_run or positional.returncode == 0:
            return positional

        # Legacy fallback for older skani variants.
        legacy = self.run(
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
            check=False,
        )
        if legacy.returncode == 0:
            return legacy

        raise CommandExecutionError(
            "skani dist failed with both positional and legacy flag syntaxes.\n"
            f"positional stderr: {positional.stderr.strip()}\n"
            f"legacy stderr: {legacy.stderr.strip()}"
        )
