from __future__ import annotations

import logging
from pathlib import Path
from typing import Mapping, Sequence

from panscape.utils.subprocess import CommandResult, run_command


class ToolRunner:
    """Base abstraction for external tools with dry-run aware execution."""

    def __init__(self, executable: str, *, logger: logging.Logger | None = None) -> None:
        self.executable = executable
        self.logger = logger

    def command(self, args: Sequence[str | Path]) -> list[str]:
        return [self.executable, *[str(arg) for arg in args]]

    def run(
        self,
        args: Sequence[str | Path],
        *,
        dry_run: bool = False,
        cwd: Path | None = None,
        env: Mapping[str, str] | None = None,
        check: bool = True,
    ) -> CommandResult:
        return run_command(
            self.command(args),
            dry_run=dry_run,
            cwd=cwd,
            env=env,
            check=check,
            logger=self.logger,
        )
