from __future__ import annotations

import shutil
from pathlib import Path

from panscape.runners.base import ToolRunner
from panscape.utils.subprocess import CommandResult


class CheckM2Runner(ToolRunner):
    """Wrapper around CheckM2 quality prediction."""

    def __init__(self, executable: str = "checkm2") -> None:
        super().__init__(executable)

    def is_available(self) -> bool:
        return shutil.which(self.executable) is not None

    def version(self, *, dry_run: bool = False) -> str:
        result = self.run(["--version"], dry_run=dry_run, check=False)
        version_line = result.stdout.strip() or result.stderr.strip()
        return version_line or "unknown"

    def predict(
        self,
        *,
        input_dir: Path,
        output_dir: Path,
        threads: int,
        database_path: Path | None,
        extension: str | None = None,
        force: bool = False,
        dry_run: bool = False,
    ) -> CommandResult:
        args: list[str | Path] = [
            "predict",
            "--input",
            input_dir,
            "--output-directory",
            output_dir,
            "--threads",
            str(threads),
        ]
        if database_path is not None:
            args.extend(["--database_path", database_path])
        if extension is not None:
            args.extend(["--extension", extension])
        if force:
            args.append("--force")
        return self.run(args, dry_run=dry_run)
