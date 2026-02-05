from __future__ import annotations

import shutil
from pathlib import Path

from panscape.runners.base import ToolRunner
from panscape.utils.subprocess import CommandResult


class MashRunner(ToolRunner):
    """Wrapper around Mash commands used by `panscape build`."""

    def __init__(self, executable: str = "mash") -> None:
        super().__init__(executable)

    def is_available(self) -> bool:
        return shutil.which(self.executable) is not None

    def version(self, *, dry_run: bool = False) -> str:
        result = self.run(["--version"], dry_run=dry_run, check=False)
        version_line = result.stdout.strip() or result.stderr.strip()
        return version_line or "unknown"

    def sketch(
        self,
        *,
        input_fastas: list[Path],
        out_prefix: Path,
        kmer_size: int,
        sketch_size: int,
        dry_run: bool = False,
    ) -> CommandResult:
        return self.run(
            [
                "sketch",
                "-k",
                str(kmer_size),
                "-s",
                str(sketch_size),
                "-o",
                str(out_prefix),
                *[str(path) for path in input_fastas],
            ],
            dry_run=dry_run,
        )

    def dist(self, *, sketch_path: Path, dry_run: bool = False) -> CommandResult:
        return self.run(["dist", str(sketch_path), str(sketch_path)], dry_run=dry_run)
