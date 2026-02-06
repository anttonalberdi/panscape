from __future__ import annotations

import shutil
from pathlib import Path
from typing import Mapping

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
        input_dir: Path | None = None,
        input_paths: list[Path] | None = None,
        output_dir: Path,
        threads: int,
        database_path: Path | None,
        genes: bool = False,
        extension: str | None = None,
        tmp_dir: Path | None = None,
        env: Mapping[str, str] | None = None,
        force: bool = False,
        dry_run: bool = False,
    ) -> CommandResult:
        if (input_dir is None) == (input_paths is None):
            raise ValueError("Provide exactly one of `input_dir` or `input_paths`.")
        if input_paths is not None and len(input_paths) == 0:
            raise ValueError("`input_paths` cannot be empty.")

        args: list[str | Path] = [
            "predict",
            "--input",
        ]
        if input_paths is not None:
            args.extend(input_paths)
        else:
            assert input_dir is not None
            args.append(input_dir)

        args.extend(
            [
                "--output-directory",
                output_dir,
                "--threads",
                str(threads),
            ]
        )
        if genes:
            args.append("--genes")
        if database_path is not None:
            args.extend(["--database_path", database_path])
        if extension is not None and input_paths is None:
            args.extend(["--extension", extension])
        if force:
            args.append("--force")

        env_payload = dict(env) if env is not None else {}
        if tmp_dir is not None:
            tmp_value = str(tmp_dir)
            env_payload.update({"TMPDIR": tmp_value, "TMP": tmp_value, "TEMP": tmp_value})

        return self.run(args, dry_run=dry_run, env=env_payload or None)
