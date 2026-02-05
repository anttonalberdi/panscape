from __future__ import annotations

import logging
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence


@dataclass(frozen=True, slots=True)
class CommandResult:
    command: list[str]
    returncode: int
    stdout: str
    stderr: str
    dry_run: bool


class CommandExecutionError(RuntimeError):
    pass


def shell_join(command: Sequence[str]) -> str:
    return shlex.join(list(command))


def run_command(
    command: Sequence[str],
    *,
    dry_run: bool = False,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
    check: bool = True,
    logger: logging.Logger | None = None,
) -> CommandResult:
    command_list = list(command)
    cmd_text = shell_join(command_list)

    if logger is not None:
        logger.debug("Executing command: %s", cmd_text)

    if dry_run:
        return CommandResult(command=command_list, returncode=0, stdout="", stderr="", dry_run=True)

    completed = subprocess.run(
        command_list,
        cwd=str(cwd) if cwd is not None else None,
        env=dict(env) if env is not None else None,
        capture_output=True,
        text=True,
        check=False,
    )

    result = CommandResult(
        command=command_list,
        returncode=completed.returncode,
        stdout=completed.stdout,
        stderr=completed.stderr,
        dry_run=False,
    )

    if check and completed.returncode != 0:
        raise CommandExecutionError(
            f"Command failed with exit code {completed.returncode}: {cmd_text}\n{completed.stderr.strip()}"
        )

    return result
