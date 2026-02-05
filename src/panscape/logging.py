from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path

from rich.logging import RichHandler


class JsonLogFormatter(logging.Formatter):
    """Emit one JSON object per log line for machine-friendly logs."""

    def format(self, record: logging.LogRecord) -> str:
        payload: dict[str, object] = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
        }
        if record.exc_info is not None:
            payload["exception"] = self.formatException(record.exc_info)
        return json.dumps(payload, ensure_ascii=True)


def _resolve_level(verbose: bool, quiet: bool) -> int:
    if quiet:
        return logging.ERROR
    if verbose:
        return logging.DEBUG
    return logging.INFO


def configure_logging(
    *,
    verbose: bool = False,
    quiet: bool = False,
    log_file: Path | None = None,
) -> None:
    """Configure console + optional JSON file logging."""

    root = logging.getLogger()
    root.handlers.clear()
    root.setLevel(_resolve_level(verbose, quiet))

    rich_handler = RichHandler(rich_tracebacks=True, show_path=False, markup=False)
    rich_handler.setFormatter(logging.Formatter("%(message)s"))
    root.addHandler(rich_handler)

    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setFormatter(JsonLogFormatter())
        root.addHandler(file_handler)


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(name)
