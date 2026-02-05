from __future__ import annotations


class PanScapeError(Exception):
    """Base class for PanScape exceptions."""

    exit_code: int = 1


class PanScapeUsageError(PanScapeError):
    """Raised when command arguments or inputs are invalid."""

    exit_code = 2
