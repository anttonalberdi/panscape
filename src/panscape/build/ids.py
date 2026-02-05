from __future__ import annotations

from typing import Callable, Iterable, TypeVar

T = TypeVar("T")


def assign_prefixed_ids(
    items: Iterable[T],
    *,
    prefix: str,
    width: int = 6,
    key: Callable[[T], tuple | str] | None = None,
) -> dict[T, str]:
    """Assign deterministic IDs such as species_000001."""

    sorted_items = sorted(items, key=key)
    assignments: dict[T, str] = {}
    for idx, item in enumerate(sorted_items, start=1):
        assignments[item] = f"{prefix}_{idx:0{width}d}"
    return assignments
