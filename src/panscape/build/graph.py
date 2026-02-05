from __future__ import annotations

from collections import defaultdict, deque
from typing import Iterable


def connected_components(nodes: Iterable[str], edges: Iterable[tuple[str, str]]) -> list[list[str]]:
    """Return deterministic connected components from an undirected graph."""

    node_set = set(nodes)
    adjacency: dict[str, set[str]] = defaultdict(set)

    for node in node_set:
        adjacency[node]

    for left, right in edges:
        if left not in node_set or right not in node_set:
            continue
        adjacency[left].add(right)
        adjacency[right].add(left)

    visited: set[str] = set()
    components: list[list[str]] = []

    for seed in sorted(node_set):
        if seed in visited:
            continue

        queue: deque[str] = deque([seed])
        visited.add(seed)
        component: list[str] = []

        while queue:
            current = queue.popleft()
            component.append(current)
            for neighbor in sorted(adjacency[current]):
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)

        components.append(sorted(component))

    components.sort(key=lambda comp: (comp[0], len(comp), tuple(comp)))
    return components
