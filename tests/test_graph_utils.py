from __future__ import annotations

from panscape.build.graph import connected_components


def test_connected_components_deterministic() -> None:
    nodes = ["g3", "g2", "g1", "g4"]
    edges = [("g1", "g2"), ("g2", "g1"), ("g3", "g4")]

    components = connected_components(nodes, edges)

    assert components == [["g1", "g2"], ["g3", "g4"]]
