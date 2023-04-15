"""Mappings of configuration file entries to community detection algorithms."""
from typing import Callable, Hashable

import networkx as nx
from algorithms import modularization

# Community detection algorithms.
ALGORITHM: dict[str, Callable[[nx.Graph, float, str], list[set[Hashable]]]] = {
    "Clauset-Newman-Moore": modularization.clauset_newman_moore,
    "Louvain": modularization.louvain,
}
