"""Mappings of configuration file entries to community detection algorithms."""
from modularization import modularization

ALGORITHM = {
    "Clauset-Newman-Moore": modularization.clauset_newman_moore,
    "Louvain": modularization.louvain,
}
