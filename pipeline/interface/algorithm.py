from ..modularization import clauset_newman_moore
from ..modularization import louvain

ALGORITHM = {
    "Clauset-Newman-Moore": clauset_newman_moore.clauset_newman_moore,
    "Louvain": louvain.louvain,
}
