from ..utilities import modularization

ALGORITHM = {
    "Clauset-Newman-Moore": modularization.clauset_newman_moore,
    "Louvain": modularization.louvain,
}
