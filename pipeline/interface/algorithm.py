from pipeline.modularization import clauset_newman_moore, louvain

ALGORITHM = {
    "Clauset-Newman-Moore": clauset_newman_moore.clauset_newman_moore,
    "Louvain": louvain.louvain,
}
