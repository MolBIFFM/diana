"""Mappings of configuration file entries to statistical tests."""
from enrichment import test

TEST = {
    "binomial": test.binomial,
    "hypergeometric": test.hypergeometric,
    "Wilcoxon": test.wilcoxon
}
