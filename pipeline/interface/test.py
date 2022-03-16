"""Mappings of configuration file entries to statistical tests."""
from analysis import test

TEST = {
    "binomial": test.binomial,
    "hypergeometric": test.hypergeometric,
    "Wilcoxon": test.wilcoxon
}
