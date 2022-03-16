"""Mappings of configuration file entries to statistical tests."""
from typing import Callable, Collection

from analysis import test

ENRICHMENT_TEST: dict[str, Callable[[int, int, int, int], float]] = {
    "binomial": test.binomial,
    "hypergeometric": test.hypergeometric,
}

LOCATION_TEST: dict[str, Callable[[Collection[float], Collection[float]],
                                  float]] = {
                                      "Wilcoxon": test.wilcoxon
                                  }
