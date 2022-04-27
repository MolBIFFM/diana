"""Mappings of configuration file entries to statistical tests."""
from typing import Callable, Collection

import scipy.stats

ENRICHMENT_TEST: dict[str, Callable[[int, int, int, int], float]] = {
    "binomial":
        lambda k, M, n, N: scipy.stats.binom.sf(k - 1, N, n / M),
    "hypergeometric":
        lambda k, M, n, N: scipy.stats.hypergeom.sf(k - 1, M, n, N)
}

LOCATION_TEST: dict[str, Callable[
    [Collection[float], Collection[float]], float]] = {
        "Welch":
            lambda x, y: scipy.stats.ttest_ind(x, y, equal_var=False).pvalue,
        "Wilcoxon":
            lambda x, y: scipy.stats.ranksums(x, y).pvalue
    }
