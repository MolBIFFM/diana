"""Mappings of configuration file entries to statistical tests."""
from typing import Callable, Collection

import scipy.stats

ENRICHMENT_TEST: dict[tuple[str, bool], Callable[
    [int, int, int, int], float]] = {
        ("binomial", False):
            lambda k, M, n, N: scipy.stats.binom.cdf(k, N, n / M),
        ("binomial", True):
            lambda k, M, n, N: scipy.stats.binom.sf(k - 1, N, n / M),
        ("hypergeometric", False):
            scipy.stats.hypergeom.cdf,
        ("hypergeometric", True):
            lambda k, M, n, N: scipy.stats.hypergeom.sf(k - 1, M, n, N),
    }

LOCATION_TEST: dict[tuple[str, bool, bool], Callable[
    [Collection[float], Collection[float]], float]] = {
        ("Welch", False, False):
            lambda x, y: scipy.stats.ttest_ind(
                x, y, equal_var=False, alternative="less").pvalue,
        ("Welch", False, True):
            lambda x, y: scipy.stats.ttest_ind([abs(xi) for xi in x],
                                               [abs(yi) for yi in y],
                                               equal_var=False,
                                               alternative="less").pvalue,
        ("Welch", True, False):
            lambda x, y: scipy.stats.ttest_ind(
                x, y, equal_var=False, alternative="greater").pvalue,
        ("Welch", True, True):
            lambda x, y: scipy.stats.ttest_ind([abs(xi) for xi in x],
                                               [abs(yi) for yi in y],
                                               equal_var=False,
                                               alternative="greater").pvalue,
        ("Wilcoxon", False, False):
            lambda x, y: scipy.stats.ranksums(x, y, alternative="less").pvalue,
        ("Wilcoxon", False, True):
            lambda x, y: scipy.stats.ranksums([abs(xi) for xi in x],
                                              [abs(yi) for yi in y],
                                              alternative="less").pvalue,
        ("Wilcoxon", True, False):
            lambda x, y: scipy.stats.ranksums(x, y, alternative="greater").
            pvalue,
        ("Wilcoxon", True, True):
            lambda x, y: scipy.stats.ranksums([abs(xi) for xi in x],
                                              [abs(yi) for yi in y],
                                              alternative="greater").pvalue,
    }
