"""Mappings of configuration file entries to statistical tests."""
from typing import Callable, Collection

import scipy.stats

# Hypothesis tests for enrichment.
ENRICHMENT_TEST: dict[tuple[str, bool], Callable[
    [int, int, int, int], float]] = {
        ("binomial", False):
            lambda k, M, n, N: float(scipy.stats.binom.cdf(k, N, n / M)),
        ("binomial", True):
            lambda k, M, n, N: float(scipy.stats.binom.sf(k - 1, N, n / M)),
        ("hypergeometric", False):
            lambda k, M, n, N: float(scipy.stats.hypergeom.cdf(k, M, n, N)),
        ("hypergeometric", True):
            lambda k, M, n, N: float(scipy.stats.hypergeom.sf(k - 1, M, n, N)),
    }

# Hypothesis tests for equality of location.
LOCATION_TEST: dict[tuple[str, bool, bool], Callable[
    [Collection[float], Collection[float]], float]] = {
        ("Mann-Whitney-Wilcoxon", False, False):
            lambda x, y: float(
                scipy.stats.mannwhitneyu(x,
                                         y,
                                         use_continuity=False,
                                         alternative="less",
                                         method="asymptotic").pvalue),
        ("Mann-Whitney-Wilcoxon", False, True):
            lambda x, y: float(
                scipy.stats.mannwhitneyu([abs(xi)
                                          for xi in x], [abs(yi) for yi in y],
                                         use_continuity=False,
                                         alternative="less",
                                         method="asymptotic").pvalue),
        ("Mann-Whitney-Wilcoxon", True, False):
            lambda x, y: float(
                scipy.stats.mannwhitneyu(x,
                                         y,
                                         use_continuity=False,
                                         alternative="greater",
                                         method="asymptotic").pvalue),
        ("Mann-Whitney-Wilcoxon", True, True):
            lambda x, y: float(
                scipy.stats.mannwhitneyu([abs(xi)
                                          for xi in x], [abs(yi) for yi in y],
                                         use_continuity=False,
                                         alternative="greater",
                                         method="asymptotic").pvalue),
        ("Welch", False, False):
            lambda x, y: float(
                scipy.stats.ttest_ind(x, y, equal_var=False, alternative="less")
                .pvalue),
        ("Welch", False, True):
            lambda x, y: float(
                scipy.stats.ttest_ind([abs(xi)
                                       for xi in x], [abs(yi) for yi in y],
                                      equal_var=False,
                                      alternative="less").pvalue),
        ("Welch", True, False):
            lambda x, y: float(
                scipy.stats.ttest_ind(
                    x, y, equal_var=False, alternative="greater").pvalue),
        ("Welch", True, True):
            lambda x, y: float(
                scipy.stats.ttest_ind([abs(xi)
                                       for xi in x], [abs(yi) for yi in y],
                                      equal_var=False,
                                      alternative="greater").pvalue),
    }
