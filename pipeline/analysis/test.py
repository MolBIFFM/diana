"""Interfaces to statistical tests."""
from typing import Collection
import scipy.stats


def binomial(k: int, M: int, n: int, N: int) -> float:
    """
    Returns the right-tail p-value for k successes in N trials with a success 
    rate of n/M.
    """
    return scipy.stats.binom.sf(k - 1, N, n / M)


def hypergeometric(k: int, M: int, n: int, N: int) -> float:
    """
    Returns the right-tail p-value for k successes in N draws for n elements 
    from a population of size M.
    """
    return scipy.stats.hypergeom.sf(k - 1, M, n, N)


def wilcoxon(x: Collection[float], y: Collection[float]) -> float:
    """
    Returns the p-value of the Wilcoxon rank-sum test of x and y.
    """
    return scipy.stats.ranksums(x, y).pvalue
