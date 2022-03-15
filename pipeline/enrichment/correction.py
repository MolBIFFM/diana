"""Multiple testing correction procedures."""
from typing import Hashable


def benjamini_hochberg(
        p_values: dict[Hashable, float]) -> dict[Hashable, float]:
    """
    Benjamini-Hochberg procedure for multiple testing correction.

    Y. Benjamini, R. Heller and D. Yekutieli, "Selective inference in complex 
    research", Philosophical Transactions of the Royal Society A: Mathematical, 
    Physical and Engineering Sciences, 2009

    Args:
        p_values: Keyed p-values.

    Returns:
        Benjamini-Hochberg-adjusted p-values.
    """
    m = len(p_values)
    sorted_p_values = sorted((list(item) for item in p_values.items()),
                             key=lambda item: item[1])

    for i in range(m - 1, 0, -1):
        if sorted_p_values[i - 1][1] / i > sorted_p_values[i][1] / (i + 1):
            sorted_p_values[i - 1][1] = sorted_p_values[i][1] * i / (i + 1)

    return {
        key: min(m * p_value / (i + 1), 1.0)
        for i, (key, p_value) in enumerate(sorted_p_values)
    }


def bonferroni(p_values: dict[Hashable, float]) -> dict[Hashable, float]:
    """
    Bonferroni procedure for multiple testing correction.

    Args:
        p_values: Keyed p-values.

    Returns:
        Bonferroni-adjusted p-values.
    """
    m = len(p_values)

    return {key: min(m * p_value, 1.0) for key, p_value in p_values.items()}
