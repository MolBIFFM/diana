"""Multiple testing correction procedures."""
from typing import Hashable, Mapping


def benjamini_hochberg(p: Mapping[Hashable, float]) -> dict[Hashable, float]:
    """
    Benjamini-Hochberg method for multiple testing correction.

    Goeman, J. J. and Solari, A (2014) Multiple hypothesis testing in genomics.
        Statistics in Medicine, 33, 1946 – 1978.

    Args:
        p: Keyed p-values.

    Returns:
        Keyed Benjamini-Hochberg-adjusted p-values.
    """
    m = len(p)

    p_sorted = dict(sorted(p.items(), key=lambda item: item[1], reverse=True))

    pj = 0.0
    for i, (key, pi) in enumerate(p_sorted.items()):
        if i and pi / (m - i) > pj / (m - i + 1):
            p_sorted[key] = pj * (m - i) / (m - i + 1)
        pj = p_sorted[key]

    return {
        key: min(m * pi / (m - i), 1.0)
        for i, (key, pi) in enumerate(p_sorted.items())
    }


def benjamini_yekutieli(p: Mapping[Hashable, float]) -> dict[Hashable, float]:
    """
    Benjamini-Yekutieli method for multiple testing correction.

    Goeman, J. J. and Solari, A (2014) Multiple hypothesis testing in genomics.
        Statistics in Medicine, 33, 1946 – 1978.

    Args:
        p: Keyed p-values.

    Returns:
        Keyed Benjamini-Yekutieli-adjusted p-values.
    """
    m = len(p)
    l = sum(1 / k for k in range(1, m + 1))

    p_sorted = dict(sorted(p.items(), key=lambda item: item[1], reverse=True))

    pj = 0.0
    for i, (key, pi) in enumerate(p_sorted.items()):
        if i and pi / (m - i) > pj / (m - i + 1):
            p_sorted[key] = pj * (m - i) / (m - i + 1)
        pj = p_sorted[key]

    return {
        key: min(m * l * pi / (m - i), 1.0)
        for i, (key, pi) in enumerate(p_sorted.items())
    }


def holm(p: Mapping[Hashable, float]) -> dict[Hashable, float]:
    """
    Holm method for multiple testing correction.

    Goeman, J. J. and Solari, A (2014) Multiple hypothesis testing in genomics.
        Statistics in Medicine, 33, 1946 – 1978.

    Args:
        p: Keyed p-values.

    Returns:
        Keyed Holm-adjusted p-values.
    """
    m = len(p)

    p_sorted = dict(sorted(p.items(), key=lambda item: item[1]))

    pj = 0.0
    for i, (key, pi) in enumerate(p_sorted.items()):
        if i and (m - i) * pi < (m - i - 1) * pj:
            p_sorted[key] = (m - i - 1) * pj / (m - i)
        pj = p_sorted[key]

    return {
        key: min((m - i) * pi, 1.0)
        for i, (key, pi) in enumerate(p_sorted.items())
    }


def hommel(p: Mapping[Hashable, float]) -> dict[Hashable, float]:
    """
    Multiple testing correction according to Hommel's inequality.

    Goeman, J. J. and Solari, A (2014) Multiple hypothesis testing in genomics.
        Statistics in Medicine, 33, 1946 – 1978.

    Args:
        p: Keyed p-values.

    Returns:
        Keyed Hommel-adjusted p-values.
    """
    m = len(p)
    l = sum(1 / k for k in range(1, m + 1))

    return {
        key: min(m * l * pi / i, 1.0) for i, (key, pi) in enumerate(
            sorted(p.items(), key=lambda item: item[1]), start=1)
    }
