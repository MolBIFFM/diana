"""Mappings of configuration file entries to conversions."""

import math
import statistics
from typing import Callable, Collection

LOGARITHM: Callable[[Collection[float]], float] = {
    None:
        lambda changes, combination: math.log2(combination(changes)),
    2:
        lambda changes, combination: math.log2(
            combination([math.pow(2.0, change) for change in changes])),
    10:
        lambda changes, combination: math.log2(
            combination([math.pow(10.0, change) for change in changes])),
}

CHANGE_CONVERSION: Callable[[Collection[float]], float] = {
    None:
        lambda change, changes: change,
    "quantile":
        lambda quantile, changes: sorted(changes)[math.floor(quantile * (len(
            changes) - 1))],
    "standard score":
        lambda standard_score, changes: standard_score * statistics.stdev(
            changes) + statistics.mean(changes)
}
