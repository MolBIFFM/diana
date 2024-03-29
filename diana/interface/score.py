"""Mappings of configuration file entries to scores."""

import math
import statistics
from typing import Callable, Collection, Optional

# Conversions from ratios or logarithms to binary logarithm.
LOGARITHM: dict[Optional[int], Callable[[float], float]] = {
    None: math.log2,
    2: lambda measurement: measurement,
    10: lambda measurement: math.log10(measurement) / math.log10(2.0)
}

# Conversions of measurements.
MEASUREMENT_SCORE: dict[Optional[str], Callable[
    [float, Collection[float]], float]] = {
        None:
            lambda measurement, _: measurement,
        "quantile":
            lambda quantile, measurements: sorted(measurements)[math.ceil(
                quantile * (len(measurements) - 1))],
        "ratio":
            lambda ratio, _: math.log2(ratio),
        "standard score":
            lambda standard_score, measurements: standard_score * statistics.
            stdev(measurements) + statistics.mean(measurements),
    }
