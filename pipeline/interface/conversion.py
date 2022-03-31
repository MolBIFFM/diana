"""Mappings of configuration file entries to conversions."""

import math
import statistics
from typing import Callable, Collection

LOGARITHM: Callable[[Collection[float]], float] = {
    None:
        lambda measurements, combination: math.log2(combination(measurements)),
    2:
        lambda measurements, combination: math.log2(
            combination(
                [math.pow(2.0, measurement) for measurement in measurements])),
    10:
        lambda measurements, combination: math.log2(
            combination(
                [math.pow(10.0, measurement) for measurement in measurements])),
}

MEASUREMENT_CONVERSION: Callable[[Collection[float]], float] = {
    None:
        lambda measurement, measurements: measurement,
    "quantile":
        lambda quantile, measurements: sorted(measurements)[math.floor(
            quantile * (len(measurements) - 1))],
    "standard score":
        lambda standard_score, measurements: standard_score * statistics.stdev(
            measurements) + statistics.mean(measurements)
}
