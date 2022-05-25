"""Mappings of configuration file entries to conversions."""

import math
import statistics
from typing import Callable, Collection, Iterable, Optional

LOGARITHM: dict[Optional[int], Callable[
    [Iterable[float], Callable[[Iterable[float]], float]], float]] = {
        None:
            lambda measurements, combination: math.log2(
                combination(measurements)),
        2:
            lambda measurements, combination: math.log2(
                combination([
                    math.pow(2.0, measurement) for measurement in measurements
                ])),
        10:
            lambda measurements, combination: math.log2(
                combination([
                    math.pow(10.0, measurement) for measurement in measurements
                ])),
    }

MEASUREMENT_CONVERSION: dict[Optional[str], Callable[
    [float, Collection[float]], float]] = {
        None:
            lambda measurement, _: measurement,
        "percentile":
            lambda percentile, measurements: sorted(measurements)[math.floor(
                percentile / 100.0 * (len(measurements) - 1))],
        "quantile":
            lambda quantile, measurements: sorted(measurements)[math.floor(
                quantile * (len(measurements) - 1))],
        "ratio":
            lambda ratio, _: math.log2(ratio),
        "standard score":
            lambda standard_score, measurements: standard_score * statistics.
            stdev(measurements) + statistics.mean(measurements),
    }
