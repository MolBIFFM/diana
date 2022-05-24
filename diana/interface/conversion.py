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
            lambda measurement, measurements: measurement,
        "ratio":
            lambda measurement, measurements: math.pow(2.0, measurement),
        "percentile":
            lambda measurement, measurements: sorted(measurements)[math.floor(
                measurement / 100.0 * (len(measurements) - 1))],
        "quantile":
            lambda measurement, measurements: sorted(measurements)[math.floor(
                measurement * (len(measurements) - 1))],
        "standard score":
            lambda measurement, measurements: measurement * statistics.stdev(
                measurements) + statistics.mean(measurements)
    }
