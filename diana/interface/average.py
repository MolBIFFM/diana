"""Mappings of configuration file entries to averages."""
import math
import statistics
from typing import Callable, Iterable, Optional

SITE_AVERAGE: dict[str, Callable[[Iterable[float]], float]] = {
    "max":
        max,
    "maxabslog":
        lambda sites: max(sites, key=lambda site: abs(math.log2(site))),
    "mean":
        statistics.mean,
    "median":
        statistics.median,
    "min":
        min,
    "minabslog":
        lambda sites: min(sites, key=lambda site: abs(math.log2(site))),
    "sum":
        sum,
    "sumabslog":
        lambda sites: sum(abs(math.log2(site)) for site in sites),
}

REPLICATE_AVERAGE: dict[str, Callable[[Iterable[float]], float]] = {
    "max":
        max,
    "maxabslog":
        lambda replicates: max(replicates,
                               key=lambda replicate: abs(math.log2(replicate))),
    "mean":
        statistics.mean,
    "median":
        statistics.median,
    "min":
        min,
    "minabslog":
        lambda replicates: min(replicates,
                               key=lambda replicate: abs(math.log2(replicate))),
    "sum":
        sum,
    "sumabslog":
        lambda replicates: sum(
            abs(math.log2(replicate)) for replicate in replicates),
}

MODULE_SIZE_AVERAGE: dict[str, Callable[[Iterable[int]], float]] = {
    "max": max,
    "mean": statistics.mean,
    "median": statistics.median,
    "min": min,
}

CONFIDENCE_SCORE_AVERAGE: dict[Optional[str], Callable[[dict[
    str, float]], float]] = {
        None: lambda scores: float(bool(scores.values())),
        "max": lambda scores: max(scores.values()),
        "mean": lambda scores: statistics.mean(scores.values()),
        "median": lambda scores: statistics.median(scores.values()),
        "min": lambda scores: min(scores.values()),
        "number": len,
        "sum": lambda scores: sum(scores.values()),
    }
