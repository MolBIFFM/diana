"""Mappings of configuration file entries to averages."""
import math
import statistics
from typing import Callable, Iterable, Optional

SITE_AVERAGE: dict[str, Callable[[Iterable[float]], float]] = {
    "maximum":
        max,
    "maximum absolute logarithm":
        lambda sites: max(sites, key=lambda site: abs(math.log2(site))),
    "mean":
        statistics.mean,
    "median":
        statistics.median,
    "minimum":
        min,
    "minimum absolute logarithm":
        lambda sites: min(sites, key=lambda site: abs(math.log2(site))),
    "sum":
        sum,
    "sum absolute logarithm":
        lambda sites: sum(abs(math.log2(site)) for site in sites),
}

REPLICATE_AVERAGE: dict[str, Callable[[Iterable[float]], float]] = {
    "maximum":
        max,
    "maximum absolute logarithm":
        lambda replicates: max(replicates,
                               key=lambda replicate: abs(math.log2(replicate))),
    "mean":
        statistics.mean,
    "median":
        statistics.median,
    "minimum":
        min,
    "minimum absolute logarithm":
        lambda replicates: min(replicates,
                               key=lambda replicate: abs(math.log2(replicate))),
    "sum":
        sum,
    "sum absolute logarithm":
        lambda replicates: sum(
            abs(math.log2(replicate)) for replicate in replicates),
}

MODULE_SIZE_AVERAGE: dict[str, Callable[[Iterable[int]], float]] = {
    "maximum": max,
    "mean": statistics.mean,
    "median": statistics.median,
    "minimum": min,
}

CONFIDENCE_SCORE_AVERAGE: dict[Optional[str], Callable[[dict[
    str, float]], float]] = {
        None: lambda scores: float(bool(scores.values())),
        "maximum": lambda scores: max(scores.values()),
        "mean": lambda scores: statistics.mean(scores.values()),
        "median": lambda scores: statistics.median(scores.values()),
        "minimum": lambda scores: min(scores.values()),
        "number": len,
        "sum": lambda scores: sum(scores.values()),
    }
