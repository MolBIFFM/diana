"""Mappings of configuration file entries to averages."""
import statistics
from typing import Callable, Iterable, Optional

SITE_AVERAGE: dict[str, Callable[[Iterable[float]], float]] = {
    "max": max,
    "maxabs": lambda sites: max(sites, key=abs),
    "mean": statistics.mean,
    "median": statistics.median,
    "mid-range": lambda sites: 0.5 * (max(sites) + min(sites)),
    "min": min,
    "minabs": lambda sites: min(sites, key=abs),
    "sum": sum,
    "sumabs": lambda sites: sum(abs(site) for site in sites),
}

REPLICATE_AVERAGE: dict[str, Callable[[Iterable[float]], float]] = {
    "max":
        max,
    "maxabs":
        lambda replicates: max(replicates, key=abs),
    "mean":
        statistics.mean,
    "median":
        statistics.median,
    "mid-range":
        lambda replicates: 0.5 * (max(replicates) + min(replicates)),
    "min":
        min,
    "minabs":
        lambda replicates: min(replicates, key=abs),
    "sum":
        sum,
    "sumabs":
        lambda replicates: sum(abs(replicate) for replicate in replicates),
}

MODULE_SIZE_AVERAGE: dict[str, Callable[[Iterable[int]], float]] = {
    "max": max,
    "mean": statistics.mean,
    "median": statistics.median,
    "mid-range": lambda sizes: 0.5 * (max(sizes) + min(sizes)),
    "min": min,
}

CONFIDENCE_SCORE_AVERAGE: dict[Optional[str], Callable[[dict[
    str, float]], float]] = {
        None:
            lambda scores: float(bool(scores.values())),
        "max":
            lambda scores: max(scores.values()),
        "mean":
            lambda scores: statistics.mean(scores.values()),
        "median":
            lambda scores: statistics.median(scores.values()),
        "mid-range":
            lambda scores: 0.5 * (max(scores.values()) + min(scores.values())),
        "min":
            lambda scores: min(scores.values()),
        "number":
            len,
        "sum":
            lambda scores: sum(scores.values()),
    }
