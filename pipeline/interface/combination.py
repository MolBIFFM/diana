"""Mappings of configuration file entries to combining functions."""
import math
import statistics
from typing import Callable, Collection

SITE_COMBINATION: dict[str, Callable[[Collection[float]], float]] = {
    "mean": statistics.mean,
    "median": statistics.median,
    "max": max,
    "absmax": lambda sites: max(sites, key=abs),
    "min": min,
    "absmin": lambda sites: min(sites, key=abs),
    "sum": math.fsum,
    "abssum": lambda sites: math.fsum(math.fabs(site) for site in sites),
}

REPLICATE_COMBINATION: dict[str, Callable[[Collection[float]], float]] = {
    "mean": statistics.mean,
    "median": statistics.median,
}

MODULE_SIZE_COMBINATION: dict[str, Callable[[Collection[int]], float]] = {
    "mean": statistics.mean,
    "median": statistics.median,
    "max": max,
    "min": min,
}

CONFIDENCE_SCORE_COMBINATION: dict[str, Callable[
    [Collection[float]], float]] = {
        "mean": lambda scores: statistics.mean(scores.values()),
        "median": lambda scores: statistics.median(scores.values()),
        "max": lambda scores: max(scores.values()),
        "min": lambda scores: min(scores.values()),
        "number": len,
        "sum": lambda scores: math.fsum(scores.values()),
        **{
            database: lambda scores, database=database: scores.get(
                database, 0.0) for database in {
                "BioGRID", "IntAct", "MINT", "Reactome", "STRING"
            }
        }, None: lambda scores: float(bool(scores.values()))
    }
