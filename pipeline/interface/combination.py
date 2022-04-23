"""Mappings of configuration file entries to combining functions."""
import math
import statistics
from typing import Callable, Collection

SITE_COMBINATION: dict[str, Callable[[Collection[float]], float]] = {
    "mean": statistics.mean,
    "median": statistics.median,
    "max": max,
    "maxabs": lambda sites: max(sites, key=abs),
    "min": min,
    "minabs": lambda sites: min(sites, key=abs),
    "sum": math.fsum,
    "sumabs": lambda sites: math.fsum(math.fabs(site) for site in sites),
    "increase": lambda sites: sum(site > 1.0 for site in sites) / len(sites),
    "decrease": lambda sites: sum(site < 1.0 for site in sites) / len(sites),
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
        None: lambda scores: float(bool(scores.values())),
        "mean": lambda scores: statistics.mean(scores.values()),
        "median": lambda scores: statistics.median(scores.values()),
        "max": lambda scores: max(scores.values()),
        "min": lambda scores: min(scores.values()),
        "number": len,
        "sum": lambda scores: math.fsum(scores.values()),
        **{
            database: lambda scores, database=database: scores.get(
                database, 0.0) for database in ("BioGRID", "CORUM", "IntAct", "MINT", "Reactome", "STRING")
        }
    }
