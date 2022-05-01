"""Mappings of configuration file entries to combining functions."""
import math
import statistics
from typing import Callable, Collection

SITE_COMBINATION: dict[str, Callable[[Collection[float]], float]] = {
    "decrease": lambda sites: sum(site < 1.0 for site in sites) / len(sites),
    "increase": lambda sites: sum(site > 1.0 for site in sites) / len(sites),
    "max": max,
    "maxabs": lambda sites: max(sites, key=abs),
    "mean": statistics.mean,
    "median": statistics.median,
    "min": min,
    "minabs": lambda sites: min(sites, key=abs),
    "sum": math.fsum,
    "sumabs": lambda sites: math.fsum(math.fabs(site) for site in sites),
}

REPLICATE_COMBINATION: dict[str, Callable[[Collection[float]], float]] = {
    "max": max,
    "maxabs": lambda sites: max(sites, key=abs),
    "mean": statistics.mean,
    "median": statistics.median,
    "min": min,
    "minabs": lambda sites: min(sites, key=abs),
    "sum": math.fsum,
    "sumabs": lambda sites: math.fsum(math.fabs(site) for site in sites),
}

MODULE_SIZE_COMBINATION: dict[str, Callable[[Collection[int]], float]] = {
    "max": max,
    "mean": statistics.mean,
    "median": statistics.median,
    "min": min,
}

CONFIDENCE_SCORE_COMBINATION: dict[str, Callable[[dict[str, float]], float]] = {
    None:
    lambda scores: float(bool(scores.values())),
    "max":
    lambda scores: max(scores.values()),
    "mean":
    lambda scores: statistics.mean(scores.values()),
    "median":
    lambda scores: statistics.median(scores.values()),
    "min":
    lambda scores: min(scores.values()),
    "number":
    len,
    "sum": lambda scores: math.fsum(scores.values(
    )),
    **{
        database: lambda scores, database=database: scores.get(database, 0.0)
        for database in ("BioGRID", "CORUM", "IntAct", "MINT", "Reactome", "STRING")
    }
}
