"""Mappings of configuration file entries to combining functions."""
import math
from typing import Callable, Collection
import statistics

SITE_COMBINATION: dict[str, Callable[[Collection[float]], float]] = {
    "mean":
        lambda changes: math.log2(
            statistics.mean(math.pow(2.0, change) for change in changes)),
    "median":
        lambda changes: math.log2(
            statistics.median(math.pow(2.0, change) for change in changes)),
    "max":
        max,
    "absmax":
        lambda changes: max(changes, key=abs),
    "min":
        min,
    "absmin":
        lambda changes: min(changes, key=abs),
    "sum":
        lambda changes: math.log2(
            sum(math.pow(2.0, change) for change in changes)),
    "abssum":
        lambda changes: math.log2(
            sum(abs(math.pow(2.0, change)) for change in changes)),
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

CONFIDENCE_SCORE_COMBINATION: dict[str, Callable[[Collection[float]], float]] = {
    "mean":
        lambda confidence_scores: statistics.mean(confidence_scores.values()),
    "median":
        lambda confidence_scores: statistics.median(confidence_scores.values()),
    "max":
        lambda confidence_scores: max(confidence_scores.values()),
    "min":
        lambda confidence_scores: min(confidence_scores.values()),
    "number":
        len,
    "sum":
        lambda confidence_scores: sum(confidence_scores.values()),
    **{
        database: lambda confidence_scores, database=database: confidence_scores.get(
            database, 0.0) for database in {
            "BioGRID", "IntAct", "MINT", "Reactome", "STRING"
        }
    }, None:
        lambda confidence_scores: float(bool(confidence_scores.values()))
}
