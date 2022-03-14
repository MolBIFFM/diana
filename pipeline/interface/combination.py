import statistics

SITE_COMBINATION = {
    "mean": statistics.mean,
    "median": statistics.median,
    "max": max,
    "absmax": lambda changes: max(changes, key=abs),
    "min": min,
    "absmin": lambda changes: min(changes, key=abs),
    "sum": sum,
    "abssum": lambda changes: sum(abs(change) for change in changes),
}

REPLICATE_COMBINATION = {
    "mean": statistics.mean,
    "median": statistics.median,
}

MODULE_SIZE_COMBINATION = {
    "mean": statistics.mean,
    "median": statistics.median,
    "max": max,
    "min": min,
}

CONFIDENCE_SCORE_COMBINATION = {
    "mean":
    lambda confidence_scores: statistics.mean(confidence_scores.values()),
    "median":
    lambda confidence_scores: statistics.median(confidence_scores.values()),
    "max":
    lambda confidence_scores: max(confidence_scores.values()),
    "min":
    lambda confidence_scores: min(confidence_scores.values()),
    "number":
    lambda confidence_scores: len(confidence_scores),
    "sum":
    lambda confidence_scores: sum(confidence_scores.values())**{
        database: lambda confidence_scores: confidence_scores.get(
            database, 0.0)
        for database in {"BioGRID", "IntAct", "MINT", "Reactome", "STRING"}
    },
    None:
    lambda confidence_scores: float(bool(confidence_scores))
}
