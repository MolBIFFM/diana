import statistics

COMBINE_CHANGES = {
    "down": lambda changes: len([change for change in changes if change < 0.0])
    / len(changes),
    "max": max,
    "maxabs": lambda changes: max(changes, key=abs),
    "mean": statistics.mean,
    "median": statistics.median,
    "min": min,
    "minabs": lambda changes: min(changes, key=abs),
    "number": lambda changes: len(changes),
    "sum": sum,
    "sumabs": lambda changes: sum(abs(change) for change in changes),
    "up": lambda changes: len([change for change in changes if change > 0.0])
    / len(changes),
}

COMBINE_MODULE_SIZES = {
    "max": max,
    "mean": statistics.mean,
    "median": statistics.median,
    "min": min,
}


COMBINE_CONFIDENCE_SCORES = {
    "max": lambda confidence_scores: max(confidence_scores.values()),
    "mean": lambda confidence_scores: statistics.mean(confidence_scores.values()),
    "median": lambda confidence_scores: statistics.median(confidence_scores.values()),
    "min": lambda confidence_scores: min(confidence_scores.values()),
    "number": lambda confidence_scores: len(confidence_scores),
    "sum": lambda confidence_scores: sum(confidence_scores.values()),
    **{
        database: lambda confidence_scores: confidence_scores.get(database, 0.0)
        for database in {"BioGRID", "CORUM", "IntAct", "Reactome", "STRING"}
    },
}
