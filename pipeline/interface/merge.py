import statistics

MERGE = {
    "max": max,
    "min": min,
    "maxabs": lambda changes: max(changes, key=abs),
    "minabs": lambda changes: min(changes, key=abs),
    "mean": statistics.mean,
    "median": statistics.median,
    "sum": sum,
}
