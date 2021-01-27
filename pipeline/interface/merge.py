import statistics

MERGE = {
    "max": max,
    "min": min,
    "max abs": lambda changes: max(changes, key=abs),
    "min abs": lambda changes: min(changes, key=abs),
    "mean": statistics.mean,
    "median": statistics.median,
}
