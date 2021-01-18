import statistics

MERGE = {
    "max": max,
    "min": min,
    "max abs": lambda changes: max(abs(change) for change in changes),
    "min abs": lambda changes: min(abs(change) for change in changes),
    "mean": statistics.mean,
    "median": statistics.median,
}
