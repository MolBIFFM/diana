import statistics

COMBINE = {
    "max": max,
    "min": min,
    "maxabs": lambda changes: max(changes, key=abs),
    "minabs": lambda changes: min(changes, key=abs),
    "mean": statistics.mean,
    "median": statistics.median,
    "sum": sum,
    "sumabs": lambda changes: sum(abs(change) for change in changes),
    "up": lambda changes: len([change for change in changes if change > 0.0])
    / len(changes),
    "down": lambda changes: len([change for change in changes if change < 0.0])
    / len(changes),
}
