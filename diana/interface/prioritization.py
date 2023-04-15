"""Mappings of configuration file entries to functions for prioritization."""
import math
from typing import Callable

# Functions determining prioritization of distinct modification sites.
SITE_PRIORITIZATION: dict[str, Callable[[float], float]] = {
    "absolute": lambda site: abs(math.log2(site)),
    "increase": lambda site: site,
    "decrease": lambda site: -site
}
