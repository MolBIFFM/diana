"""Mappings of configuration file entries to functions to define order."""
import math
from typing import Callable, Union

# Functions determining order of distinct modification sites.
SITE_ORDER: dict[str, Callable[[tuple[int, float]], Union[int, float]]] = {
    "absolute measurement": lambda site: abs(math.log2(site[1])),
    "measurement": lambda site: site[1],
    "position": lambda site: site[0]
}
