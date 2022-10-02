"""Mappings of configuration file entries to functions for prioritization."""

from typing import Callable

SITE_PRIORITIZATION: dict[str, Callable[[float], float]] = {
    "abs": abs,
    "increase": lambda site: site,
    "decrease": lambda site: -site
}
