"""Mappings of configuration file entries to multiple testing corrections."""
from typing import Callable, Hashable
from analysis import correction

CORRECTION: dict[str,
                 Callable[[dict[Hashable, float]], dict[Hashable, float]]] = {
                     "Benjamini-Hochberg": correction.benjamini_hochberg,
                     "Bonferroni": correction.bonferroni
                 }
