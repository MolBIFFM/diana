"""Mappings of configuration file entries to multiple testing corrections."""
from typing import Callable, Hashable, Mapping

from analysis import correction

CORRECTION: dict[str, Callable[[Mapping[Hashable, float]],
                               dict[Hashable, float]]] = {
                                   "Benjamini-Hochberg":
                                       correction.benjamini_hochberg,
                                   "Bonferroni":
                                       correction.bonferroni,
                               }
