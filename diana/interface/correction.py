"""Mappings of configuration file entries to multiple testing corrections."""
from typing import Callable, Hashable, Mapping

from algorithms import correction

# Multiple testing correction procedures.
CORRECTION: dict[str, Callable[[Mapping[Hashable, float]],
                               dict[Hashable, float]]] = {
                                   "Benjamini-Hochberg":
                                       correction.benjamini_hochberg,
                                   "Benjamini-Yekutieli":
                                       correction.benjamini_yekutieli,
                                   "Holm":
                                       correction.holm,
                                   "Hommel":
                                       correction.hommel,
                               }
