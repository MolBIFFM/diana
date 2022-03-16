"""Mappings of bases to combination support logarithmic representation."""

import math
from typing import Callable, Collection

LOGARITHM: Callable[[Collection[float]], float] = {
    None:
        lambda changes, combination: math.log2(combination(changes)),
    2:
        lambda changes, combination: math.log2(
            combination([math.pow(2.0, change) for change in changes])),
    10:
        lambda changes, combination: math.log2(
            combination([math.pow(10.0, change) for change in changes])),
}
