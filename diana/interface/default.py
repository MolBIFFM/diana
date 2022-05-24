"""Mappings of configuration file entries to dependent default settings."""

from typing import Optional

MEASUREMENT_RANGE: dict[Optional[str], tuple[float, float]] = {
    None: (-1.0, 1.0),
    "ratio": (0.5, 2.0),
    "percentile": (25.0, 75.0),
    "quantile": (0.25, 0.75),
    "standard score": (-1.0, 1.0)
}
