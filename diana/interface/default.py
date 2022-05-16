"""Mappings of configuration file entries to dependent default settings."""

from typing import Optional

MEASUREMENT_RANGE: dict[Optional[str], tuple[float, float]] = {
    None: (-1.0, 1.0),
    "percentile": (2.5, 97.5),
    "quantile": (0.025, 0.975),
    "standard score": (-2.0, 2.0)
}
