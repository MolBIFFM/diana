"""Mappings of configuration file entries to dependent default settings."""

CHANGE_RANGE: dict[str, tuple[float, float]] = {
    None: (-1.0, 1.0),
    "quantile": (0.025, 0.975),
    "standard score": (-2.0, 2.0)
}