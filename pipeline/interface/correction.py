"""Mappings of configuration file entries to multiple testing corrections."""
from enrichment import correction

CORRECTION = {
    "Benjamini-Hochberg": correction.benjamini_hochberg,
    "Bonferroni": correction.bonferroni
}
