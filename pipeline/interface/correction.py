from enrichment import correction

CORRECTION = {
    "Benjamini-Hochberg": correction.benjamini_hochberg,
    "Bonferroni": correction.bonferroni
}
