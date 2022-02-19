def benjamini_hochberg(p_values):
    # Benjamini, Heller, Yekutieli (2009)
    m = len(p_values)
    adjusted_p_values = sorted([list(item) for item in p_values.items()],
                               key=lambda item: item[1])

    return {
        key:
        min(min(m * adjusted_p_values[j][1] / (j + 1) for j in range(i, m)),
            1.0)
        for i, (key, _) in enumerate(adjusted_p_values)
    }


def bonferroni(p_values):
    m = len(p_values)

    return {key: min(m * p_value, 1.0) for key, p_value in p_values.items()}
