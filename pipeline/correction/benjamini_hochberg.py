def benjamini_hochberg(p_values):
    # Benjamini-Hochberg p-value adjustment (Benjamini, Heller, Yekutieli (2009))
    adjusted_p_values = sorted([list(item) for item in p_values.items()],
                               key=lambda item: item[1])
    m = len(adjusted_p_values)

    for i in range(m):
        adjusted_p_values[i][1] = min(
            min(m * adjusted_p_values[j][1] / (j + 1) for j in range(i, m)),
            1.0)

    return {
        key: adjusted_p_value
        for key, adjusted_p_value in adjusted_p_values
    }
