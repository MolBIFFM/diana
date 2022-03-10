def benjamini_hochberg(p_values):
    # Benjamini, Heller, Yekutieli (2009)
    m = len(p_values)
    sorted_p_values = sorted((list(item) for item in p_values.items()),
                             key=lambda item: item[1])

    for i in range(m - 1, 0, -1):
        if sorted_p_values[i - 1][1] / i > sorted_p_values[i][1] / (i + 1):
            sorted_p_values[i - 1][1] = sorted_p_values[i][1] * i / (i + 1)

    return {
        key: min(m * p_value / (i + 1), 1.0)
        for i, (key, p_value) in enumerate(sorted_p_values)
    }


def bonferroni(p_values):
    m = len(p_values)

    return {key: min(m * p_value, 1.0) for key, p_value in p_values.items()}
