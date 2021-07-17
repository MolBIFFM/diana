def benjamini_hochberg(p_values):
    # Benjamini-Hochberg p-value adjustment (Benjamini, Heller, Yekutieli (2009))
    adj_p_values = sorted([[key, p_value]
                           for key, p_value in p_values.items()],
                          key=lambda item: item[1])
    m = len(adj_p_values)

    for i in range(m):
        adj_p_values[i][1] = min(
            min(m * adj_p_values[j][1] / (j + 1) for j in range(i, m)), 1.0)

    return {key: adj_p_value for key, adj_p_value in adj_p_values}
