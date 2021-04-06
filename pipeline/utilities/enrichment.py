def benjamini_hochberg(p_values):
    # Benjamini-Hochberg p-value adjustment (Yekutieli & Benjamini (1999))
    adj_p_values = sorted(
        [[key, p_value] for key, p_value in p_values.items()], key=lambda item: item[1]
    )

    for i in range(len(adj_p_values)):
        adj_p_values[i][1] = min(len(adj_p_values) * adj_p_values[i][1] / (i + 1), 1.0)
        for j in range(i, 0, -1):
            if adj_p_values[j - 1][1] <= adj_p_values[j][1]:
                break
            adj_p_values[j - 1][1] = adj_p_values[j][1]

    return {key: adj_p_value for key, adj_p_value in adj_p_values}
