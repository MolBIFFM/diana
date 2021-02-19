import networkx as nx


def louvain(G, weight="weight"):
    # Blondel et al. (2008)
    name = list(G.nodes())
    communities = [[set([i]) for i in range(G.number_of_nodes())]]

    change = True
    while change:
        A = nx.linalg.graphmatrix.adjacency_matrix(G, weight=weight)
        n = G.number_of_nodes()
        k = [A[i].sum() for i in range(n)]
        m = sum(k) / 2.0

        sigma_tot = [k[ci] for ci in range(n)]
        sigma_in = [A[ci, ci] for ci in range(n)]

        k_in = A.toarray()

        community = [i for i in range(n)]
        communities.append([set([i]) for i in range(n)])

        change = False

        for i in range(n):
            deltaQ = {}
            for j in range(n):
                if A[i, j]:
                    deltaQ[j] = (
                        (
                            (sigma_in[community[j]] + k_in[i, community[j]]) / (2 * m)
                            - ((sigma_tot[community[j]] + k[i]) / (2 * m)) ** 2
                        )
                        - (
                            sigma_in[community[j]] / (2 * m)
                            - (sigma_tot[community[j]] / (2 * m)) ** 2
                            - (k[i] / (2 * m)) ** 2
                        )
                    ) - (
                        (
                            (sigma_in[community[i]] + k_in[i, community[i]]) / (2 * m)
                            - ((sigma_tot[community[i]] + k[i]) / (2 * m)) ** 2
                        )
                        - (
                            sigma_in[community[i]] / (2 * m)
                            - (sigma_tot[community[i]] / (2 * m)) ** 2
                            - (k[i] / (2 * m)) ** 2
                        )
                    )

            if deltaQ:
                maxj = max(deltaQ.items(), key=lambda item: item[1])[0]
                if deltaQ[maxj] > 0.0:
                    sigma_tot[community[i]] -= A[i].sum()
                    sigma_tot[community[maxj]] += A[i].sum()

                    sigma_in[community[i]] -= sum(
                        [A[i, l] for l in range(n) if community[l] == community[i]]
                    )
                    sigma_in[community[maxj]] += sum(
                        [A[i, l] for l in range(n) if community[l] == community[maxj]]
                    )

                    for l in range(n):
                        k_in[l, community[i]] -= A[l, i]
                        k_in[l, community[maxj]] += A[l, i]

                    communities[1][community[i]].remove(i)
                    communities[1][community[maxj]].add(i)

                    community[i] = community[maxj]

                    change = True

        communities[1] = [community for community in communities[1] if community]

        if change:
            for ci in range(len(communities[1])):
                communities[1][ci] = set.union(
                    *[
                        set(i for i in communities[0][node])
                        for node in communities[1][ci]
                    ]
                )

            communities = communities[1:2]

            G = nx.Graph()
            G.add_nodes_from(range(len(communities[0])))

            weights = {}
            for i in range(n):
                for j in range(n):
                    if community[i] not in weights:
                        weights[community[i]] = {}
                    if community[j] not in weights[community[i]]:
                        weights[community[i]][community[j]] = 0.0
                    weights[community[i]][community[j]] += k_in[i, community[j]]

            for ci in range(len(communities[0])):
                for cj in range(len(communities[0])):
                    if weights.get(ci, {}).get(cj, 0.0):
                        G.add_edge(ci, cj, weight=weights[ci][cj])

    return [set(name[node] for node in community) for community in communities[0]]
