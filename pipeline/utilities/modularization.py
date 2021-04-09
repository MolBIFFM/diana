import networkx as nx


def louvain(G, weight="weight"):
    # Blondel, Guillaume, Lambiotte, Lefebvre (2008)
    name = list(G.nodes())
    communities = [[set([i]) for i in range(G.number_of_nodes())]]

    while True:
        A = nx.linalg.graphmatrix.adjacency_matrix(G, weight=weight)
        n = G.number_of_nodes()
        k = [A[i].sum() for i in range(n)]
        m = sum(k) / 2.0

        sigma_tot = [A[ci].sum() for ci in range(n)]
        sigma_in = [A[ci, ci] for ci in range(n)]

        k_in = A.toarray()

        community = [i for i in range(n)]
        communities.append([set([i]) for i in range(n)])

        community_aggregation = False
        modularity_optimization = True

        while modularity_optimization:
            modularity_optimization = False

            for i in range(n):
                sigma_tot[community[i]] -= A[i].sum()
                sigma_in[community[i]] -= sum([A[i, l] for l in communities[1][i]])

                k_in[i, community[i]] -= A[i, i]

                deltaQ = {}
                for j in range(n):
                    if A[i, j] and community[i] != community[j]:

                        k_in[j, community[i]] -= A[j, i]

                        deltaQ[j] = (
                            (
                                (sigma_in[community[j]] + k_in[i, community[j]])
                                / (2 * m)
                                - ((sigma_tot[community[j]] + k[i]) / (2 * m)) ** 2
                            )
                            - (
                                sigma_in[community[j]] / (2 * m)
                                - (sigma_tot[community[j]] / (2 * m)) ** 2
                                - (k[i] / (2 * m)) ** 2
                            )
                        ) - (
                            (
                                (sigma_in[community[i]] + k_in[i, community[i]])
                                / (2 * m)
                                - ((sigma_tot[community[i]] + k[i]) / (2 * m)) ** 2
                            )
                            - (
                                sigma_in[community[i]] / (2 * m)
                                - (sigma_tot[community[i]] / (2 * m)) ** 2
                                - (k[i] / (2 * m)) ** 2
                            )
                        )

                        k_in[j, community[i]] += A[j, i]

                if deltaQ:
                    max_j = max(deltaQ.items(), key=lambda item: item[1])[0]
                    if deltaQ[max_j] > 0.0:
                        for l in range(n):
                            if l != i:
                                k_in[l, community[i]] -= A[l, i]
                            k_in[l, community[max_j]] += A[l, i]

                        sigma_tot[community[max_j]] += A[i].sum()
                        sigma_in[community[max_j]] += sum(
                            [A[i, l] for l in communities[1][max_j]]
                        )

                        communities[1][community[i]].remove(i)
                        communities[1][community[max_j]].add(i)

                        community[i] = community[max_j]

                        modularity_optimization = True
                        community_aggregation = True

                    else:
                        k_in[i, community[i]] += A[i, i]
                        sigma_tot[community[i]] += A[i].sum()
                        sigma_in[community[i]] += sum(
                            [A[i, l] for l in communities[1][i]]
                        )

                else:
                    k_in[i, community[i]] += A[i, i]
                    sigma_tot[community[i]] += A[i].sum()
                    sigma_in[community[i]] += sum([A[i, l] for l in communities[1][i]])

        if community_aggregation:
            for ci in range(len(communities[1])):
                if communities[1][ci]:
                    communities[1][ci] = set.union(
                        *[
                            set(i for i in communities[0][node])
                            for node in communities[1][ci]
                        ]
                    )

            communities = communities[1:2]

            weights = [
                [0.0 for cj in range(len(communities[0]))]
                for ci in range(len(communities[0]))
            ]

            for i in range(n):
                for j in range(n):
                    weights[community[i]][community[j]] += k_in[i, community[j]]

            weights = [
                [
                    weights[ci][cj]
                    for cj in range(len(communities[0]))
                    if communities[0][cj]
                ]
                for ci in range(len(communities[0]))
                if communities[0][ci]
            ]

            communities[0] = [community for community in communities[0] if community]

            G = nx.Graph()
            G.add_nodes_from(range(len(communities[0])))

            for ci in range(len(communities[0])):
                if weights[ci][ci]:
                    G.add_edge(ci, ci, weight=2.0 * weights[ci][ci])
                for cj in range(ci):
                    if weights[ci][cj] or weights[cj][ci]:
                        G.add_edge(ci, cj, weight=weights[ci][cj] + weights[cj][ci])
        else:
            break

    return [set(name[node] for node in community) for community in communities[0]]
