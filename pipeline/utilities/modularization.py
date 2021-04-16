import networkx as nx


def clauset_newman_moore(G, weight="weight"):
    # Clauset, Newman, Moore (2004), Newman (2004)
    A = nx.linalg.graphmatrix.adjacency_matrix(G, weight=weight)
    n = G.number_of_nodes()
    k = [A[i].sum() for i in range(n)]
    m = sum(k) / 2.0
    a = [k[i] / (2.0 * m) for i in range(len(k))]

    communities = [[node] for node in G.nodes()]
    connected = [[bool(A[i, j]) for j in range(i)] for i in range(n)]

    deltaQ = [[0.0 for j in range(i)] for i in range(n)]
    for i in range(n):
        for j in range(i):
            if A[i, j]:
                deltaQ[i][j] = 1 / (2 * m) - k[i] * k[j] / ((2 * m) ** 2)

    max_entry = -1.0
    for i in range(n):
        for j in range(i):
            if deltaQ[i][j] > max_entry:
                max_entry = deltaQ[i][j]
                max_i = i
                max_j = j

    while deltaQ[max_i][max_j] > 0.0:
        communities[max_j].extend(communities[max_i])
        communities[max_i] = []

        deltaQprime = [deltaQ[i][:] for i in range(len(deltaQ))]

        for k in range(n):
            if communities[k] and k != max_i and k != max_j:
                if k < max_i and k < max_j:
                    if connected[max_i][k] and connected[max_j][k]:
                        deltaQprime[max_j][k] = deltaQ[max_i][k] + deltaQ[max_j][k]
                    elif connected[max_j][k]:
                        deltaQprime[max_j][k] = deltaQ[max_j][k] - 2 * a[max_i] * a[k]
                    elif connected[max_i][k]:
                        deltaQprime[max_j][k] = deltaQ[max_i][k] - 2 * a[max_j] * a[k]
                    connected[max_j][k] = connected[max_i][k] or connected[max_j][k]

                elif k < max_i:
                    if connected[max_i][k] and connected[k][max_j]:
                        deltaQprime[k][max_j] = deltaQ[max_i][k] + deltaQ[k][max_j]
                    elif connected[k][max_j]:
                        deltaQprime[k][max_j] = deltaQ[k][max_j] - 2 * a[max_i] * a[k]
                    elif connected[max_i][k]:
                        deltaQprime[k][max_j] = deltaQ[max_i][k] - 2 * a[max_j] * a[k]
                    connected[k][max_j] = connected[max_i][k] or connected[k][max_j]

                else:
                    if connected[k][max_i] and connected[k][max_j]:
                        deltaQprime[k][max_j] = deltaQ[k][max_i] + deltaQ[k][max_i]
                    elif connected[k][max_j]:
                        deltaQprime[k][max_j] = deltaQ[k][max_j] - 2 * a[max_i] * a[k]
                    elif connected[k][max_i]:
                        deltaQprime[k][max_j] = deltaQ[k][max_i] - 2 * a[max_j] * a[k]
                    connected[k][max_j] = connected[k][max_i] or connected[k][max_j]

        deltaQ = deltaQprime

        for i in range(n):
            if communities[i]:
                if i < max_i:
                    deltaQ[max_i][i] = 0.0
                    connected[max_i][i] = False

                elif i > max_i:
                    deltaQ[i][max_i] = 0.0
                    connected[i][max_i] = False

        a[max_j] += a[max_i]
        a[max_i] = 0.0

        max_entry = -1.0
        for i in range(n):
            if communities[i]:
                for j in range(i):
                    if communities[j]:
                        if deltaQ[i][j] > max_entry:
                            max_entry = deltaQ[i][j]
                            max_i = i
                            max_j = j

    return [set(community) for community in communities if community]


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

                k_in[i, community[i]] += A[i, i]

                if deltaQ:
                    max_j = max(deltaQ.items(), key=lambda item: item[1])[0]

                    if deltaQ[max_j] > 0.0:
                        sigma_tot[community[max_j]] += A[i].sum()
                        sigma_in[community[max_j]] += sum(
                            [A[i, l] for l in communities[1][max_j]]
                        )

                        for l in range(n):
                            k_in[l, community[i]] -= A[l, i]
                            k_in[l, community[max_j]] += A[l, i]

                        communities[1][community[i]].remove(i)
                        communities[1][community[max_j]].add(i)

                        community[i] = community[max_j]

                        modularity_optimization = True
                        community_aggregation = True

                    else:
                        sigma_tot[community[i]] += A[i].sum()
                        sigma_in[community[i]] += sum(
                            [A[i, l] for l in communities[1][i]]
                        )
                        k_in[i, community[i]] += A[i, i]

                else:
                    sigma_tot[community[i]] += A[i].sum()
                    sigma_in[community[i]] += sum([A[i, l] for l in communities[1][i]])
                    k_in[i, community[i]] += A[i, i]

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
