import networkx as nx


def louvain(G, weight="weight"):
    # Blondel et al. (2008)
    A = nx.linalg.graphmatrix.adjacency_matrix(G, weight=weight)
    n = G.number_of_nodes()
    k = [A[i].sum() for i in range(n)]
    m = sum(k) / 2.0

    sigma_tot = {ci: k[ci] for ci in range(n)}
    sigma_in = {ci: A[ci, ci] for ci in range(n)}

    k_in = {ni: {cj: A[ni, cj] for cj in range(n)} for ni in range(n)}

    community = {i: i for i in range(n)}
    communities = [set([i]) for i in range(n)]

    change = False

    for i in range(n):
        deltaQ = {}
        for j in range(n):
            if A[i, j]:
                deltaQ[j] = (
                    (sigma_in[community[j]] + k_in[i][community[j]]) / (2 * m)
                    - ((sigma_tot[community[j]] + k[i]) / (2 * m)) ** 2
                ) - (
                    sigma_in[community[j]] / (2 * m)
                    - (sigma_tot[community[j]] / (2 * m)) ** 2
                    - (k[i] / (2 * m)) ** 2
                )

        if deltaQ:
            max_j = max(deltaQ.items(), key=lambda item: item[1])[0]
            if deltaQ[max_j] > 0.0:
                sigma_tot[community[i]] -= A[i].sum()
                sigma_tot[community[max_j]] += A[i].sum()

                sigma_in[community[i]] -= sum(
                    [A[i, l] for l in range(n) if community[l] == community[i]]
                )
                sigma_in[community[max_j]] += sum(
                    [A[i, l] for l in range(n) if community[l] == community[max_j]]
                )

                for n in k_in:
                    k_in[n][community[i]] -= A[n, i]
                    k_in[n][community[max_j]] += A[n, i]

                communities[community[i]].remove(i)
                communities[community[max_j]].add(i)

                community[i] = community[max_j]

                change = True

    deleted = 0
    for i in range(len(communities)):
        if not communities[i - deleted]:
            del communities[i - deleted]
            deleted += 1

    if change:
        H = nx.Graph()
        H.add_nodes_from(range(len(communities)))

        weights = {}
        for i in range(n):
            for j in range(n):
                if community[i] not in weights:
                    weights[community[i]] = {}
                if community[j] not in weights[community[i]]:
                    weights[community[i]][community[j]] = 0.0
                weights[community[i]][community[j]] += k_in[i][community[j]]

        for i in range(len(communities)):
            for j in range(len(communities)):
                if weights.get(i, {}).get(j, 0.0):
                    H.add_edge(i, j, weight=weights[i][j])

        
        expanded_communities = []
        for merged_community in louvain(H, weight=weight):
            if merged_community:
                expanded_communities.append(set())
                for node in merged_community:
                    expanded_communities[-1] |= communities[node]

        nodes = list(G.nodes())
        for i in range(len(expanded_communities)):
            expanded_communities[i] = set(nodes[n] for n in expanded_communities[i]) 

        return expanded_communities

    else:
        return communities
