"""Community detection algorithms."""
import math
from typing import Hashable

import networkx as nx


def clauset_newman_moore(network: nx.Graph,
                         resolution: float = 1.0,
                         weight: str = "weight") -> list[set[Hashable]]:
    """
    Clauset-Newman-Moore community detection algorithm for undirected, weighted
    networks with parameterized modularity.

    A. Clauset, M. E. J. Newman and C. Moore, "Finding community structure
    in very large networks", Physical Review E, 2004.

    M. E. J. Newman, "Analysis of weighted networks", Physical Review E, 2004.

    M. E. J. Newman, "Equivalence between modularity optimization and maximum
    likelihood methods for community detection", Physical Review E, 2016.

    Args:
        network: An undirected, weighted graph.
        resolution: The resolution parameter for modularity.
        weight: The edge attribute to be utilized as edge weight.

    Returns:
        Communities of the network
    """
    A = nx.linalg.graphmatrix.adjacency_matrix(network, weight=weight)
    n = network.number_of_nodes()
    k = [math.fsum(A[i, l] for l in range(n)) for i in range(n)]
    m = math.fsum(k) / 2.0
    a = [k[i] / (2.0 * m) for i in range(n)]

    communities = [{node} for node in network.nodes()]
    connected = [[bool(A[i, j]) for j in range(i)] for i in range(n)]

    delta_q = [[0.0 for j in range(i)] for i in range(n)]
    for i in range(n):
        for j in range(i):
            if A[i, j]:
                delta_q[i][j] = 1 / (2.0 * m) - resolution * k[i] * k[j] / (
                    (2.0 * m)**2.0)

    max_entry = -1.0
    for i in range(n):
        for j in range(i):
            if delta_q[i][j] > max_entry:
                max_entry = delta_q[i][j]
                max_i = i
                max_j = j

    while delta_q[max_i][max_j] > 0.0:
        communities[max_j].update(communities[max_i])
        communities[max_i] = set()

        delta_q_prime = [delta_q[i][:] for i in range(len(delta_q))]

        for k in range(n):
            if communities[k] and k != max_i and k != max_j:
                if k < max_i and k < max_j:
                    if connected[max_i][k] and connected[max_j][k]:
                        delta_q_prime[max_j][
                            k] = delta_q[max_i][k] + delta_q[max_j][k]
                    elif connected[max_j][k]:
                        delta_q_prime[max_j][k] = delta_q[max_j][
                            k] - 2.0 * resolution * a[max_i] * a[k]
                    elif connected[max_i][k]:
                        delta_q_prime[max_j][k] = delta_q[max_i][
                            k] - 2.0 * resolution * a[max_j] * a[k]
                    connected[max_j][
                        k] = connected[max_i][k] or connected[max_j][k]

                elif k < max_i:
                    if connected[max_i][k] and connected[k][max_j]:
                        delta_q_prime[k][
                            max_j] = delta_q[max_i][k] + delta_q[k][max_j]
                    elif connected[k][max_j]:
                        delta_q_prime[k][max_j] = delta_q[k][
                            max_j] - 2.0 * resolution * a[max_i] * a[k]
                    elif connected[max_i][k]:
                        delta_q_prime[k][max_j] = delta_q[max_i][
                            k] - 2.0 * resolution * a[max_j] * a[k]
                    connected[k][
                        max_j] = connected[max_i][k] or connected[k][max_j]

                else:
                    if connected[k][max_i] and connected[k][max_j]:
                        delta_q_prime[k][
                            max_j] = delta_q[k][max_i] + delta_q[k][max_i]
                    elif connected[k][max_j]:
                        delta_q_prime[k][max_j] = delta_q[k][
                            max_j] - 2.0 * resolution * a[max_i] * a[k]
                    elif connected[k][max_i]:
                        delta_q_prime[k][max_j] = delta_q[k][
                            max_i] - 2.0 * resolution * a[max_j] * a[k]
                    connected[k][
                        max_j] = connected[k][max_i] or connected[k][max_j]

        delta_q = delta_q_prime

        for i in range(n):
            if communities[i]:
                if i < max_i:
                    delta_q[max_i][i] = 0.0
                    connected[max_i][i] = False

                elif i > max_i:
                    delta_q[i][max_i] = 0.0
                    connected[i][max_i] = False

        a[max_j] += a[max_i]
        a[max_i] = 0.0

        max_entry = -1.0
        for i in range(n):
            if communities[i]:
                for j in range(i):
                    if communities[j]:
                        if delta_q[i][j] > max_entry:
                            max_entry = delta_q[i][j]
                            max_i = i
                            max_j = j

    return [community for community in communities if community]


def louvain(network: nx.Graph,
            resolution: float = 1.0,
            weight: str = "weight") -> list[set[Hashable]]:
    """
    Louvain community detection algorithm for undirected, weighted networks with
    parameterized modularity.

    V. D. Blondel, J.-L. networkuillaume, R. Lambiotte and E. Lefebvre, "Fast
    unfolding of communities in large networks", Journal of Statistical
    Mechanics: Theory and Experiment, 2008.

    Args:
        network: An undirected, weighted graph.
        resolution: The resolution parameter for modularity.
        weight: The edge attribute to be utilized as edge weight.

    Returns:
        Communities of the network
    """
    name = list(network.nodes())
    communities = [[{i} for i in range(network.number_of_nodes())]]

    community_aggregation = True

    while community_aggregation:
        community_aggregation = False

        A = nx.linalg.graphmatrix.adjacency_matrix(network, weight=weight)
        n = network.number_of_nodes()
        k = [math.fsum(A[i, l] for l in range(n)) for i in range(n)]
        m = math.fsum(k) / 2.0

        sigma_tot = [A[ci].sum() for ci in range(n)]
        sigma_in = [A[ci, ci] for ci in range(n)]

        k_in = A.toarray()

        community = list(range(n))
        communities.append([{i} for i in range(n)])

        modularity_optimization = True

        while modularity_optimization:
            modularity_optimization = False

            for i in range(n):
                sigma_tot[community[i]] -= math.fsum(A[i, l] for l in range(n))
                sigma_in[community[i]] -= math.fsum(
                    [A[i, l] for l in range(n) if community[l] == community[i]])

                k_in[i, community[i]] -= A[i, i]

                delta_q = {}
                for j in range(n):
                    if A[i, j] and community[i] != community[j]:
                        delta_q[j] = (
                            ((sigma_in[community[j]] + k_in[i, community[j]]) /
                             (2 * m) - resolution *
                             ((sigma_tot[community[j]] + k[i]) / (2 * m))**2) -
                            (sigma_in[community[j]] / (2 * m) - resolution *
                             (sigma_tot[community[j]] /
                              (2 * m))**2 - resolution * (k[i] / (2 * m))**2)
                        ) - (((sigma_in[community[i]] + k_in[i, community[i]]) /
                              (2 * m) - resolution *
                              ((sigma_tot[community[i]] + k[i]) / (2 * m))**2) -
                             (sigma_in[community[i]] / (2 * m) - resolution *
                              (sigma_tot[community[i]] /
                               (2 * m))**2 - resolution * (k[i] / (2 * m))**2))

                sigma_tot[community[i]] += math.fsum(A[i, l] for l in range(n))
                sigma_in[community[i]] += math.fsum(
                    A[i, l] for l in range(n) if community[l] == community[i])

                k_in[i, community[i]] += A[i, i]

                if delta_q:
                    max_j = max(delta_q.items(), key=lambda item: item[1])[0]

                    if delta_q[max_j] > 0.0:
                        modularity_optimization = True
                        community_aggregation = True

                        sigma_tot[community[i]] -= math.fsum(
                            A[i, l] for l in range(n))
                        sigma_tot[community[max_j]] += math.fsum(
                            A[i, l] for l in range(n))

                        sigma_in[community[i]] -= math.fsum(
                            A[i, l]
                            for l in range(n)
                            if community[l] == community[i])
                        sigma_in[community[max_j]] += math.fsum(
                            A[i, l]
                            for l in range(n)
                            if community[l] == community[max_j])

                        for l in range(n):
                            k_in[l, community[i]] -= A[l, i]
                            k_in[l, community[max_j]] += A[l, i]

                        communities[1][community[i]].remove(i)
                        communities[1][community[max_j]].add(i)

                        community[i] = community[max_j]

        if community_aggregation:
            for ci in range(len(communities[1])):
                if communities[1][ci]:
                    communities[1][ci] = set.union(*[
                        set(i
                            for i in communities[0][node])
                        for node in communities[1][ci]
                    ])

            communities = communities[1:2]

            weights = [[0.0
                        for _ in range(len(communities[0]))]
                       for _ in range(len(communities[0]))]

            for i in range(n):
                weights[community[i]][community[i]] += A[i, i]

                for j in range(i):
                    weights[community[i]][community[j]] += A[i, j]
                    weights[community[j]][community[i]] += A[j, i]

            weights = [[
                weights[ci][cj]
                for cj in range(len(communities[0]))
                if communities[0][cj]
            ]
                       for ci in range(len(communities[0]))
                       if communities[0][ci]]

            communities[0] = [
                community for community in communities[0] if community
            ]

            network = nx.Graph()
            network.add_nodes_from(range(len(communities[0])))

            for ci in range(len(communities[0])):
                if weights[ci][ci]:
                    network.add_edge(ci, ci, weight=weights[ci][ci])

                for cj in range(ci):
                    if weights[ci][cj]:
                        network.add_edge(ci, cj, weight=weights[ci][cj])

    return [
        set(name[node] for node in community) for community in communities[0]
    ]
