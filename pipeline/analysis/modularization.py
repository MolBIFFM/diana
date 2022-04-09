"""Community detection algorithms."""
from typing import Hashable

import networkx as nx


def clauset_newman_moore(network: nx.Graph,
                         resolution: float = 1.0,
                         weight: str = "weight") -> list[set[Hashable]]:
    """
    Clauset-Newman-Moore community detection algorithm for undirected, weighted
    networks with parameterized modularity.

    Clauset, A. et al. (2004) Finding community structure in very large
        networks, Phys. Rev. E, 70.

    Newman, M. E. J. (2004) Analysis of weighted networks, Phys. Rev. E, 70.

    Newman, M. E. J. (2016) Equivalence between modularity optimization and
        maximum likelihood methods for community detection, Phys. Rev. E, 94.

    Args:
        network: An undirected, weighted graph.
        resolution: The resolution parameter for modularity.
        weight: The edge attribute to be utilized as edge weight.

    Returns:
        The communities of the network.
    """
    communities = {i: {i} for i in network.nodes()}

    A = {
        i: {j: network.edges[i, j][weight] for j in nx.neighbors(network, i)
           } for i in network.nodes()
    }

    connected = {i: set(A[i]) - {i} for i in A}

    k = {i: sum(A[i].values()) for i in A}
    m = sum(k.values()) / 2.0

    delta_q = {
        i: {
            j: (1.0 / (2.0 * m) - resolution * k[i] * k[j] /
                ((2.0 * m)**2.0) if j in A[i] else 0.0) for j in A if j < i
        } for i in A
    }

    a = {i: k[i] / (2.0 * m) for i in A}

    max_entry = -1.0
    for i in delta_q:
        for j in delta_q[i]:
            if delta_q[i][j] > max_entry:
                max_entry = delta_q[i][j]
                max_i, max_j = i, j

    while delta_q[max_i][max_j] > 0.0:
        communities[max_j].update(communities[max_i])
        del communities[max_i]

        delta_q_prime = {
            i: {j: delta_q[i][j] for j in delta_q[i]} for i in delta_q
        }

        for n in connected[max_i] & connected[max_j]:
            if n < max_j:
                delta_q_prime[max_j][n] = delta_q[max_i][n] + delta_q[max_j][n]
            elif n < max_i:
                delta_q_prime[n][max_j] = delta_q[max_i][n] + delta_q[n][max_j]
            else:
                delta_q_prime[n][max_j] = delta_q[n][max_i] + delta_q[n][max_j]

        for n in connected[max_i] - connected[max_j] - {max_j}:
            if n < max_j:
                delta_q_prime[max_j][
                    n] = delta_q[max_i][n] - 2.0 * resolution * a[max_j] * a[n]
            elif n < max_i:
                delta_q_prime[n][max_j] = delta_q[max_i][
                    n] - 2.0 * resolution * a[max_j] * a[n]
            else:
                delta_q_prime[n][max_j] = delta_q[n][
                    max_i] - 2.0 * resolution * a[max_j] * a[n]

            connected[max_j].add(n)
            connected[n].add(max_j)

        for n in connected[max_j] - connected[max_i] - {max_i}:
            if n < max_j:
                delta_q_prime[max_j][
                    n] = delta_q[max_j][n] - 2.0 * resolution * a[max_i] * a[n]
            else:
                delta_q_prime[n][max_j] = delta_q[n][
                    max_j] - 2.0 * resolution * a[max_i] * a[n]

        delta_q = delta_q_prime

        for n in list(delta_q):
            if n > max_i:
                del delta_q[n][max_i]
            connected[n].discard(max_i)
        del delta_q[max_i]

        a[max_j] += a[max_i]

        max_entry = -1.0
        for i in delta_q:
            for j in delta_q[i]:
                if delta_q[i][j] > max_entry:
                    max_entry = delta_q[i][j]
                    max_i, max_j = i, j

    return list(communities.values())


def louvain(network: nx.Graph,
            resolution: float = 1.0,
            weight: str = "weight") -> list[set[Hashable]]:
    """
    Louvain community detection algorithm for undirected, weighted networks with
    parameterized modularity.

    Blondel, V. D. et al. (2008) Fast unfolding of communities in large
        networks, J. Stat. Mech.: Theory Exp.

    Clauset, A. et al. (2004) Finding community structure in very large
        networks, Phys. Rev. E, 70.

    Newman, M. E. J. (2016) Equivalence between modularity optimization and
        maximum likelihood methods for community detection, Phys. Rev. E, 94.

    Args:
        network: An undirected, weighted graph.
        resolution: The resolution parameter for modularity.
        weight: The edge attribute to be used as edge weight.

    Returns:
        The communities of the network.
    """
    communities = [{i: {i} for i in network.nodes()}]

    community_aggregation = True
    while community_aggregation:
        community_aggregation = False

        community = {i: i for i in network.nodes()}
        communities.append({i: {i} for i in network.nodes()})

        A = {
            i:
            {j: network.edges[i, j][weight] for j in nx.neighbors(network, i)}
            for i in network.nodes()
        }

        sigma_in = {i: A[i].get(i, 0.0) for i in network.nodes()}
        sigma_tot = {i: sum(A[i].values()) for i in network.nodes()}

        k = {i: sum(A[i].values()) for i in A}

        k_in = {
            i: {
                **{
                    community[i]: 0.0
                },
                **{community[j]: A[i][j] for j in A[i]}
            } for i in A
        }

        m = sum(k.values()) / 2.0

        modularity_optimization = True

        while modularity_optimization:
            modularity_optimization = False

            for i in A:
                for n in A[i]:
                    sigma_tot[community[i]] -= A[i][n]

                    if community[n] == community[i]:
                        sigma_in[community[i]] -= A[i][n]

                k_in[i][community[i]] -= A[i].get(i, 0.0)

                delta_q = {}
                for j in A[i]:
                    if community[i] != community[j]:
                        delta_q[j] = (
                            ((sigma_in[community[j]] + k_in[i][community[j]]) /
                             (2.0 * m) - resolution *
                             ((sigma_tot[community[j]] + k[i]) /
                              (2.0 * m))**2.0) -
                            (sigma_in[community[j]] / (2.0 * m) - resolution *
                             (sigma_tot[community[j]] /
                              (2.0 * m))**2.0 - resolution * (k[i] /
                                                              (2.0 * m))**2.0)
                        ) - (((sigma_in[community[i]] + k_in[i][community[i]]) /
                              (2.0 * m) - resolution *
                              ((sigma_tot[community[i]] + k[i]) /
                               (2.0 * m))**2.0) -
                             (sigma_in[community[i]] / (2.0 * m) - resolution *
                              (sigma_tot[community[i]] /
                               (2.0 * m))**2.0 - resolution * (k[i] /
                                                               (2.0 * m))**2.0))

                for n in A[i]:
                    sigma_tot[community[i]] += A[i][n]

                    if community[n] == community[i]:
                        sigma_in[community[i]] += A[i][n]

                k_in[i][community[i]] += A[i].get(i, 0.0)

                if delta_q:
                    max_j = max(delta_q.items(), key=lambda item: item[1])[0]

                    if delta_q[max_j] > 0.0:
                        modularity_optimization = True
                        community_aggregation = True

                        for n in A[i]:
                            sigma_tot[community[i]] -= A[i][n]
                            sigma_tot[community[max_j]] += A[i][n]

                            if community[n] == community[i]:
                                sigma_in[community[i]] -= A[i][n]
                            elif community[n] == community[max_j]:
                                sigma_in[community[max_j]] += A[i][n]

                            k_in[n][community[i]] -= A[n][i]

                            if community[max_j] not in k_in[n]:
                                k_in[n][community[max_j]] = 0.0

                            k_in[n][community[max_j]] += A[n][i]

                        communities[1][community[i]].remove(i)
                        communities[1][community[max_j]].add(i)

                        community[i] = community[max_j]

        if community_aggregation:
            for ci in communities[1]:
                if communities[1][ci]:
                    communities[1][ci] = set.union(
                        *(communities[0][node] for node in communities[1][ci]))
            communities = communities[1:2]

            weights = {ci: {} for ci in community.values()}
            for i in A:
                for cj in k_in[i]:
                    if cj not in weights[community[i]]:
                        weights[community[i]][cj] = 0.0
                    weights[community[i]][cj] += k_in[i][cj]

            network = nx.Graph()
            for ci in weights:
                if communities[0][ci]:
                    network.add_node(ci)

                if weights[ci][ci]:
                    network.add_edge(ci, ci, weight=weights[ci][ci])

                for cj in weights[ci]:
                    if weights[ci][cj]:
                        network.add_edge(ci, cj, weight=weights[ci][cj])

    return [community for community in communities[0].values() if community]
