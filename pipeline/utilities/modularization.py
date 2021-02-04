import networkx as nx
import numpy as np


def clauset_newman_moore(G, weight="weight"):
    # Clauset, Newman, Moore (2004), Newman (2004)
    A = nx.linalg.graphmatrix.adjacency_matrix(G, weight=weight).toarray()
    n = G.number_of_nodes()
    m = np.concatenate(A).sum() / 2.0
    k = np.sum(A, axis=0)
    a = k / (2.0 * m)

    deltaQ = np.zeros_like(A, dtype=float)
    connected = np.array(A, dtype=bool)
    for i in range(n):
        for j in range(i):
            if A[i, j]:
                deltaQ[i, j] = 1 / (2 * m) - k[i] * k[j] / ((2 * m) ** 2)

    communities = [[node] for node in G.nodes()]

    i, j = np.unravel_index(np.argmax(deltaQ), deltaQ.shape)
    while deltaQ[i, j] > 0.0:
        deltaQprime = np.copy(deltaQ)
        for k in range(n):
            if k != i and k != j:
                if k < i and k < j:
                    if connected[i, k] and connected[j, k]:
                        deltaQprime[j, k] = deltaQ[i, k] + deltaQ[j, k]
                    elif connected[j, k]:
                        deltaQprime[j, k] = deltaQ[j, k] - 2 * a[i] * a[k]
                    elif connected[i, k]:
                        deltaQprime[j, k] = deltaQ[i, k] - 2 * a[j] * a[k]
                    connected[j, k] = connected[i, k] or connected[j, k]

                elif k < i:
                    if connected[i, k] and connected[k, j]:
                        deltaQprime[k, j] = deltaQ[i, k] + deltaQ[k, j]
                    elif connected[k, j]:
                        deltaQprime[k, j] = deltaQ[k, j] - 2 * a[i] * a[k]
                    elif connected[i, k]:
                        deltaQprime[k, j] = deltaQ[i, k] - 2 * a[j] * a[k]
                    connected[k, j] = connected[i, k] or connected[k, j]

                else:
                    if connected[k, i] and connected[k, j]:
                        deltaQprime[k, j] = deltaQ[k, i] + deltaQ[k, i]
                    elif connected[k, j]:
                        deltaQprime[k, j] = deltaQ[k, j] - 2 * a[i] * a[k]
                    elif connected[k, i]:
                        deltaQprime[k, j] = deltaQ[k, i] - 2 * a[j] * a[k]
                    connected[k, j] = connected[k, i] or connected[k, j]

        deltaQprime[i, :] = 0.0
        deltaQprime[:, i] = 0.0
        deltaQ = deltaQprime

        a[j] += a[i]
        a[i] = 0

        connected[i, :] = False
        connected[:, i] = False

        communities[j].extend(communities[i])
        communities[i] = []

        i, j = np.unravel_index(np.argmax(deltaQ), deltaQ.shape)

    return [frozenset(community) for community in communities if community]
