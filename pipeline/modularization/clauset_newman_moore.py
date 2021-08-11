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
                deltaQ[i][j] = 1 / (2 * m) - k[i] * k[j] / ((2 * m)**2)

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
                        deltaQprime[max_j][
                            k] = deltaQ[max_i][k] + deltaQ[max_j][k]
                    elif connected[max_j][k]:
                        deltaQprime[max_j][
                            k] = deltaQ[max_j][k] - 2 * a[max_i] * a[k]
                    elif connected[max_i][k]:
                        deltaQprime[max_j][
                            k] = deltaQ[max_i][k] - 2 * a[max_j] * a[k]
                    connected[max_j][
                        k] = connected[max_i][k] or connected[max_j][k]

                elif k < max_i:
                    if connected[max_i][k] and connected[k][max_j]:
                        deltaQprime[k][
                            max_j] = deltaQ[max_i][k] + deltaQ[k][max_j]
                    elif connected[k][max_j]:
                        deltaQprime[k][
                            max_j] = deltaQ[k][max_j] - 2 * a[max_i] * a[k]
                    elif connected[max_i][k]:
                        deltaQprime[k][
                            max_j] = deltaQ[max_i][k] - 2 * a[max_j] * a[k]
                    connected[k][
                        max_j] = connected[max_i][k] or connected[k][max_j]

                else:
                    if connected[k][max_i] and connected[k][max_j]:
                        deltaQprime[k][
                            max_j] = deltaQ[k][max_i] + deltaQ[k][max_i]
                    elif connected[k][max_j]:
                        deltaQprime[k][
                            max_j] = deltaQ[k][max_j] - 2 * a[max_i] * a[k]
                    elif connected[k][max_i]:
                        deltaQprime[k][
                            max_j] = deltaQ[k][max_i] - 2 * a[max_j] * a[k]
                    connected[k][
                        max_j] = connected[k][max_i] or connected[k][max_j]

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
