import networkx as nx
import scipy.sparse

def clauset_newman_moore(G, weight="weight"):
    # Clauset, Newman, Moore (2004), Newman (2004)
    A = nx.linalg.graphmatrix.adjacency_matrix(G, weight=weight).toarray().tolist()
    n = G.number_of_nodes()
    k = [sum(row) for row in A]
    m = sum(k) / 2.0
    a = [k[i] / (2.0 * m) for i in range(len(k))]
    
    communities = [[node] for node in G.nodes()]
    connected = [[bool(A[i][j]) for j in range(i)] for i in range(n)]

    deltaQ = [[0.0 for j in range(i)] for i in range(n)]
    for i in range(n):
        for j in range(i):
            if A[i][j]:
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
            if k != i and k != j:
                if k < i and k < j:
                    if connected[i][k] and connected[j][k]:
                        deltaQprime[j][k] = deltaQ[i][k] + deltaQ[j][k]
                    elif connected[j][k]:
                        deltaQprime[j][k] = deltaQ[j][k] - 2 * a[i] * a[k]
                    elif connected[i][k]:
                        deltaQprime[j][k] = deltaQ[i][k] - 2 * a[j] * a[k]
                    connected[j][k] = connected[i][k] or connected[j][k]

                elif k < i:
                    if connected[i][k] and connected[k][j]:
                        deltaQprime[k][j] = deltaQ[i][k] + deltaQ[k][j]
                    elif connected[k][j]:
                        deltaQprime[k][j] = deltaQ[k][j] - 2 * a[i] * a[k]
                    elif connected[i][k]:
                        deltaQprime[k][j] = deltaQ[i][k] - 2 * a[j] * a[k]
                    connected[k][j] = connected[i][k] or connected[k][j]

                else:
                    if connected[k][i] and connected[k][j]:
                        deltaQprime[k][j] = deltaQ[k][i] + deltaQ[k][i]
                    elif connected[k][j]:
                        deltaQprime[k][j] = deltaQ[k][j] - 2 * a[i] * a[k]
                    elif connected[k][i]:
                        deltaQprime[k][j] = deltaQ[k][i] - 2 * a[j] * a[k]
                    connected[k][j] = connected[k][i] or connected[k][j]

        deltaQ = deltaQprime

        for i in range(n):
            if i < max_i:
                deltaQ[max_i][i] = 0.0
                connected[max_i][i] = False

            elif i > max_i:
                deltaQ[i][max_i] = 0.0
                connected[i][max_i] = False

        a[j] += a[i]
        a[i] = 0.0

        max_entry = -1.0
        for i in range(n):
            for j in range(i):
                if deltaQ[i][j] > max_entry:
                    max_entry = deltaQ[i][j]
                    max_i = i
                    max_j = j

    return [frozenset(community) for community in communities if community]
