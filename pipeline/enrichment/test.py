import scipy.stats


def binomial(k, M, n, N):
    return scipy.stats.binom.sf(k - 1, N, n / M)


def hypergeometric(k, M, n, N):
    return scipy.stats.hypergeom.sf(k - 1, M, n, N)
