# scripts/util/misc.sage

def distance_matrix_from_points(points):
    n = len(points)
    D = matrix(RR, n, n)
    for i in range(n):
        D[i,i] = RR(0)
        xi, yi = points[i]
        for j in range(i+1, n):
            xj, yj = points[j]
            d = sqrt((xi-xj)^2 + (yi-yj)^2)
            D[i,j] = d
            D[j,i] = d
    return D


def distance_matrix_from_array(X):
    n = int(X.shape[0])
    D = matrix(RR, n, n)
    for i in range(n):
        D[i,i] = RR(0)
        xi = X[i]
        for j in range(i+1, n):
            dij = float(np.linalg.norm(xi - X[j]))
            D[i,j] = RR(dij)
            D[j,i] = RR(dij)
    return D


def n_points(D):
    return int(D.nrows()) if hasattr(D, 'nrows') else int(len(D))


def dist(D, i, j):
    if hasattr(D, 'nrows'):
        return RR(D[int(i), int(j)])
    return RR(D[int(i)][int(j)])


def all_points_set(n):
    return frozenset(range(int(n)))


def _sorted_neighbor_indices(D, i):
    n = n_points(D)
    i = int(i)
    arr = [(dist(D, i, j), j) for j in range(n) if j != i]
    arr.sort(key=lambda t: (t[0], t[1]))
    return [j for (_, j) in arr]


def knn_set(D, i, k):
    n = n_points(D)
    k = max(0, min(int(k), n - 1))
    i = int(i)
    if k == 0:
        return frozenset([i])
    nbrs = _sorted_neighbor_indices(D, i)[:k]
    return frozenset([i] + [int(j) for j in nbrs])


def knn_radius(D, i, k_density):
    n = n_points(D)
    if n <= 1:
        return None
    k = max(1, min(int(k_density), n - 1))
    i = int(i)
    nbrs = _sorted_neighbor_indices(D, i)
    return dist(D, i, nbrs[k - 1])


def point_density_knn_radius(D, k_density=7, eps=1e-12):
    n = n_points(D)
    r = [None] * n
    rho = [RR(0)] * n
    for i in range(n):
        ri = knn_radius(D, i, k_density)
        r[i] = ri
        rho[i] = RR(0) if ri is None else RR(1) / (RR(ri) + RR(eps))
    return rho, r


def discrete_minimal_opens(n):
    n = int(n)
    return [frozenset([i]) for i in range(n)]


def minimal_nonempty_intersections(C, n):
    n = int(n)
    X = frozenset(range(n))

    G = []
    seen = set()
    for S in C:
        if S is None:
            continue
        fs = frozenset(S)
        if len(fs) == 0:
            continue
        if fs not in seen:
            seen.add(fs)
            G.append(fs)

    if X not in seen:
        G.append(X)

    contains = [[] for _ in range(n)]
    for fs in G:
        for x in fs:
            xi = int(x)
            if 0 <= xi < n:
                contains[xi].append(fs)

    U = [None] * n
    for x in range(n):
        lst = contains[x]
        inter = lst[0]
        for fs in lst[1:]:
            inter = inter.intersection(fs)
            if len(inter) == 0:
                break
        U[x] = X if len(inter) == 0 else frozenset(inter)

    return U


def distinct_minimal_opens(U, drop_X=True):
    if U is None or len(U) == 0:
        return []
    n = len(U)
    X = frozenset(range(n))
    seen = set()
    V = []
    for S in U:
        fs = frozenset(S)
        if drop_X and fs == X:
            continue
        if fs not in seen:
            seen.add(fs)
            V.append(fs)
    return V


def continuity_ok(prev_U, curr_U, max_examples=3):
    if prev_U is None or curr_U is None:
        return True, 0, []
    n = min(len(prev_U), len(curr_U))
    bad = 0
    ex = []
    for x in range(n):
        union_prev = set()
        for y in curr_U[x]:
            union_prev.update(prev_U[int(y)])
        union_prev = frozenset(union_prev)
        if union_prev != curr_U[x]:
            bad += 1
            if len(ex) < int(max_examples):
                ex.append((x, len(curr_U[x]), len(union_prev)))
    return (bad == 0), bad, ex
