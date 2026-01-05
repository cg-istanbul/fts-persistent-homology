############################################
# fts_basis.sage
############################################

import numpy as np
from itertools import combinations
from collections import defaultdict
from sage.graphs.graph import Graph

############################################
# 1. Distance matrix
############################################

def pairwise_distances(points):
    n = len(points)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            D[i, j] = np.linalg.norm(points[i] - points[j])
    return D

def trivial_topology_check(D, alpha):
    n = len(D)

    # smallest nonzero distance
    min_d = min(
        D[i][j] for i in range(n) for j in range(n) if i != j
    )
    max_d = max(
        D[i][j] for i in range(n) for j in range(n)
    )

    if alpha < min_d:
        return "discrete"

    if alpha >= max_d:
        return "indiscrete"

    return None


############################################
# 2. Maximal α-cliques (Rips-like open sets)
############################################

def rips_graph(points, alpha):
    n = len(points)
    G = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if np.linalg.norm(points[i] - points[j]) <= alpha:
                G[i].add(j)
                G[j].add(i)
    return G

def maximal_cliques_from_graph(G):
    edges = []
    for v in G:
        for u in G[v]:
            if u > v:
                edges.append((v,u))
    graph = Graph(edges)
    return [set(c) for c in graph.cliques_maximal()]

############################################
# 3. Minimal basis from open sets
############################################

def minimal_basis(open_sets, n):
    X = frozenset(range(n))
    open_sets = list(open_sets)
    if X not in open_sets:
        open_sets.append(X)

    basis = {}
    for x in range(n):
        containing = [U for U in open_sets if x in U]
        # now containing is NEVER empty
        Ux = set(containing[0])
        for U in containing[1:]:
            Ux &= set(U)
        basis[x] = frozenset(Ux)
    return basis


############################################
# 4. T0-identification
############################################

def t0_identification(basis):
    n = len(basis)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[ry] = rx

    for x in range(n):
        for y in range(n):
            if x in basis[y] and y in basis[x]:
                union(x, y)

    classes = defaultdict(set)
    for x in range(n):
        classes[find(x)].add(x)

    return list(classes.values())


############################################
# 5. Quotient minimal basis
############################################

def quotient_basis(basis, classes):
    class_index = {}
    for i, C in enumerate(classes):
        for x in C:
            class_index[x] = i

    new_basis = defaultdict(set)
    for x, Ux in basis.items():
        cx = class_index[x]
        for y in Ux:
            new_basis[cx].add(class_index[y])

    return {i: frozenset(U) for i, U in new_basis.items()}


############################################
# 6. Specialization poset
############################################

def specialization_order(basis):
    order = set()
    for i, Ui in basis.items():
        for j, Uj in basis.items():
            if Ui <= Uj:
                order.add((i, j))
    return order


############################################
# 7. Crosscut complex
############################################

def crosscut_2_skeleton(basis):
    atoms = list(basis.keys())
    vertices = [(i,) for i in atoms]
    edges = []
    triangles = []

    # edges
    for i, j in combinations(atoms, 2):
        if set(basis[i]) & set(basis[j]):
            edges.append((i, j))

    # triangles
    for i, j, k in combinations(atoms, 3):
        inter = set(basis[i]) & set(basis[j]) & set(basis[k])
        if inter:
            triangles.append((i, j, k))

    return vertices + edges + triangles


############################################

def fts_pipeline(X, alpha):
    """
    Full FTS pipeline up to crosscut complex.
    Returns:
        B0      : T0-reduced minimal basis
        poset   : specialization order
        crosscut : 2-skeleton of crosscut complex
    """
    D = pairwise_distances(X)
    case = trivial_topology_check(D, alpha)

    if case == "discrete":
        n = len(X)
        basis = {i: frozenset({i}) for i in range(n)}
        return basis, set(), [(i,) for i in range(n)]

    if case == "indiscrete":
        basis = {0: frozenset(range(len(X)))}
        return basis, {(0, 0)}, [(0,)]

    G = rips_graph(X, alpha)
    open_sets = maximal_cliques_from_graph(G)
    B = minimal_basis(open_sets, len(X))
    classes = t0_identification(B)
    B0 = quotient_basis(B, classes)
    poset = specialization_order(B0)
    crosscut = crosscut_2_skeleton(B0)

    return B0, poset, crosscut

