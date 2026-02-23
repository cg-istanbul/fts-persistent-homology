# scripts/topology/order_complex.sage

def order_complex_2skeleton(poset):
    n0 = int(poset["n0"])
    leq = poset["leq"]

    V = list(range(n0))

    lt = [[False]*n0 for _ in range(n0)]
    for a in range(n0):
        for b in range(n0):
            if a != b and leq[a][b] and (not leq[b][a]):
                lt[a][b] = True

    E = []
    out_nbr = [set() for _ in range(n0)]
    for a in range(n0):
        for b in range(n0):
            if lt[a][b]:
                E.append((a, b) if a < b else (b, a))
                out_nbr[a].add(b)

    E = sorted(list(set(E)))

    T = set()
    for a in range(n0):
        for b in out_nbr[a]:
            for c in out_nbr[b]:
                if c in out_nbr[a]:
                    tri = tuple(sorted((a, b, c)))
                    if len(set(tri)) == 3:
                        T.add(tri)

    T = sorted(list(T))

    tri_edges = set()
    for (a, b, c) in T:
        aa, bb, cc = sorted((int(a), int(b), int(c)))
        tri_edges.add((aa, bb))
        tri_edges.add((aa, cc))
        tri_edges.add((bb, cc))

    E_set = set((min(int(u), int(v)), max(int(u), int(v))) for (u, v) in E if int(u) != int(v))
    E_set |= tri_edges
    E = sorted(E_set)

    return dict(V=V, E=E, T=T)
