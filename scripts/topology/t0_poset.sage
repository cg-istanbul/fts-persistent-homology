# scripts/topology/t0_poset.sage

from sage.all import RR

def t0_poset_from_minimal_opens(U):
    n = len(U)

    key_to_cid = {}
    cls_of = [-1] * n
    rep_of = []
    U_of = []

    for x in range(n):
        ux = frozenset(U[x])
        cid = key_to_cid.get(ux, None)
        if cid is None:
            cid = len(rep_of)
            key_to_cid[ux] = cid
            rep_of.append(int(x))
            U_of.append(ux)
        cls_of[x] = cid

    n0 = len(rep_of)

    members_of = [[] for _ in range(n0)]
    for x in range(n):
        members_of[int(cls_of[x])].append(int(x))
    for c in range(n0):
        members_of[c].sort()
        
    leq = [[False] * n0 for _ in range(n0)]
    for a in range(n0):
        Ua = U_of[a]
        for b in range(n0):
            if Ua.issubset(U_of[b]):
                leq[a][b] = True

    E = []
    for a in range(n0):
        for b in range(n0):
            if a != b and leq[a][b]:
                E.append((a, b))

    return dict(
        n=int(n),
        n0=int(n0),
        cls_of=list(cls_of),
        rep_of=list(rep_of),
        U_of=list(U_of),
        members_of=members_of,
        leq=leq,
        E=E,
    )


def poset_class_map(P_i, P_j):
    rep_i = P_i["rep_of"]
    cls_j = P_j["cls_of"]
    n0_i = int(P_i["n0"])
    f = [0] * n0_i
    for c in range(n0_i):
        x = int(rep_i[c])
        f[c] = int(cls_j[x])
    return f


def check_poset_map_monotone(P_i, P_j, f, max_examples=5):
    leq_i = P_i.get("leq", None)
    leq_j = P_j.get("leq", None)
    if leq_i is None or leq_j is None:
        return False, {"reason": "missing leq matrix on one of the posets"}

    n0_i = int(P_i["n0"])
    bad = []
    for a in range(n0_i):
        fa = int(f[a])
        row = leq_i[a]
        for b in range(n0_i):
            if not row[b]:
                continue
            fb = int(f[b])
            if not leq_j[fa][fb]:
                bad.append((int(a), int(b), int(fa), int(fb)))
                if len(bad) >= int(max_examples):
                    return False, {"bad_count_at_least": len(bad), "examples": bad}

    return True, {"bad_count": 0, "examples": []}
