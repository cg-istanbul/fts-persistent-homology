# fts_pipeline.sage

from sage.all import SimplicialComplex, RR, floor, QQ
import os

# ----------------------------
# Basic helpers
# ----------------------------

def n_points(D):
    return int(D.nrows()) if hasattr(D, "nrows") else int(len(D))

def dist(D, i, j):
    if hasattr(D, "nrows"):
        return RR(D[i, j])
    return RR(D[i][j])

def all_points_set(n):
    return frozenset(range(n))

def _sorted_neighbor_indices(D, i):
    n = n_points(D)
    arr = [(dist(D, i, j), j) for j in range(n) if j != i]
    arr.sort(key=lambda t: (t[0], t[1]))
    return [j for (_, j) in arr]

def knn_set(D, i, k):
    n = n_points(D)
    k = max(0, min(int(k), n-1))
    if k == 0:
        return frozenset([i])
    nbrs = _sorted_neighbor_indices(D, i)[:k]
    return frozenset([i] + nbrs)

def knn_radius(D, i, k_density):
    n = n_points(D)
    if n <= 1:
        return None
    k = max(1, min(int(k_density), n-1))
    nbrs = _sorted_neighbor_indices(D, i)
    return dist(D, i, nbrs[k-1])

def point_density_knn_radius(D, k_density=7, eps=1e-12):
    """
    rho(i) = 1 / (r_k(i) + eps)
    """
    n = n_points(D)
    r = [None]*n
    rho = [RR(0)]*n
    for i in range(n):
        ri = knn_radius(D, i, k_density)
        r[i] = ri
        if ri is None:
            rho[i] = RR(0)
        else:
            rho[i] = RR(1) / (RR(ri) + RR(eps))
    return rho, r

# ----------------------------
# Quantile utilities
# ----------------------------

def _quantile_value(sorted_vals, q):
    n = len(sorted_vals)
    if n == 0:
        return RR(0)
    if q <= 0:
        return RR(sorted_vals[0])
    if q >= 1:
        return RR(sorted_vals[-1])
    idx = int(floor(RR(q) * RR(n-1)))
    idx = max(0, min(idx, n-1))
    return RR(sorted_vals[idx])

def taus_from_rho_quantiles(rho, m, q_start=0.02, q_end=0.90, mode="linear_in_q"):
    assert m >= 2
    q_start = float(q_start); q_end = float(q_end)
    q_start = max(0.0, min(q_start, 1.0))
    q_end   = max(0.0, min(q_end,   1.0))
    if q_end < q_start:
        q_start, q_end = q_end, q_start

    vals = sorted([RR(x) for x in rho])
    t_start = _quantile_value(vals, q_start)
    t_end   = _quantile_value(vals, q_end)

    taus = []
    if mode == "linear_in_tau":
        if t_end == t_start:
            return [t_start for _ in range(m)]
        for i in range(m):
            a = RR(i) / RR(m-1)
            taus.append(t_start + a*(t_end - t_start))
    else:
        if m == 1:
            return [t_start]
        for i in range(m):
            a = float(i) / float(m-1)
            q = q_start + a*(q_end - q_start)
            taus.append(_quantile_value(vals, q))

    for i in range(1, m):
        if taus[i] < taus[i-1]:
            taus[i] = taus[i-1]

    uniq = []
    for t in taus:
        if len(uniq) == 0 or t != uniq[-1]:
            uniq.append(t)

    if len(uniq) >= m:
        return uniq[:m]

    if t_end == t_start:
        return [t_start for _ in range(m)]
    taus2 = []
    for i in range(m):
        a = RR(i) / RR(m-1)
        taus2.append(t_start + a*(t_end - t_start))
    return taus2


def taus_linear(tau_min, tau_max, m):
    assert m >= 2
    t0 = RR(tau_min)
    t1 = RR(tau_max)
    if t1 < t0:
        t0, t1 = t1, t0
    if t1 == t0:
        return [t0 for _ in range(m)]
    taus = []
    for i in range(int(m)):
        a = RR(i) / RR(int(m)-1)
        taus.append(t0 + a*(t1 - t0))
    return taus

# ----------------------------
# Anchor neighborhoods
# ----------------------------

def anchors_at_threshold(rho, tau):
    tau = RR(tau)
    return set([i for i,ri in enumerate(rho) if RR(ri) >= tau])

def anchor_neighborhood(D, a, r_k_a, lambda_scale=1.5):
    n = n_points(D)
    if r_k_a is None:
        return frozenset([a])
    rad = RR(lambda_scale) * RR(r_k_a)
    return frozenset([x for x in range(n) if dist(D, a, x) <= rad])

# ----------------------------
# Minimal basis from generators
# ----------------------------

def minimal_basis_from_generators(generators, n, include_X=False):
    X = all_points_set(n)
    gens = [frozenset(G) for G in generators if G is not None]
    if include_X:
        gens.append(X)

    U = []
    for x in range(n):
        containing = [G for G in gens if x in G]
        if len(containing) == 0:
            U.append(X)
        else:
            inter = set(containing[0])
            for G in containing[1:]:
                inter.intersection_update(G)
            U.append(frozenset(inter) if len(inter) > 0 else X)
    return U

# ----------------------------
# T0 reduction
# ----------------------------

def t0_reduce_by_equal_minopens(U):
    n = len(U)
    mp = {}
    classes = []
    class_of = [-1]*n
    for x in range(n):
        key = U[x]
        if key not in mp:
            mp[key] = len(classes)
            classes.append([])
        cid = mp[key]
        classes[cid].append(x)
        class_of[x] = cid

    repU = [None]*len(classes)
    for key,cid in mp.items():
        repU[cid] = key
    return classes, class_of, repU

def compute_class_minopens(classes, class_of, U_rep):
    return [frozenset({class_of[x] for x in U_rep[c]}) for c in range(len(classes))]

def maximal_elements_by_minopens(Ubar):
    q = len(Ubar)
    is_max = [True]*q
    for c in range(q):
        for d in range(q):
            if c == d:
                continue
            if Ubar[c].issubset(Ubar[d]) and Ubar[c] != Ubar[d]:
                is_max[c] = False
                break
    return [c for c in range(q) if is_max[c]]

# ----------------------------
# Crosscut 2-skeleton on maxima
# ----------------------------

def crosscut_2skeleton_from_maxima(maxima, Ubar):
    V = sorted(list(maxima))
    edges = []
    triangles = []

    for i in range(len(V)):
        for j in range(i+1, len(V)):
            a, b = V[i], V[j]
            if len(Ubar[a].intersection(Ubar[b])) > 0:
                edges.append((a,b))

    for i in range(len(V)):
        for j in range(i+1, len(V)):
            for k in range(j+1, len(V)):
                a, b, c = V[i], V[j], V[k]
                if len(Ubar[a].intersection(Ubar[b]).intersection(Ubar[c])) > 0:
                    triangles.append((a,b,c))

    S_count = len(V) + len(edges) + len(triangles)
    return V, edges, triangles, S_count

def betti_0_1_from_2skeleton(V, edges, triangles):
    if len(V) == 0:
        return 0, 0

    facets = []
    for t in triangles:
        facets.append(list(t))

    tri_edges = set()
    for (a,b,c) in triangles:
        tri_edges.add(tuple(sorted((a,b))))
        tri_edges.add(tuple(sorted((a,c))))
        tri_edges.add(tuple(sorted((b,c))))
    for e in edges:
        ee = tuple(sorted(e))
        if ee not in tri_edges:
            facets.append(list(ee))

    verts_in_edges = set()
    for (a,b) in edges:
        verts_in_edges.add(a); verts_in_edges.add(b)
    for v in V:
        if v not in verts_in_edges:
            facets.append([v])

    K = SimplicialComplex(facets)
    H = K.homology(reduced=False, base_ring=QQ)
    b0 = int(H[0].rank()) if 0 in H else 0
    b1 = int(H[1].rank()) if 1 in H else 0
    return b0, b1

# ----------------------------
# Continuity check
# ----------------------------

def continuity_ok(prev_U, curr_U, max_examples=3):
    if prev_U is None or curr_U is None:
        return True, 0, []
    n = min(len(prev_U), len(curr_U))
    bad = 0
    ex = []
    for x in range(n):
        union_prev = set()
        for y in curr_U[x]:
            union_prev.update(prev_U[y])
        union_prev = frozenset(union_prev)
        if union_prev != curr_U[x]:
            bad += 1
            if len(ex) < int(max_examples):
                ex.append((x, len(curr_U[x]), len(union_prev)))
    return (bad == 0), bad, ex


def _stats_minmedmax(sizes):
    if len(sizes) == 0:
        return (0,0,0)
    ss = sorted([int(s) for s in sizes])
    return (ss[0], ss[len(ss)//2], ss[-1])

def _rep_point_of_class(classes, cid):
    return int(min(classes[cid]))

def _pairwise_rep_distances(D, reps):
    arr = []
    for i in range(len(reps)):
        for j in range(i+1, len(reps)):
            arr.append(dist(D, reps[i], reps[j]))
    arr.sort()
    return arr

def _eps_from_rep_distances(D, reps, eps_q):
    arr = _pairwise_rep_distances(D, reps)
    if len(arr) == 0:
        return RR(0)
    return _quantile_value(arr, float(eps_q))

def merge_maximals_generators(D, U, merge_k=1, eps_q=0.10):
    classes, class_of, U_rep = t0_reduce_by_equal_minopens(U)
    Ubar = compute_class_minopens(classes, class_of, U_rep)
    maxima = maximal_elements_by_minopens(Ubar)

    if len(maxima) <= 1:
        return [], 0

    reps = [_rep_point_of_class(classes, cid) for cid in maxima]
    eps = _eps_from_rep_distances(D, reps, eps_q)

    cand = []
    for i in range(len(maxima)):
        for j in range(i+1, len(maxima)):
            ai = maxima[i]; bj = maxima[j]
            di = dist(D, reps[i], reps[j])
            cand.append((di, ai, bj))
    cand.sort(key=lambda t: (t[0], t[1], t[2]))

    merged = []
    used = set()
    merged_count = 0
    for (di, a, b) in cand:
        if merged_count >= int(merge_k):
            break
        if di > eps:
            break
        if a in used or b in used:
            continue
        used.add(a); used.add(b)
        merged.append(frozenset(set(U_rep[a]).union(set(U_rep[b]))))
        merged_count += 1

    return merged, merged_count

def pipeline_step_density_anchor_functorial(
    D, rho, r_k,
    tau,
    prev_U,
    step_index=2,
    k_density=7,
    lambda_scale=1.5,
    include_X=False,
    tiny_k=0,
    uncovered_k_every=3,
    uncovered_k_max=6,
    do_saturation=True,
    saturation_passes=1,
    anchor_cap_set=None,
    do_merge_maximals=False,
    merge_k=1,
    eps_q=0.10,
    keep_generators=True,
):
    n = n_points(D)
    X = all_points_set(n)

    anchors = anchors_at_threshold(rho, tau)
    if anchor_cap_set is not None:
        anchors = anchors.intersection(set(anchor_cap_set))
    anchor_list = sorted(list(anchors))
    gens_anchor = []
    for a in anchor_list:
        gens_anchor.append(anchor_neighborhood(D, a, r_k[a], lambda_scale=lambda_scale))

    if prev_U is None:
        prev_U = [frozenset([i]) for i in range(n)]

    covered = 0
    uncovered = 0
    raw_Ax = [None]*n
    for x in range(n):
        containing = [G for G in gens_anchor if x in G]
        if len(containing) == 0:
            S = set(prev_U[x])
            uncovered += 1
            if int(tiny_k) > 0:
                every = max(1, int(uncovered_k_every))
                k = int(tiny_k) + max(0, int(step_index) - 2) // every
                k = min(int(uncovered_k_max), k)
                k_eff = k if len(prev_U[x]) == 1 else min(2, k)
                if k_eff > 0:
                    S.update(knn_set(D, x, k_eff))
            raw_Ax[x] = frozenset(S)
        else:
            covered += 1
            Gmax = max(containing, key=lambda G: (len(G), tuple(sorted(G))))
            raw_Ax[x] = frozenset(Gmax)

    V = minimal_basis_from_generators(raw_Ax, n, include_X=bool(include_X))

    W = [None]*n
    for x in range(n):
        uu = set()
        for y in V[x]:
            uu.update(prev_U[y])
        W[x] = frozenset(uu) if len(uu) > 0 else X

    U_new = minimal_basis_from_generators(W, n, include_X=bool(include_X))
    if bool(do_saturation):
        for _ in range(max(0, int(saturation_passes)-1)):
            W2 = []
            for x in range(n):
                uu = set()
                for y in U_new[x]:
                    uu.update(prev_U[y])
                W2.append(frozenset(uu) if len(uu) > 0 else X)
            U_new = minimal_basis_from_generators(W2, n, include_X=bool(include_X))

    merged_count = 0
    if bool(do_merge_maximals):
        merged_gens, merged_count = merge_maximals_generators(
            D, U_new, merge_k=int(merge_k), eps_q=float(eps_q)
        )
        if len(merged_gens) > 0:
            if bool(keep_generators):
                gens2 = list(set([frozenset(G) for G in raw_Ax] + [frozenset(G) for G in merged_gens]))
            else:
                gens2 = list(set([frozenset(G) for G in merged_gens]))
            U_new = minimal_basis_from_generators(gens2, n, include_X=bool(include_X))

            if bool(do_saturation):
                Wm = []
                for x in range(n):
                    uu = set()
                    for y in U_new[x]:
                        uu.update(prev_U[y])
                    Wm.append(frozenset(uu) if len(uu) > 0 else X)
                U_new = minimal_basis_from_generators(Wm, n, include_X=bool(include_X))

    classes, class_of, U_rep = t0_reduce_by_equal_minopens(U_new)
    Ubar = compute_class_minopens(classes, class_of, U_rep)
    maxima = maximal_elements_by_minopens(Ubar)
    Vv, Ee, Tt, S_count = crosscut_2skeleton_from_maxima(maxima, Ubar)
    b0, b1 = betti_0_1_from_2skeleton(Vv, Ee, Tt)

    out = {
        "U": U_new,
        "V": V,
        "tau": RR(tau),
        "#T0": int(len(classes)),
        "#max": int(len(maxima)),
        "beta0": int(b0),
        "beta1": int(b1),
        "avg|Ux|": RR(sum([len(u) for u in U_new]))/RR(n) if n>0 else RR(0),
        "avg|Vx|": RR(sum([len(v) for v in V]))/RR(n) if n>0 else RR(0),
        "#anchors": int(len(anchors)),
        "covered": int(covered),
        "uncovered": int(uncovered),
        "#V": int(len(Vv)),
        "#E": int(len(Ee)),
        "#T": int(len(Tt)),
        "#S": int(S_count),
        "merged": int(merged_count),
        "A_stats": _stats_minmedmax([len(raw_Ax[i]) for i in range(n)]),
        "V_stats": _stats_minmedmax([len(V[i]) for i in range(n)]),
    }
    return out


def pipeline_step_join_radius_functorial(
    D,
    radius,
    prev_U,
    step_index=2,
    include_X=False,

    do_saturation=True,
    saturation_passes=1,
):
    n = n_points(D)
    X = all_points_set(n)
    if prev_U is None:
        prev_U = [frozenset([i]) for i in range(n)]

    r = RR(radius)
    raw_Gx = [None]*n
    for x in range(n):
        Sx = [y for y in range(n) if dist(D, x, y) <= r]
        uu = set()
        for y in Sx:
            uu.update(prev_U[y])
        raw_Gx[x] = frozenset(uu) if len(uu) > 0 else X

    V = minimal_basis_from_generators(raw_Gx, n, include_X=bool(include_X))

    W = [None]*n
    for x in range(n):
        uu = set()
        for y in V[x]:
            uu.update(prev_U[y])
        W[x] = frozenset(uu) if len(uu) > 0 else X

    U_new = minimal_basis_from_generators(W, n, include_X=bool(include_X))
    if bool(do_saturation):
        for _ in range(max(0, int(saturation_passes)-1)):
            W2 = []
            for x in range(n):
                uu = set()
                for y in U_new[x]:
                    uu.update(prev_U[y])
                W2.append(frozenset(uu) if len(uu) > 0 else X)
            U_new = minimal_basis_from_generators(W2, n, include_X=bool(include_X))

    classes, class_of, U_rep = t0_reduce_by_equal_minopens(U_new)
    Ubar = compute_class_minopens(classes, class_of, U_rep)
    maxima = maximal_elements_by_minopens(Ubar)
    Vv, Ee, Tt, S_count = crosscut_2skeleton_from_maxima(maxima, Ubar)
    b0, b1 = betti_0_1_from_2skeleton(Vv, Ee, Tt)

    out = {
        "U": U_new,
        "V": V,
        "tau": None,
        "#T0": int(len(classes)),
        "#max": int(len(maxima)),
        "beta0": int(b0),
        "beta1": int(b1),
        "avg|Ux|": RR(sum([len(u) for u in U_new]))/RR(n) if n>0 else RR(0),
        "avg|Vx|": RR(sum([len(v) for v in V]))/RR(n) if n>0 else RR(0),
        "#anchors": 0,
        "covered": 0,
        "uncovered": 0,
        "#V": int(len(Vv)),
        "#E": int(len(Ee)),
        "#T": int(len(Tt)),
        "#S": int(S_count),
        "merged": 0,
        "radius": RR(radius),
    }
    return out

# ----------------------------
# Run filtration (table)
# ----------------------------

def run_density_filtration(
    D,
    m=25,
    k_density=7,
    lambda_scale=2.2,
    schedule="quantile",              
    tau_q_start=0.02,
    tau_q_end=0.90,
    tau_min=None,
    tau_max=None,
    include_X=False,
    tiny_k=0,
    uncovered_k_every=3,
    uncovered_k_max=6,
    do_saturation=True,
    saturation_passes=1,
    coarsen_mode="none",              
    coarsen_start=10**9,              
    radius_start=0.0,
    radius_end=1.0,
    merge_maximals=False,
    merge_every=1,
    merge_k=1,
    eps_q_start=0.05,
    eps_q_end=0.90,
    keep_generators=True,
    check_continuity=True,
    cont_report_examples=3,
    anchor_cap_topA=None,
    anchor_cap_until_step=None,   
):
    n = n_points(D)
    rho, r_k = point_density_knn_radius(D, k_density=k_density)

    anchor_cap_set = None
    if anchor_cap_topA is not None:
        A = int(anchor_cap_topA)
        if A > 0:
            order = sorted(range(n), key=lambda i: (RR(rho[i]), -RR(i)), reverse=True)
            anchor_cap_set = set(order[:min(A, n)])

    if schedule == "linear":
        if tau_min is None or tau_max is None:
            vals = sorted([RR(x) for x in rho])
            t_start = _quantile_value(vals, float(tau_q_start))
            t_end = _quantile_value(vals, float(tau_q_end))
            tau_min_eff = t_start
            tau_max_eff = t_end
        else:
            tau_min_eff = RR(tau_min)
            tau_max_eff = RR(tau_max)
        taus = taus_linear(tau_min_eff, tau_max_eff, m)
    else:
        tau_mode = "linear_in_q"
        if schedule == "linear_in_tau":
            tau_mode = "linear_in_tau"
        taus = taus_from_rho_quantiles(rho, m=m, q_start=tau_q_start, q_end=tau_q_end, mode=tau_mode)

    print("")
    print("=== Density-anchor filtration ===")
    print("n = {} | m = {} | k_density = {} | lambda = {} | schedule = {} | include_X = {} | tiny_k = {}".format(
        n, m, k_density, RR(lambda_scale), schedule, bool(include_X), int(tiny_k)
    ))
    if anchor_cap_set is not None:
        if anchor_cap_until_step is None:
            print("Anchor cap: global top-A by density, A = {}".format(int(anchor_cap_topA)))
        else:
            print("Anchor cap: global top-A by density, A = {} (active through step i<= {})".format(int(anchor_cap_topA), int(anchor_cap_until_step)))
    if schedule == "linear":
        print("Tau linear: min = {} | max = {}".format(RR(taus[0]), RR(taus[-1])))
    else:
        print("Tau quantiles: start = {} | end = {}".format(RR(tau_q_start), RR(tau_q_end)))
    if str(coarsen_mode) == "join_radius":
        _re = ""
        try:
            _re = str(radius_end).lower()
        except Exception:
            _re = ""
        if (radius_end is None) or (_re == "auto"):
            diam = RR(0)
            for _i in range(n):
                for _j in range(_i+1, n):
                    dij = dist(D, _i, _j)
                    if dij > diam:
                        diam = dij
            radius_end = RR(1.01) * diam

    if str(coarsen_mode) != "none":
        print("Late-stage coarsen: mode = {} | start_i = {} | radius: {} -> {}".format(
            str(coarsen_mode), int(coarsen_start), RR(radius_start), RR(radius_end)
        ))

    if bool(merge_maximals):
        print("Merge maximals: True | every = {} | merge_k = {} | eps_q schedule: {} -> {} | do_saturation = {} | saturation_passes = {}".format(
            int(merge_every), int(merge_k), RR(eps_q_start), RR(eps_q_end), bool(do_saturation), int(saturation_passes)
        ))
    else:
        print("Merge maximals: False | do_saturation = {} | saturation_passes = {}".format(bool(do_saturation), int(saturation_passes)))

    print("")
    header = (
        "i\t"
        "tau\t"
        "#T0\t#max\tbeta0\tbeta1\t"
        "avg|Ux|\tavg|Vx|\t"
        "#anchors\tcovered\tuncovered\t"
        "#V\t#E\t#T\t#S\t"
    )
    if bool(merge_maximals):
        header += "merged\t"
    header += "cont_ok?"
    print(header)
    print("-"*160)

    prev_U = None
    outs = []

    for i in range(1, int(m)+1):
        tau = taus[i-1]

        eps_q = float(eps_q_start)
        if bool(merge_maximals) and int(m) > 1:
            a = float(i-1)/float(m-1)
            eps_q = float(eps_q_start) + a*(float(eps_q_end) - float(eps_q_start))
            eps_q = max(0.0, min(eps_q, 1.0))

        do_merge_now = bool(merge_maximals) and (int(merge_every) > 0) and ((i % int(merge_every)) == 0)

        if i == 1:
            prev_U = [frozenset([x]) for x in range(n)]
            classes, class_of, U_rep = t0_reduce_by_equal_minopens(prev_U)
            line = "{}\t-\t{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
                i,
                len(classes),
                len(classes),
                n, 0,
                1.0, 1.0,
                0, n, 0,
                n, 0, 0, n
            )
            if bool(merge_maximals):
                line += "0\t"
            line += "-"
            print(line)
            outs.append({"U": prev_U, "tau": None})
            continue

        if (str(coarsen_mode) == "join_radius") and (int(i) >= int(coarsen_start)):
            denom = max(1, int(m) - int(coarsen_start))
            a = float(int(i) - int(coarsen_start)) / float(denom)
            a = max(0.0, min(a, 1.0))
            radius = RR(radius_start) + RR(a) * (RR(radius_end) - RR(radius_start))
            res = pipeline_step_join_radius_functorial(
                D,
                radius=radius,
                prev_U=prev_U,
                step_index=i,
                include_X=include_X,
                do_saturation=do_saturation,
                saturation_passes=saturation_passes,
            )
            res["tau"] = RR(radius)
        else:
            cap_now = anchor_cap_set if (anchor_cap_set is not None and (anchor_cap_until_step is None or i <= anchor_cap_until_step) and i >= 2) else None
            res = pipeline_step_density_anchor_functorial(
                D, rho, r_k,
                tau=tau,
                prev_U=prev_U,
                step_index=i,
                k_density=k_density,
                lambda_scale=lambda_scale,
                include_X=include_X,
                tiny_k=tiny_k,
                uncovered_k_every=uncovered_k_every,
                uncovered_k_max=uncovered_k_max,
                do_saturation=do_saturation,
                saturation_passes=saturation_passes,
                anchor_cap_set=cap_now,
                do_merge_maximals=do_merge_now,
                merge_k=merge_k,
                eps_q=eps_q,
                keep_generators=keep_generators,
            )

        curr_U = res["U"]

        cont_ok, bad, ex = continuity_ok(prev_U, curr_U, max_examples=cont_report_examples) if check_continuity else (True, 0, [])
        cont_str = "ok" if cont_ok else "BAD({})".format(int(bad))

        line = (
            "{}\t{:.6f}\t{}\t{}\t{}\t{}\t"
            "{:.2f}\t{:.2f}\t"
            "{}\t{}\t{}\t"
            "{}\t{}\t{}\t{}\t"
        ).format(
            i, float(res["tau"]),
            int(res["#T0"]), int(res["#max"]),
            int(res["beta0"]), int(res["beta1"]),
            float(res["avg|Ux|"]), float(res["avg|Vx|"]),
            int(res["#anchors"]), int(res["covered"]), int(res["uncovered"]),
            int(res["#V"]), int(res["#E"]), int(res["#T"]), int(res["#S"]),
        )
        if bool(merge_maximals):
            line += "{}\t".format(int(res["merged"]))
        line += cont_str
        print(line)

        if check_continuity and (not cont_ok) and len(ex) > 0:
            print("  continuity examples (x,|Ux|,|union_prev|): " + ", ".join(["({}, {}, {})".format(a,b,c) for (a,b,c) in ex]))

        prev_U = curr_U
        outs.append(res)

    return outs
