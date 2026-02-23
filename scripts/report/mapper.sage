# scripts/report/mapper.sage

def _linspace(a, b, m):
    a = RR(a); b = RR(b)
    if m <= 1:
        return [a]
    step = (b - a) / RR(m-1)
    return [a + RR(i)*step for i in range(m)]

def _cover_intervals(vals, n_intervals=10, overlap=0.30):
    vals = [RR(v) for v in vals]
    vmin = min(vals); vmax = max(vals)
    if vmin == vmax:
        return [(vmin, vmax)]
    n = int(n_intervals)
    ov = float(overlap)
    width = (float(vmax - vmin) / float(n)) if n > 0 else float(vmax - vmin)
    step = width * (1.0 - ov) if n > 1 else width
    intervals = []
    L = float(vmin)
    for _ in range(n):
        R = L + width
        intervals.append((RR(L), RR(R)))
        L = L + step
    return intervals

def _eps_global_quantile(D, q=0.10):
    n = int(D.nrows())
    ds = []
    for i in range(n):
        for j in range(i+1, n):
            ds.append(float(D[i,j]))
    ds.sort()
    if not ds:
        return RR(0)
    idx = int(round(float(q) * float(len(ds)-1)))
    idx = max(0, min(idx, len(ds)-1))
    return RR(ds[idx])

def _components_eps_graph(D, subset, eps):
    subset = list(subset)
    if not subset:
        return []
    S = set(subset)
    eps = RR(eps)
    seen = set()
    comps = []
    for v in subset:
        if v in seen:
            continue
        stack = [v]
        seen.add(v)
        comp = [v]
        while stack:
            a = stack.pop()
            for b in subset:
                if b in seen:
                    continue
                if D[a,b] <= eps:
                    seen.add(b)
                    stack.append(b)
                    comp.append(b)
        comps.append(comp)
    return comps

def run_mapper_density_baseline(
    D,
    k_density=7,
    n_intervals=10,
    overlap=0.30,
    eps_q=0.10,
    eps=None,
    min_cluster_size=2,
):
    n = int(D.nrows())

    k = int(k_density)
    lens = []
    for i in range(n):
        row = [float(D[i,j]) for j in range(n) if j != i]
        row.sort()
        kk = max(1, min(k, len(row)))
        lens.append(RR(row[kk-1]) if row else RR(0))

    intervals = _cover_intervals(lens, n_intervals=int(n_intervals), overlap=float(overlap))

    if eps is None:
        eps = _eps_global_quantile(D, q=float(eps_q))
    eps = RR(eps)

    nodes = []
    node_id = 0
    for t, (L, R) in enumerate(intervals):
        pull = [i for i in range(n) if (lens[i] >= L and lens[i] <= R)]
        comps = _components_eps_graph(D, pull, eps)
        for C in comps:
            if len(C) < int(min_cluster_size):
                continue
            nodes.append({
                "id": int(node_id),
                "interval_idx": int(t),
                "members": frozenset(C),
                "rep": int(min(C)),
                "L": RR(L),
                "R": RR(R),
            })
            node_id += 1

    edges = []
    for i in range(len(nodes)):
        Ai = nodes[i]["members"]
        for j in range(i+1, len(nodes)):
            Aj = nodes[j]["members"]
            if Ai.intersection(Aj):
                edges.append((int(nodes[i]["id"]), int(nodes[j]["id"])))

    # Graph invariants (no plotting)
    try:
        G = Graph()
        G.add_vertices([nd["id"] for nd in nodes])
        G.add_edges(edges)
        b0 = int(G.connected_components_number()) if G.num_verts() > 0 else 0
        b1 = int(G.num_edges() - G.num_verts() + b0) if G.num_verts() > 0 else 0
    except Exception:
        G = None
        b0 = 0
        b1 = 0

    return {
        "lens": lens,
        "intervals": intervals,
        "eps": eps,
        "eps_q": float(eps_q),
        "nodes": nodes,
        "edges": edges,
        "G": G,               
        "beta0": int(b0),
        "beta1": int(b1),
        "params": {
            "k_density": int(k_density),
            "n_intervals": int(n_intervals),
            "overlap": float(overlap),
            "min_cluster_size": int(min_cluster_size),
        }
    }


def _jaccard(A, B):
    A = set(A); B = set(B)
    if not A and not B:
        return 1.0
    if not A or not B:
        return 0.0
    return float(len(A & B)) / float(len(A | B))

def _best_matches(A, regions, k=1):
    out = []
    for j, B in enumerate(regions):
        out.append((_jaccard(A, B), int(j), int(len(B))))
    out.sort(reverse=True, key=lambda t: t[0])
    return out[:int(k)] if int(k) > 0 else []

def _get_out_at_step(outs, step):
    s = int(step)
    idx = s - 1
    if outs is None or idx < 0 or idx >= len(outs):
        return None
    return outs[idx]

def _poset_from_out(out):
    if not isinstance(out, dict):
        return None
    if isinstance(out.get("P", None), dict):
        return out["P"]
    if isinstance(out.get("poset", None), dict):
        return out["poset"]
    K = out.get("K", None)
    if isinstance(K, dict):
        P = K.get("poset", None)
        if isinstance(P, dict):
            return P
        P = K.get("P", None)
        if isinstance(P, dict):
            return P
    return None

def _t0_classes_from_U(U):
    n = len(U)
    key_to_cid = {}
    cls_of = [-1] * n
    U_of = []
    for x in range(n):
        ux = frozenset(U[x])
        cid = key_to_cid.get(ux, None)
        if cid is None:
            cid = len(U_of)
            key_to_cid[ux] = cid
            U_of.append(ux)
        cls_of[x] = cid
    return cls_of, U_of

def _support_of_component_classes(comp_classes, U_of):
    supp = set()
    for c in comp_classes:
        cc = int(c)
        if 0 <= cc < len(U_of):
            supp |= set(U_of[cc])
    return supp

def _components_from_edges(nV, E):
    nV = int(nV)
    adj = [set() for _ in range(nV)]
    for (a, b) in (E or []):
        a = int(a); b = int(b)
        if 0 <= a < nV and 0 <= b < nV and a != b:
            adj[a].add(b)
            adj[b].add(a)

    seen = [False]*nV
    comps = []
    for v in range(nV):
        if seen[v]:
            continue
        stack = [v]
        seen[v] = True
        cc = set([v])
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if not seen[w]:
                    seen[w] = True
                    stack.append(w)
                    cc.add(w)
        comps.append(cc)
    return comps

def filtration_regions_at_step(outs, step=2, topk=None):
    out = _get_out_at_step(outs, step)
    if out is None:
        return []

    K = out.get("K", None)
    if not isinstance(K, dict):
        return []

    V = K.get("V", None)
    E = K.get("E", None)

    if V is None:
        return []
    nV = len(V)

    U = out.get("U", None)
    if U is None:
        P = _poset_from_out(out)
        if isinstance(P, dict) and "U_of" in P:
            U_of = list(P["U_of"])
        else:
            return []
    else:
        _, U_of = _t0_classes_from_U(U)

    comps = _components_from_edges(nV, E)

    regions = []
    for comp in comps:
        supp = _support_of_component_classes(comp, U_of)
        if supp:
            regions.append(supp)

    regions.sort(key=lambda S: len(S), reverse=True)
    if topk is not None:
        regions = regions[:int(topk)]
    return regions

# ------------------------------------------------------------
# Mapper regions 
# ------------------------------------------------------------

def mapper_regions(mapper_out, topk=None):
    nodes = mapper_out.get("nodes", []) if isinstance(mapper_out, dict) else []
    regions = [set(nd.get("members", [])) for nd in nodes]
    regions = [S for S in regions if len(S) > 0]
    regions.sort(key=lambda S: len(S), reverse=True)
    if topk is not None:
        regions = regions[:int(topk)]
    return regions

def mapper_component_regions(mapper_out, topk=None):
    nodes = mapper_out.get("nodes", []) if isinstance(mapper_out, dict) else []
    edges = mapper_out.get("edges", []) if isinstance(mapper_out, dict) else []

    node_ids = [int(nd.get("id", 0)) for nd in nodes]
    adj = {u: set() for u in node_ids}
    for (a, b) in edges:
        a = int(a); b = int(b)
        if a in adj and b in adj:
            adj[a].add(b)
            adj[b].add(a)

    members = {}
    for nd in nodes:
        u = int(nd.get("id", 0))
        members[u] = set(int(x) for x in nd.get("members", []))

    seen = set()
    comps = []
    for u in node_ids:
        if u in seen:
            continue
        stack = [u]
        seen.add(u)
        cc = set([u])
        while stack:
            v = stack.pop()
            for w in adj.get(v, []):
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
                    cc.add(w)
        comps.append(cc)

    regions = []
    for cc in comps:
        U = set()
        for u in cc:
            U |= members.get(u, set())
        if U:
            regions.append(U)

    regions.sort(key=lambda S: len(S), reverse=True)
    if topk is not None:
        regions = regions[:int(topk)]
    return regions

# ------------------------------------------------------------
# Compare / sweep
# ------------------------------------------------------------

def compare_filtration_beta0_vs_mapper(outs, mapper_out, step=2, topk=10, mapper_mode="nodes"):
    filtration = filtration_regions_at_step(outs, step=int(step), topk=None)

    if mapper_mode == "nodes":
        mpr = mapper_regions(mapper_out, topk=None)
    else:
        mpr = mapper_component_regions(mapper_out, topk=None)

    topk = int(topk)

    best_mapper_for_filtration = []
    for r, A in enumerate(filtration[:topk]):
        bm = _best_matches(A, mpr, k=1)
        if bm:
            j, jidx, msz = bm[0]
            best_mapper_for_filtration.append({
                "filtration_rank": int(r),
                "filtration_size": int(len(A)),
                "bestJ": float(j),
                "mapper_rank": int(jidx),
                "mapper_size": int(msz),
            })
        else:
            best_mapper_for_filtration.append({
                "filtration_rank": int(r),
                "filtration_size": int(len(A)),
                "bestJ": 0.0,
                "mapper_rank": -1,
                "mapper_size": 0,
            })

    best_filtration_for_mapper = []
    for r, A in enumerate(mpr[:topk]):
        bm = _best_matches(A, filtration, k=1)
        if bm:
            j, jidx, fsz = bm[0]
            best_filtration_for_mapper.append({
                "map_rank": int(r),
                "map_size": int(len(A)),
                "bestJ": float(j),
                "filtration_rank": int(jidx),
                "filtration_size": int(fsz),
            })
        else:
            best_filtration_for_mapper.append({
                "map_rank": int(r),
                "map_size": int(len(A)),
                "bestJ": 0.0,
                "filtration_rank": -1,
                "filtration_size": 0,
            })

    raw_nodes = int(len(mapper_out.get("nodes", []))) if isinstance(mapper_out, dict) else 0

    return {
        "step": int(step),
        "tau": (outs[int(step)-1].get("tau", None) if 0 <= int(step)-1 < len(outs) else None),
        "filtration_regions_count": int(len(filtration)),
        "mapper_mode": str(mapper_mode),

        "mapper_regions_count": int(len(mpr)),
        "mapper_nodes_count": int(len(mpr)),
        "mapper_raw_nodes_count": raw_nodes,

        "filtration_top_sizes": [int(len(S)) for S in filtration[:topk]],
        "mapper_top_sizes": [int(len(S)) for S in mpr[:topk]],
        "best_mapper_for_filtration": best_mapper_for_filtration,
        "best_filtration_for_mapper": best_filtration_for_mapper,
    }

def sweep_best_step_by_overlap(outs, mapper_out, step_min=2, step_max=None, topk=5, mapper_mode="nodes", min_filtration_size=5):
    if step_max is None:
        step_max = len(outs)

    if mapper_mode == "nodes":
        mpr = mapper_regions(mapper_out, topk=None)
    else:
        mpr = mapper_component_regions(mapper_out, topk=None)

    step_min = int(step_min)
    step_max = int(step_max)
    topk = int(topk)

    best = {"step": None, "tau": None, "score": -1.0, "filtration_sizes": []}

    for step in range(step_min, step_max + 1):
        filtration = filtration_regions_at_step(outs, step=step, topk=None)
        filtration = [S for S in filtration if len(S) >= int(min_filtration_size)]
        if not filtration or not mpr:
            continue

        scores = []
        for A in filtration[:topk]:
            bm = _best_matches(A, mpr, k=1)
            scores.append(bm[0][0] if bm else 0.0)

        while len(scores) < topk:
            scores.append(0.0)

        score = sum(scores) / float(topk)

        if score > best["score"]:
            best = {
                "step": int(step),
                "tau": outs[step-1].get("tau", None),
                "score": float(score),
                "filtration_sizes": [int(len(S)) for S in filtration[:topk]],
            }

    best["mapper_mode"] = str(mapper_mode)
    best["min_filtration_size"] = int(min_filtration_size)
    return best
