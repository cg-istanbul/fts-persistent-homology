# scripts/algebra/homology1.sage

from sage.all import GF, matrix

# ---------------------------
# Small helpers
# ---------------------------

def _edge(u, v):
    u = int(u); v = int(v)
    return (u, v) if u < v else (v, u)

def _tri(a, b, c):
    return tuple(sorted((int(a), int(b), int(c))))

def _edge_list(K):
    E = sorted(set(_edge(a, b) for (a, b) in (K.get("E", []) or [])))
    return E

def _edge_index(K):
    E = _edge_list(K)
    idx = {e: i for i, e in enumerate(E)}
    return E, idx

def _tri_list(K):
    T = sorted(set(_tri(a, b, c) for (a, b, c) in (K.get("T", []) or [])))
    return T

def _tri_index(K):
    T = _tri_list(K)
    idx = {t: i for i, t in enumerate(T)}
    return T, idx

def _mat_from_cols_GF2(cols, nrows):
    F = GF(2)
    k = len(cols)
    if k == 0:
        return matrix(F, nrows, 0, 0)

    return matrix(F, [c.list() for c in cols]).transpose()

def _gf2_rank_bitrows(rows):
    basis = {}   
    rank = 0
    for r in rows:
        r = int(r)
        while r:
            b = r.bit_length() - 1
            if b in basis:
                r ^^= basis[b]
            else:
                basis[b] = r
                rank += 1
                break
    return int(rank)

def beta1_from_2skeleton(nV, E, T):
    nV = int(nV)
    E = list(E)
    T = list(T)

    m1 = len(E)  
    m2 = len(T)  

    b0 = beta0_from_graph(nV, E)
    rank_d1 = nV - b0  

    dim_ker_d1 = m1 - rank_d1

    if m2 == 0 or m1 == 0:
        return int(dim_ker_d1)

    edge_index = {}
    for j, (u, v) in enumerate(E):
        uu, vv = (int(u), int(v))
        if uu == vv:
            continue
        if uu > vv:
            uu, vv = vv, uu
        edge_index[(uu, vv)] = j

    rows = []
    for (a, b, c) in T:
        aa, bb, cc = sorted((int(a), int(b), int(c)))
        if aa == bb or aa == cc or bb == cc:
            continue

        bits = 0

        for e in ((aa, bb), (aa, cc), (bb, cc)):
            j = edge_index.get(e)
            if j is None:
                continue
            bits ^^= (1 << int(j))

        if bits != 0:
            rows.append(bits)
            
    rank_d2 = _gf2_rank_bitrows(rows)

    b1 = dim_ker_d1 - rank_d2
    return int(b1) if b1 >= 0 else 0  


def beta1_diagnostics_from_2skeleton(nV, E, T):
    nV = int(nV)
    E = list(E)
    T = list(T)

    m1 = len(E)
    m2 = len(T)

    b0 = beta0_from_graph(nV, E)
    rank_d1 = nV - b0
    dim_ker_d1 = m1 - rank_d1

    if m2 == 0 or m1 == 0:
        return {
            "b0": int(b0),
            "rank_d1": int(rank_d1),
            "dim_ker_d1": int(dim_ker_d1),
            "rank_d2": 0,
            "b1": int(dim_ker_d1),
            "m1": int(m1),
            "m2": int(m2),
        }

    edge_index = {}
    for j, (u, v) in enumerate(E):
        uu, vv = int(u), int(v)
        if uu == vv:
            continue
        if uu > vv:
            uu, vv = vv, uu
        edge_index[(uu, vv)] = j

    rows = []
    missing_edge_hits = 0
    for (a, b, c) in T:
        aa, bb, cc = sorted((int(a), int(b), int(c)))
        bits = 0
        for e in [(aa, bb), (aa, cc), (bb, cc)]:
            j = edge_index.get(e, None)
            if j is None:
                missing_edge_hits += 1
                continue
            bits ^^= (1 << j)
        if bits != 0:
            rows.append(bits)

    rank_d2 = _gf2_rank_bitrows(rows)
    b1 = dim_ker_d1 - rank_d2
    if b1 < 0:
        b1 = 0

    return {
        "b0": int(b0),
        "rank_d1": int(rank_d1),
        "dim_ker_d1": int(dim_ker_d1),
        "rank_d2": int(rank_d2),
        "b1": int(b1),
        "m1": int(m1),
        "m2": int(m2),
        "rows": int(len(rows)),
        "missing_edge_hits": int(missing_edge_hits),
    }


# ---------------------------
# Standard simplicial map checks
# ---------------------------

def simplicial_map_diagnostics_2skeleton(Ki, Kj, vmap, max_examples=3):
    Ej = set(_edge_list(Kj))
    Tj = set(_tri_list(Kj))

    missing_edges = 0
    missing_tris = 0

    ex_edges = []
    ex_tris = []

    for (u, v) in (Ki.get("E", []) or []):
        u = int(u); v = int(v)
        fu = int(vmap[u]); fv = int(vmap[v])
        if fu == fv:
            continue  
        eimg = _edge(fu, fv)
        if eimg not in Ej:
            missing_edges += 1
            if len(ex_edges) < int(max_examples):
                ex_edges.append(((u, v), (fu, fv)))

    for (a, b, c) in (Ki.get("T", []) or []):
        a = int(a); b = int(b); c = int(c)
        fa = int(vmap[a]); fb = int(vmap[b]); fc = int(vmap[c])
        imgs = sorted(set([fa, fb, fc]))
        if len(imgs) <= 1:
            continue  
        if len(imgs) == 2:
            eimg = _edge(imgs[0], imgs[1])
            if eimg not in Ej:
                missing_tris += 1
                if len(ex_tris) < int(max_examples):
                    ex_tris.append(((a, b, c), (fa, fb, fc), "needs_edge", eimg))
        else:
            timg = tuple(imgs)
            if timg not in Tj:
                missing_tris += 1
                if len(ex_tris) < int(max_examples):
                    ex_tris.append(((a, b, c), (fa, fb, fc), "needs_tri", timg))

    ok = (missing_edges == 0 and missing_tris == 0)
    diag = dict(
        missing_edges=int(missing_edges),
        missing_tris=int(missing_tris),
        examples_edges=ex_edges,
        examples_tris=ex_tris,
    )
    return ok, diag


# ---------------------------
# Boundary matrices over GF(2)
# ---------------------------

def boundary_d1_GF2(K):
    F = GF(2)
    V = K.get("V", [])
    nV = len(V)
    E, _ = _edge_index(K)
    d1 = matrix(F, nV, len(E), 0)
    for j, (u, v) in enumerate(E):
        d1[int(u), j] = 1
        d1[int(v), j] = 1
    return d1, E


def boundary_d2_GF2(K, E_index=None, strict=True, warn=True, max_examples=5):
    F = GF(2)
    E, eidx = _edge_index(K) if E_index is None else E_index
    T = _tri_list(K)
    d2 = matrix(F, len(E), len(T), 0)

    bad = []
    for j, (a, b, c) in enumerate(T):
        edges = [_edge(a, b), _edge(a, c), _edge(b, c)]
        missing = [e for e in edges if e not in eidx]
        if missing:
            bad.append({"tri": (int(a), int(b), int(c)), "missing_edges": [(int(x), int(y)) for (x,y) in missing]})
            if warn and len(bad) <= int(max_examples):
                try:
                    print("[d2][WARN] triangle missing edges:", bad[-1])
                except Exception:
                    pass
            if strict:
                raise RuntimeError(f"boundary_d2_GF2: triangle missing edges (first): {bad[-1]}")
            continue

        for e in edges:
            d2[eidx[e], j] = 1

    if (not strict) and len(bad) > 0:
        K.setdefault("warnings", []).append({"where": "d2", "bad_triangles": bad[:int(max_examples)], "bad_count": len(bad)})

    return d2, T


# ---------------------------
# Chain maps induced by vmap (standard simplicial)
# ---------------------------

def C1_map_matrix_GF2(Ki, Kj, vmap, strict_missing=True):
    F = GF(2)
    Ei, _ = _edge_index(Ki)
    Ej, idx_j = _edge_index(Kj)
    Ej_set = set(Ej)

    M = matrix(F, len(Ej), len(Ei), 0)
    missing = 0

    for col, (u, v) in enumerate(Ei):
        fu = int(vmap[int(u)])
        fv = int(vmap[int(v)])
        if fu == fv:
            continue  
        eimg = _edge(fu, fv)
        if eimg not in Ej_set:
            missing += 1
            if strict_missing:
                pass
            continue
        row = idx_j[eimg]
        M[row, col] = 1

    diag = dict(
        num_edges_i=int(len(Ei)),
        num_edges_j=int(len(Ej)),
        missing_edges=int(missing),
    )
    return M, diag


# ---------------------------
# H1 basis data and induced maps
# ---------------------------

def H1_basis_data_GF2(K):
    F = GF(2)

    d1, E = boundary_d1_GF2(K)
    d2, T = boundary_d2_GF2(K, E_index=(E, {e:i for i,e in enumerate(E)}), strict=True, warn=True)

    Z = d1.right_kernel().basis_matrix().transpose()  

    if d2.ncols() == 0:
        B = matrix(F, d2.nrows(), 0, 0)
    else:
        B = d2.column_space().basis_matrix().transpose()  # |E| x bdim

    cols_B = [B.column(j) for j in range(B.ncols())]
    cols_Z = [Z.column(j) for j in range(Z.ncols())]

    basis_cols = list(cols_B)
    H_cols = []

    cur = _mat_from_cols_GF2(basis_cols, Z.nrows())
    cur_rank = int(cur.rank())

    for zcol in cols_Z:
        test = cur.augment(_mat_from_cols_GF2([zcol], Z.nrows()))
        r = int(test.rank())
        if r > cur_rank:
            H_cols.append(zcol)
            cur = test
            cur_rank = r

    H = _mat_from_cols_GF2(H_cols, Z.nrows())
    Z_basis = B.augment(H)


    if Z_basis.ncols() == 0:
        pivot_rows = []
    else:
        A = Z_basis.transpose()               
        Aech = A.echelon_form()
        pivot_rows = list(Aech.pivots())      

    return dict(
        field=F,
        d1=d1, d2=d2,
        E=E, T=T,
        B=B, H=H,
        Z_basis=Z_basis,
        pivot_rows=pivot_rows,
        beta1=int(H.ncols()),
        bdim=int(B.ncols()),
        zdim=int(Z.ncols()),
    )


def _coords_in_Z_basis_GF2(Z_basis, pivot_rows, v):
    F = GF(2)
    v = v.column() if hasattr(v, "is_vector") and v.is_vector() else v  
    k = int(Z_basis.ncols())
    if k == 0:
        return matrix(F, 0, 1, [])
    rows = pivot_rows

    if len(rows) != k:
        raise ValueError("pivot_rows size mismatch; cannot solve coordinates reliably.")
    A = Z_basis.matrix_from_rows(rows)        
    b = v.matrix_from_rows(rows)              
    x = A.solve_right(b)                      
    return x


def H1_induced_map_GF2(Ki, Kj, vmap, cache_i=None, cache_j=None, strict_missing=True):
    Fi = GF(2)

    Di = cache_i if cache_i is not None else H1_basis_data_GF2(Ki)
    Dj = cache_j if cache_j is not None else H1_basis_data_GF2(Kj)

    beta1_i = int(Di["beta1"])
    beta1_j = int(Dj["beta1"])

    if beta1_i == 0 or beta1_j == 0:
        return matrix(Fi, beta1_j, beta1_i, 0), dict(beta1_i=beta1_i, beta1_j=beta1_j, missing_edges=0, rank=0)

    f1, d_f1 = C1_map_matrix_GF2(Ki, Kj, vmap, strict_missing=strict_missing)

    H_i = Di["H"]
    Zb_j = Dj["Z_basis"]
    piv_j = Dj["pivot_rows"]
    bdim_j = int(Dj["bdim"])

    M_H = matrix(Fi, beta1_j, beta1_i, 0)

    for col in range(beta1_i):
        h = H_i.column(col)                     
        img = f1 * h                            

        x = _coords_in_Z_basis_GF2(Zb_j, piv_j, img)   
        for r in range(beta1_j):
            M_H[r, col] = x[bdim_j + r, 0]

    diag = dict(
        beta1_i=beta1_i,
        beta1_j=beta1_j,
        missing_edges=int(d_f1.get("missing_edges", 0)),
        rank=int(M_H.rank()),
        num_edges_i=int(d_f1.get("num_edges_i", 0)),
        num_edges_j=int(d_f1.get("num_edges_j", 0)),
    )
    return M_H, diag


def build_H1_induced_maps_GF2(
    outs,
    start_step=2,
    strict_missing=True,
    debug=False,
    store_matrices=True,
    max_matrix_entries=20000
):
    steps = [o for o in outs
             if int(o.get("step", 0)) >= int(start_step)
             and o.get("P") is not None and o.get("K") is not None]
    steps.sort(key=lambda o: int(o["step"]))

    cache = {}
    for o in steps:
        s = int(o["step"])
        D = H1_basis_data_GF2(o["K"])
        cache[s] = D

        o["H1_cache_summary"] = dict(
            beta1=int(D["beta1"]),
            bdim=int(D["bdim"]),
            zdim=int(D["zdim"]),
            num_edges=int(D["d1"].ncols()),
            num_tris=int(D["d2"].ncols()),
        )

    for t in range(len(steps) - 1):
        oi, oj = steps[t], steps[t + 1]
        si, sj = int(oi["step"]), int(oj["step"])
        vmap = oi.get("P_class_map", None)

        if vmap is None:
            oi["H1_map_diag_to_next"] = dict(
                beta1_i=int(cache[si]["beta1"]),
                beta1_j=int(cache[sj]["beta1"]),
                rank=None,
                invalid=True,
                reason="missing_P_class_map",
                to_step=int(sj),
            )
            if debug:
                print(f"[H1_map][WARN] missing P_class_map at step {si} (to {sj})")
            continue

        ok, sdiag = simplicial_map_diagnostics_2skeleton(oi["K"], oj["K"], vmap, max_examples=3)
        oi["simplicial_diag_to_next"] = sdiag

        if (not ok) and strict_missing:
            oi["H1_map_diag_to_next"] = dict(
                beta1_i=int(cache[si]["beta1"]),
                beta1_j=int(cache[sj]["beta1"]),
                rank=None,
                missing_edges=int(sdiag.get("missing_edges", 0)),
                num_edges_i=int(sdiag.get("num_edges_i", 0)),
                num_edges_j=int(sdiag.get("num_edges_j", 0)),
                invalid=True,
                reason="non_simplicial",
                to_step=int(sj),
            )
            if debug:
                print(f"[simplicial][WARN] step {si}->{sj}: BAD missing_edges={sdiag.get('missing_edges')} missing_tris={sdiag.get('missing_tris')}")
            break  

        M, hdiag = H1_induced_map_GF2(
            oi["K"], oj["K"], vmap,
            cache_i=cache[si],
            cache_j=cache[sj],
            strict_missing=strict_missing
        )

        hdiag = dict(hdiag)
        hdiag["invalid"] = False
        hdiag["to_step"] = int(sj)

        oi["H1_map_diag_to_next"] = hdiag

        bi = int(hdiag.get("beta1_i", 0))
        bj = int(hdiag.get("beta1_j", 0))
        entries = bi * bj

        if store_matrices and entries <= int(max_matrix_entries):
            oi["H1_map_to_next"] = M
            oi["H1_map_diag_to_next"]["stored_matrix"] = True
        else:
            oi["H1_map_to_next"] = None
            oi["H1_map_diag_to_next"]["stored_matrix"] = False
            oi["H1_map_diag_to_next"]["matrix_entries"] = int(entries)

        if debug:
            print(f"[H1_map] {si}->{sj}: beta1 {bi}→{bj}, rank={hdiag.get('rank')}, missing_edges={hdiag.get('missing_edges')}, stored={oi['H1_map_diag_to_next']['stored_matrix']}")

    return outs


def compute_beta1_barcodes_GF2(outs, start_step=2, min_step_life=1, strict_missing=True, debug=False):
    def _warn(o, msg, payload=None):
        try:
            print("[beta1][WARN]", msg, "" if payload is None else str(payload))
        except Exception:
            pass
        if o is not None:
            o.setdefault("warnings", []).append({
                "where": "beta1",
                "msg": str(msg),
                "payload": payload,
            })

    def _collect_step_warnings(o):
        try:
            return list(o.get("warnings", [])) if o is not None else []
        except Exception:
            return []

    def _attach_bar_warnings(meta, birth_o=None, end_o=None):
        if birth_o is not None:
            w = _collect_step_warnings(birth_o)
            if len(w) > 0:
                meta.setdefault("warnings", {})
                meta["warnings"]["birth_step"] = w
        if end_o is not None:
            w = _collect_step_warnings(end_o)
            if len(w) > 0:
                meta.setdefault("warnings", {})
                meta["warnings"]["end_step"] = w


    def _cycle_meta_from_basis_index(step_o, step_s, basis_index, max_points=30, max_edges=2000):
        meta = {}
        
        c = cache.get(int(step_s), None)
        if c is None:
            _warn(step_o, f"Missing H1 cache for step {step_s}; cannot attach cycle metadata.")
            return meta

        E = c.get("E", None)
        H = c.get("H", None)

        if E is None or H is None:
            _warn(step_o, f"H1 cache at step {step_s} lacks E/H; cannot attach cycle metadata.",
                  payload={"has_E": E is not None, "has_H": H is not None})
            return meta

        j = int(basis_index)
        if j < 0 or j >= int(H.ncols()):
            _warn(step_o, f"Requested basis_index out of range at step {step_s}.",
                  payload={"basis_index": j, "beta1": int(H.ncols())})
            return meta

        cycle_edges = []
        support_classes = set()

        mE = int(H.nrows())
        cap = min(mE, int(max_edges))

        for r in range(cap):
            if int(H[r, j]) == 1:
                u, v = E[r]
                u = int(u); v = int(v)
                cycle_edges.append([u, v])
                support_classes.add(u)
                support_classes.add(v)

        if mE > cap:
            _warn(step_o, f"cycle metadata truncated: too many edges ({mE}); capped at {cap}.",
                  payload={"step": int(step_s), "basis_index": j})

        meta["cycle_edges"] = cycle_edges
        meta["cycle_support_classes"] = sorted([int(x) for x in support_classes])

        pts = []
        try:
            members_of = step_o["P"].get("members_of", None)
            if members_of is not None:
                for c_id in meta["cycle_support_classes"]:
                    c_id = int(c_id)

                    if isinstance(members_of, dict):
                        arr = members_of.get(c_id, [])
                    elif isinstance(members_of, (list, tuple)):
                        arr = members_of[c_id] if (0 <= c_id < len(members_of)) else []
                    else:
                        arr = []

                    for p in arr:
                        pts.append(int(p))
                        if len(pts) >= int(max_points):
                            break
                    if len(pts) >= int(max_points):
                        break

        except Exception as e:
            _warn(step_o, "Failed to build cycle_support_points_prefix.", payload={"err": str(e)})

        meta["cycle_support_points_prefix"] = pts
        return meta

    def _cycle_meta_from_h1_vector(step_o, step_s, h1_vec, max_points=30, max_edges=2000):
        meta = {}
        c = cache.get(int(step_s), None)
        if c is None:
            _warn(step_o, f"Missing H1 cache for step {step_s}; cannot attach cycle metadata.")
            return meta
        E = c.get("E", None)
        H = c.get("H", None)
        if E is None or H is None:
            _warn(step_o, f"H1 cache at step {step_s} lacks E/H; cannot attach cycle metadata.",
                  payload={"has_E": E is not None, "has_H": H is not None})
            return meta

        try:
            chain = H * h1_vec   
        except Exception as e:
            _warn(step_o, "Failed to compute chain representative H*h1_vec.", payload={"err": str(e)})
            return meta

        cycle_edges = []
        support_classes = set()

        mE = int(chain.nrows())
        cap = min(mE, int(max_edges))
        for r in range(cap):
            if int(chain[r, 0]) == 1:
                u, v = E[r]
                u = int(u); v = int(v)
                cycle_edges.append([u, v])
                support_classes.add(u); support_classes.add(v)

        if mE > cap:
            _warn(step_o, f"cycle metadata truncated: too many edges ({mE}); capped at {cap}.",
                  payload={"step": int(step_s)})

        meta["cycle_edges"] = cycle_edges
        meta["cycle_support_classes"] = sorted([int(x) for x in support_classes])

        pts = []
        try:
            members_of = step_o["P"].get("members_of", None)
            if members_of is not None:
                for c_id in meta["cycle_support_classes"]:
                    c_id = int(c_id)

                    if isinstance(members_of, dict):
                        arr = members_of.get(c_id, [])
                    elif isinstance(members_of, (list, tuple)):
                        arr = members_of[c_id] if (0 <= c_id < len(members_of)) else []
                    else:
                        arr = []

                    for p in arr:
                        pts.append(int(p))
                        if len(pts) >= int(max_points):
                            break
                    if len(pts) >= int(max_points):
                        break

        except Exception as e:
            _warn(step_o, "Failed to build cycle_support_points_prefix.", payload={"err": str(e)})

        meta["cycle_support_points_prefix"] = pts
        return meta

    F = GF(2)

    steps = [o for o in outs if int(o.get("step",0)) >= int(start_step) and o.get("P") is not None and o.get("K") is not None]
    steps.sort(key=lambda o: int(o["step"]))
    if len(steps) == 0:
        return []

    cache = {}
    for o in steps:
        cache[int(o["step"])] = H1_basis_data_GF2(o["K"])

    active = []   
    bars = []
    next_id = 0

    o0 = steps[0]
    s0 = int(o0["step"])
    t0 = float(o0["tau"])
    b0 = int(cache[s0]["beta1"])

    for j in range(b0):
        v = matrix(F, b0, 1, 0)
        v[j, 0] = 1
        meta_birth = dict(
            q=int(o0["P"]["n0"]),
            birth_type="initial",
            h1_basis_index=int(j),
        )
        meta_birth.update(_cycle_meta_from_basis_index(o0, s0, j))
        w = _collect_step_warnings(o0)
        if len(w) > 0:
            meta_birth["warnings_at_birth_step"] = w

        active.append(dict(
            id=int(next_id),
            birth_step=int(s0),
            birth_time=float(t0),
            vec=v,
            meta_birth=meta_birth,
        ))
        next_id += 1

    def _pivot_columns(A):
        if A.ncols() == 0:
            return []
        Aech = A.echelon_form()
        return sorted([int(c) for c in Aech.pivots()])

    def _complement_basis_vectors(span_cols, dim):
        cols = []
        cur = span_cols
        cur_rank = int(cur.rank()) if cur.ncols() > 0 else 0

        for r in range(int(dim)):
            e = matrix(F, dim, 1, 0)
            e[r, 0] = 1
            test = cur.augment(e) if cur.ncols() > 0 else e
            if int(test.rank()) > cur_rank:
                cols.append(e)
                cur = test
                cur_rank += 1
            if cur_rank == int(dim):
                break
        return cols

    for idx in range(len(steps) - 1):
        oi = steps[idx]
        oj = steps[idx + 1]
        si, sj = int(oi["step"]), int(oj["step"])
        ti, tj = float(oi["tau"]), float(oj["tau"])

        vmap = oi.get("P_class_map", None)
        if vmap is None:
            _warn(oi, f"Missing P_class_map for step {si}->{sj}. Aborting β1 barcode computation.",
                  payload={"step_from": si, "step_to": sj})
            break

        ok_mono, mono_diag = check_poset_map_monotone(oi["P"], oj["P"], vmap, max_examples=5)
        if not ok_mono:
            _warn(oi, f"Poset class map is not order-preserving at {si}->{sj}. Aborting β1 barcode computation.",
                  payload={"step_from": si, "step_to": sj, "diag": mono_diag})
            break

        ok, sdiag = simplicial_map_diagnostics_2skeleton(oi["K"], oj["K"], vmap, max_examples=3)
        if (not ok) and strict_missing:
            _warn(oi, f"Non-simplicial map detected at {si}->{sj}. Aborting β1 barcode computation.",
                  payload={"step_from": si, "step_to": sj, "diag": sdiag})
            break

        pre_diag = oi.get("H1_map_diag_to_next", None)
        pre_M = oi.get("H1_map_to_next", None)

        if isinstance(pre_diag, dict) and pre_diag.get("invalid", False):
            _warn(oi, f"Precomputed H1 map marked invalid at {si}->{sj}. Aborting β1 barcode computation.",
                  payload={"step_from": si, "step_to": sj, "diag": pre_diag})
            break

        if isinstance(pre_diag, dict) and (pre_M is not None):
            M = pre_M
            hdiag = pre_diag
        else:
            M, hdiag = H1_induced_map_GF2(
                oi["K"], oj["K"], vmap,
                cache_i=cache[si], cache_j=cache[sj],
                strict_missing=strict_missing
            )
        bi = int(hdiag["beta1_i"])
        bj = int(hdiag["beta1_j"])

        if bi == 0:
            for a in active:
                bar = dict(
                    dim=1,
                    birth_step=int(a["birth_step"]),
                    death_step=int(sj),
                    birth_time=float(a["birth_time"]),
                    death_time=float(tj),
                    meta=dict(
                        birth=a["meta_birth"],
                        end=dict(step=int(sj), tau=float(tj)),
                    ),
                )
                if (bar["death_step"] - bar["birth_step"]) >= int(min_step_life):
                    bars.append(bar)
            active = []

            for j in range(bj):
                v = matrix(F, bj, 1, 0); v[j,0] = 1
                meta_birth = dict(
                    q=int(oj["P"]["n0"]),
                    birth_type="initial",          
                    h1_basis_index=int(j),
                )
                meta_birth.update(_cycle_meta_from_basis_index(oj, sj, j))
                w = _collect_step_warnings(oj)
                if len(w) > 0:
                    meta_birth["warnings_at_birth_step"] = w
                active.append(dict(id=int(next_id), birth_step=int(sj), birth_time=float(tj), vec=v, meta_birth=meta_birth))
                next_id += 1
            continue

        k = len(active)
        if k != bi:
            _warn(oi,
                  f"H1 basis dimension mismatch at step {si}: active={k} but beta1={bi}. "
                  f"This indicates incorrect bookkeeping. Aborting β1 barcode computation.",
                  payload={"step": si, "active_dim": k, "beta1": bi})
            break


        imgs = []
        for a in active:
            imgs.append(M * a["vec"])   

        A = matrix(F, bj, k, 0)
        for col in range(k):
            A.set_column(col, imgs[col].column(0))

        piv = _pivot_columns(A)
        piv_set = set(piv)

        survivors = []
        deaths = []
        for col, a in enumerate(active):
            if col in piv_set:
                a2 = dict(a)
                a2["vec"] = imgs[col]
                survivors.append(a2)
            else:
                deaths.append(a)

        for a in deaths:
            meta = dict(
                birth=a["meta_birth"],
                end=dict(step=int(sj), tau=float(tj)),
            )
            _attach_bar_warnings(meta, birth_o=None, end_o=oj)  

            bar = dict(
                dim=1,
                birth_step=int(a["birth_step"]),
                death_step=int(sj),
                birth_time=float(a["birth_time"]),
                death_time=float(tj),
                meta=meta,
            )
            if (bar["death_step"] - bar["birth_step"]) >= int(min_step_life):
                bars.append(bar)

        Im = A.matrix_from_columns(piv) if len(piv) > 0 else matrix(F, bj, 0, 0)
        new_cols = _complement_basis_vectors(Im, bj)

        for e in new_cols:
            meta_birth = dict(
                q=int(oj["P"]["n0"]),
                birth_type="cokernel",
                h1_coord_support=[int(r) for r in range(int(bj)) if int(e[r,0]) == 1],
            )
            meta_birth.update(_cycle_meta_from_h1_vector(oj, sj, e))

            w = _collect_step_warnings(oj)
            if len(w) > 0:
                meta_birth["warnings_at_birth_step"] = w

            survivors.append(dict(
                id=int(next_id),
                birth_step=int(sj),
                birth_time=float(tj),
                vec=e,
                meta_birth=meta_birth,
            ))
            next_id += 1

        active = survivors

        if len(active) != bj:
            _warn(oj, f"After update at step {si}->{sj}, active basis size {len(active)} != beta1 {bj}. Aborting.",
                  payload={"step_from": si, "step_to": sj, "active_dim": len(active), "beta1_j": bj})
            break


        if debug:
            print(f"[beta1] {si}->{sj}: bi={bi}, bj={bj}, rank(M)={hdiag['rank']}, deaths={len(deaths)}, births={len(new_cols)}")

    last = steps[-1]
    sl = int(last["step"])
    tl = float(last["tau"])
    for a in active:
        meta = dict(
            birth=a["meta_birth"],
            end=dict(step=int(sl), tau=float(tl)),
        )
        _attach_bar_warnings(meta, birth_o=None, end_o=last)

        bar = dict(
            dim=1,
            birth_step=int(a["birth_step"]),
            death_step=None,
            birth_time=float(a["birth_time"]),
            death_time=None,
            meta=meta,
        )
        bars.append(bar)

    def _pers(b):
        bt = float(b.get("birth_time", 0.0))
        dt = b.get("death_time", None)
        return 10**18 if dt is None else float(dt) - bt

    bars.sort(key=lambda b: (int(b.get("birth_step", 0)), -_pers(b)))
    return bars
