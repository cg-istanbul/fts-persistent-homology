# scripts/algebra/homology0.sage

from sage.all import RR

def _uf_init(n):
    return list(range(n)), [0]*n

def _uf_find(p, x):
    while p[x] != x:
        p[x] = p[p[x]]
        x = p[x]
    return x

def _uf_union(p, r, a, b):
    ra = _uf_find(p, a)
    rb = _uf_find(p, b)
    if ra == rb:
        return
    if r[ra] < r[rb]:
        p[ra] = rb
    elif r[ra] > r[rb]:
        p[rb] = ra
    else:
        p[rb] = ra
        r[ra] += 1

def beta0_from_graph(nV, E):
    nV = int(nV)
    p, r = _uf_init(nV)
    for (u, v) in E:
        _uf_union(p, r, int(u), int(v))
    roots = set(_uf_find(p, i) for i in range(nV))
    return int(len(roots))


def _class_components(q, E):
    q = int(q)
    parent = list(range(q))
    size = [1] * q
    for (u, v) in (E or []):
        u = int(u); v = int(v)
        if 0 <= u < q and 0 <= v < q and u != v:
            _uf_union(parent, size, u, v)

    root_to_cid = {}
    comp_id_of_class = [-1] * q
    comps_classes = []

    for a in range(q):
        r = _uf_find(parent, a)
        cid = root_to_cid.get(r, None)
        if cid is None:
            cid = len(comps_classes)
            root_to_cid[r] = cid
            comps_classes.append([])
        comp_id_of_class[a] = cid
        comps_classes[cid].append(a)

    for C in comps_classes:
        C.sort()
    return comp_id_of_class, comps_classes

def K_relabel_by_Uof(P, K):
    U_of = P["U_of"]
    n0 = int(P["n0"])

    old_to_key = [frozenset(U_of[v]) for v in range(n0)]

    key_to_old = {}
    for v in range(n0):
        key = old_to_key[v]
        if key not in key_to_old or v < key_to_old[key]:
            key_to_old[key] = v

    Eg = []
    usedV = set()

    for (a, b) in K.get("E", []):
        ka = old_to_key[int(a)]
        kb = old_to_key[int(b)]
        if ka == kb:
            continue
        if kb < ka:
            ka, kb = kb, ka
        Eg.append((ka, kb))
        usedV.add(ka); usedV.add(kb)

    Eg = sorted(set(Eg), key=lambda e: (len(e[0]) + len(e[1]),
                                       (len(e[0]), tuple(sorted(e[0]))),
                                       (len(e[1]), tuple(sorted(e[1])))))

    Tg = []
    for (a, b, c) in K.get("T", []):
        ka = old_to_key[int(a)]
        kb = old_to_key[int(b)]
        kc = old_to_key[int(c)]
        if len({ka, kb, kc}) < 3:
            continue
        tri = tuple(sorted((ka, kb, kc),
                           key=lambda s: (len(s), tuple(sorted(s)))))
        Tg.append(tri)
        usedV.update(tri)

    Tg = sorted(set(Tg), key=lambda t: (sum(len(x) for x in t),
                                        (len(t[0]), tuple(sorted(t[0]))),
                                        (len(t[1]), tuple(sorted(t[1]))),
                                        (len(t[2]), tuple(sorted(t[2])))))

    Vg = sorted(usedV, key=lambda s: (len(s), tuple(sorted(s))))

    return dict(V=Vg, E=Eg, T=Tg, old_to_key=old_to_key, key_to_old=key_to_old)


def global_vertex_map(P_i, P_j):
    cls_j = P_j["cls_of"]
    U_of_j = P_j["U_of"]

    key_map = {}
    for Ui in P_i["U_of"]:
        Ui = frozenset(Ui)
        if len(Ui) == 0:
            continue
        x = min(Ui)  
        cj = int(cls_j[int(x)])
        Uj = frozenset(U_of_j[cj])
        key_map[Ui] = Uj
    return key_map


def compute_beta0_barcodes(outs, nX, start_step=2, max_support_points=0):
    nX = int(nX)

    steps = []
    for o in outs:
        if int(o.get("step", 0)) < int(start_step):
            continue
        if o.get("P") is None or o.get("K") is None:
            continue
        steps.append(o)
    steps.sort(key=lambda o: int(o["step"]))
    if not steps:
        return []

    def _components_on_class_graph(P, K):
        n0 = int(P["n0"])
        E = K.get("E", [])

        p, r = _uf_init(n0)
        for (a, b) in E:
            _uf_union(p, r, int(a), int(b))

        comps = {}
        for c in range(n0):
            root = _uf_find(p, c)
            comps.setdefault(root, []).append(int(c))
        for root in comps:
            comps[root].sort()

        roots = sorted(comps.keys(), key=lambda x: (len(comps[x]), comps[x]))
        comp_index = {root: i for i, root in enumerate(roots)}
        comp_id_of_class = [0]*n0
        comp_classes = []
        for root in roots:
            i = comp_index[root]
            classes = comps[root]
            comp_classes.append(classes)
            for c in classes:
                comp_id_of_class[c] = i

        return comp_id_of_class, comp_classes

    def _support_points(P, classes):
        members = P["members_of"]
        supp = []
        for c in classes:
            supp.extend(members[int(c)])
        supp = sorted(set(int(x) for x in supp))
        return supp

    def _birth_meta(o, P, comp_classes):
        supp = _support_points(P, comp_classes)
        rep_point = int(supp[0]) if supp else 0

        cls_of = P["cls_of"]
        U_of = P.get("U_of", [])
        rep_class = int(cls_of[rep_point])

        mb = dict(
            q=int(P["n0"]),
            comp_num_classes=int(len(comp_classes)),
            suppX_size=int(len(supp)),
            rep_point=int(rep_point),
            rep_class=int(rep_class),
            rep_class_size=int(len(U_of[rep_class]) if 0 <= rep_class < len(U_of) else 0),
        )
        if int(max_support_points) > 0:
            mb["suppX_points_prefix"] = supp[:int(max_support_points)]
        return mb, supp

    o0 = steps[0]
    P0, K0 = o0["P"], o0["K"]
    tau0 = float(o0["tau"])
    step0 = int(o0["step"])

    comp_id0, comp_classes0 = _components_on_class_graph(P0, K0)

    bars = []
    active = {} 

    for cid, classes in enumerate(comp_classes0):
        mb, supp = _birth_meta(o0, P0, classes)
        rep_point = int(mb["rep_point"])

        bar = dict(
            dim=0,
            birth_step=step0,
            death_step=None,
            birth_time=float(tau0),
            death_time=None,
            meta=dict(birth=mb, end=None),
        )
        active[rep_point] = dict(bar=bar, cid=int(cid), classes=list(classes), supp=list(supp))

    for t in range(len(steps) - 1):
        oi, oj = steps[t], steps[t + 1]
        Pi, Kj = oi["P"], oj["K"]
        Pj = oj["P"]

        tau_j = float(oj["tau"])
        step_j = int(oj["step"])

        f = oi.get("P_class_map", None)
        if f is None:
            f = poset_class_map(Pi, Pj)

        comp_id_j, comp_classes_j = _components_on_class_graph(Pj, Kj)

        target_of = {}
        for rep_point, st in active.items():
            classes_i = st["classes"]
            imgs = [int(f[int(c)]) for c in classes_i]
            targets = {}
            for cc in imgs:
                targets[int(comp_id_j[cc])] = targets.get(int(comp_id_j[cc]), 0) + 1

            if len(targets) == 1:
                target_of[rep_point] = next(iter(targets.keys()))
            else:
                best = sorted(targets.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
                target_of[rep_point] = int(best)
                st["bar"].setdefault("meta", {}).setdefault("warnings", []).append(
                    {"type": "class_component_splits", "step_from": int(oi["step"]), "targets": sorted(list(targets.keys()))}
                )

        bucket = {}
        for rep_point, tgt in target_of.items():
            bucket.setdefault(tgt, []).append(int(rep_point))

        new_active = {}

        for tgt_cid, reps in bucket.items():
            elder = min(reps)  

            for rp in reps:
                if rp == elder:
                    continue
                bar = active[rp]["bar"]
                bar["death_step"] = int(step_j)
                bar["death_time"] = float(tau_j)
                bar["meta"]["end"] = dict(step=int(step_j), tau=float(tau_j), merged_into=int(elder))
                bars.append(bar)

            classes_new = comp_classes_j[int(tgt_cid)]
            mb, supp = _birth_meta(oj, Pj, classes_new)
            st = active[elder]
            st["cid"] = int(tgt_cid)
            st["classes"] = list(classes_new)
            st["supp"] = list(supp)
            new_active[elder] = st

        occupied = set(st["cid"] for st in new_active.values())
        for cid, classes in enumerate(comp_classes_j):
            if int(cid) in occupied:
                continue
            mb, supp = _birth_meta(oj, Pj, classes)
            rp = int(mb["rep_point"])
            if rp in new_active:
                rp = int(10**9 + cid)

            bar = dict(
                dim=0,
                birth_step=int(step_j),
                death_step=None,
                birth_time=float(tau_j),
                death_time=None,
                meta=dict(birth=mb, end=None),
            )
            new_active[rp] = dict(bar=bar, cid=int(cid), classes=list(classes), supp=list(supp))

        active = new_active

    last = steps[-1]
    for rp, st in active.items():
        bar = st["bar"]
        bar["meta"]["end"] = dict(step=int(last["step"]), tau=float(last["tau"]))
        bars.append(bar)

    def _pers(b):
        bt = float(b.get("birth_time", 0.0))
        dt = b.get("death_time", None)
        return 10**18 if dt is None else float(dt) - bt

    bars.sort(key=lambda b: (float(b.get("birth_time", 0.0)), -_pers(b)))
    return bars
