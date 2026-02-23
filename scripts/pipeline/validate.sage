# scripts/pipeline/validate.sage

def sanity_check_Uof_map_from_poset(outs, start_step=2):
    steps = [o for o in outs if int(o.get("step",0)) >= int(start_step) and o.get("P") is not None]
    steps.sort(key=lambda o: int(o["step"]))

    for t in range(len(steps)-1):
        oi, oj = steps[t], steps[t+1]
        si, sj = int(oi["step"]), int(oj["step"])
        Pi, Pj = oi["P"], oj["P"]
        f = oi.get("P_class_map", None)
        if f is None:
            print(f"[Uof_map] missing P_class_map at step {si} (to {sj})")
            continue

        U_i = Pi["U_of"]
        U_j = Pj["U_of"]

        ok = True
        for c in range(int(Pi["n0"])):
            if int(f[c]) < 0 or int(f[c]) >= int(Pj["n0"]):
                ok = False
                break

        if ok:
            print(f"[Uof_map] step {si}->{sj}: OK")
        else:
            print(f"[Uof_map] step {si}->{sj}: BAD")


def sanity_check_poset_map(outs, start_step=2, max_bad=5):
    steps = [o for o in outs if int(o.get("step",0)) >= int(start_step) and o.get("P") is not None]
    steps.sort(key=lambda o: int(o["step"]))
    bad = 0

    for t in range(len(steps)-1):
        oi, oj = steps[t], steps[t+1]
        si, sj = int(oi["step"]), int(oj["step"])
        Pi, Pj = oi["P"], oj["P"]
        f = oi.get("P_class_map", None)
        if f is None:
            print(f"[poset_map] missing at step {si} (to {sj})")
            bad += 1
            if bad >= max_bad: break
            continue

        leq_j = Pj["leq"]
        E_i = Pi.get("E", [])  
        violations = 0

        max_edges_check = oi.get("debug_max_edges_check", None)
        if max_edges_check is not None:
            E_i = E_i[:int(max_edges_check)]

        for (a, b) in E_i:
            a = int(a); b = int(b)
            fa = int(f[a]); fb = int(f[b])
            if not leq_j[fa][fb]:
                violations += 1
                if violations <= 3:
                    print(f"[poset_map] {si}->{sj} violates: {a}<= {b} but f(a)={fa} not <= f(b)={fb}")


        if violations == 0:
            print(f"[poset_map] step {si}->{sj}: OK")
        else:
            print(f"[poset_map] step {si}->{sj}: BAD {violations}")
            bad += 1
            if bad >= max_bad:
                break

    return (bad == 0)

