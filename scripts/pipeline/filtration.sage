# scripts/pipeline/filtration.sage

from sage.all import RR

def run_filtration_density_guided(D, **params):
    n = n_points(D)
    X = frozenset(range(n))

    m = int(params.get("m", 25))
    k_density = int(params.get("k_density", 12))
    lambda_scale = RR(params.get("lambda_scale", 1.5))

    schedule = str(params.get("schedule", "quantile")).strip().lower()
    tau_q_start = float(params.get("tau_q_start", 0.10))
    tau_q_end   = float(params.get("tau_q_end", 0.98))
    tau_q_exponent = float(params.get("tau_q_exponent", 2.0))
    tau_min = params.get("tau_min", None)
    tau_max = params.get("tau_max", None)

    debug = bool(params.get("debug", False))
    debug_schedule = bool(params.get("debug_schedule", debug))
    debug_Uof_vmap = bool(params.get("debug_Uof_vmap", debug))
    debug_poset_map = bool(params.get("debug_poset_map", debug))


    rho, r_k = point_density_knn_radius(D, k_density=k_density)

    # ------------------------------------------------------------
    # Schedule
    # ------------------------------------------------------------
    if schedule == "quantile":
        taus = taus_from_rho_quantiles(
            rho, m - 1,
            q_start=tau_q_start,
            q_end=tau_q_end,
            mode="linear_in_q",
            q_exponent=tau_q_exponent
        )
        if debug_schedule:
            vals = sorted([RR(x) for x in rho])
            q_list = [tau_q_start, 0.20, 0.50, 0.90, tau_q_end]
            s = ", ".join([f"q={q:.3f}→τ={_quantile_value(vals, q):.6f}" for q in q_list])
            print(f"[schedule_debug] ρ-quantiles: {s}")
    else:
        vals = sorted([RR(x) for x in rho])
        tmin = RR(tau_min) if tau_min is not None else RR(_quantile_value(vals, tau_q_start))
        tmax = RR(tau_max) if tau_max is not None else RR(_quantile_value(vals, tau_q_end))
        taus = taus_linear(tmin, tmax, m - 1)

    # ------------------------------------------------------------
    # Table header
    # ------------------------------------------------------------
    if bool(params.get("print_table", True)):
        print_filtration_header_density_guided(params, n)

    outs = []

    # ------------------------------------------------------------
    # Step 1: discrete topology
    # ------------------------------------------------------------
    U_prev = discrete_minimal_opens(n)
    V_prev = distinct_minimal_opens(U_prev, drop_X=True)

    out = dict(
        step=1,
        tau=None,
        n=n,
        A_sets=[],
        B_sets=[],
        U=U_prev,
        V=V_prev,
        beta0=n,
        beta1=0,
        diagnostics={"anchors": 0, "covered": n, "uncovered": 0, "nV": n},
    )

    outs.append(out)
            
    print_filtration_row_density_guided(out)

    # ------------------------------------------------------------
    # Anchor plan for steps 2..m 
    # ------------------------------------------------------------
    AV_list = build_A_V_over_taus(
        D=D, rho=rho, r_k=r_k, taus=taus, lambda_scale=lambda_scale,
        drop_X_from_V=True
    )

    AV_by_step = {}
    covered_points_by_step = {1: list(range(n))}
    uncovered_points_by_step = {1: []}
    for k, av in enumerate(AV_list):
        step = k + 2
        AV_by_step[step] = av
        covered_points_by_step[step] = av["covered_set"]
        uncovered_points_by_step[step] = av["uncovered_set"]


    # ------------------------------------------------------------
    # Steps 2..m
    # ------------------------------------------------------------
    for i in range(2, m + 1):
        av = AV_by_step[i]
        tau = av["tau"]
        A_sets = av["A_sets"]
        U_A = av["U_A"]
        V_A = av["V_A"]
        covered_points = covered_points_by_step[i]
        uncovered_points = uncovered_points_by_step[i]


        if i == 2:
            B_by_point = {int(x): frozenset([int(x)]) for x in uncovered_points}
            seen = set()
            B_sets = []
            for S in B_by_point.values():
                if S not in seen:
                    seen.add(S)
                    B_sets.append(S)

        else:
            prev_av = AV_by_step[i - 1]
            B_by_point, B_sets = build_B_sets(
                D=D,
                prev_U=U_prev,
                uncovered_points=uncovered_points,
                prev_uncovered_points=uncovered_points_by_step[i - 1],
                prev_UA=prev_av["U_A"],
                curr_VA=av["V_A"],
                curr_covered_set=covered_points,
            )

            
        topo_out = minimal_basis_of_topology_density_guided(
            n=n,
            A_sets=A_sets,
            B_sets=B_sets,
            drop_X_from_V=True
        )
        
        U_i = topo_out["U"]
        V_i = topo_out["V"]
        diag_topo = topo_out.get("diagnostics", {}) or {}

        if debug:
            print(f"[step {i}] built U,V; |V_i|={len(V_i)}; start poset")
            
        # ------------------------------------------------------------
        # Faithfulness: every anchor-only minimal open (distinct) appears among U_i(x)
        # ------------------------------------------------------------
        Uvals = set(frozenset(S) for S in U_i)
        missing = [S for S in (V_A or []) if frozenset(S) not in Uvals]

        if len(missing) > 0:
            ex = [(len(S), sorted(list(S))[:10]) for S in missing[:3]]
            print(f"[faithfulness] step {i}: FAIL missing={len(missing)}/{len(V_A)} examples(size, first10pts)={ex}")

            

        # ------------------------------------------------------------
        # Continuity check
        # ------------------------------------------------------------
        ok, bad, ex = continuity_ok(
            prev_U=U_prev,
            curr_U=U_i,
            max_examples=params.get("continuity_max_examples", 3)
        )

        if not ok:
            print(f"[continuity_ok] step {i}: bad={bad}/{n}, examples={ex}")


        # ------------------------------------------------------------
        # T0 → poset → order complex (2-skeleton) and Betti
        # ------------------------------------------------------------
        P = t0_poset_from_minimal_opens(U_i)            
        K = order_complex_2skeleton(P)

        E = K["E"]
        T = K["T"]

        nV0 = len(K["V"])
        beta0 = beta0_from_graph(nV0, K["E"])
        beta1 = beta1_from_2skeleton(nV0, K["E"], K["T"])

        diagnostics = (av.get("diagnostics", {}) or {}).copy()
        diagnostics.update(diag_topo)


        out = dict(
            step=int(i),
            tau=RR(tau),
            n=n,
            A_sets=A_sets,
            U_A=U_A,
            V_A=V_A,
            B_sets=B_sets,
            U=U_i,
            V=V_i,
            P=P,
            K=K,
            beta0=beta0,
            beta1=beta1,
            diagnostics=diagnostics,
        )
        outs.append(out)

        if bool(params.get("print_table", True)):
            print_filtration_row_density_guided(out)
    
        U_prev = U_i
        V_prev = V_i

    
    # ------------------------------------------------------------
    # Canonical class maps between consecutive steps (for barcodes / functorial tracking)
    # ------------------------------------------------------------
    outs.sort(key=lambda o: int(o.get("step", 0)))

    for t in range(len(outs) - 1):
        oi, oj = outs[t], outs[t + 1]
        if oi.get("P") is None or oj.get("P") is None:
            continue
        oi["P_class_map"] = poset_class_map(oi["P"], oj["P"])
    
    if debug_Uof_vmap:
        sanity_check_Uof_map_from_poset(outs, start_step=2)
    if debug_poset_map:
        sanity_check_poset_map(outs, start_step=2)


    return outs


