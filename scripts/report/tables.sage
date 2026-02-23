# scripts/report/tables.sage

from sage.all import RR


def _fmt6(x, dash="-"):
    if x is None:
        return dash
    try:
        return "{:.6f}".format(float(x))
    except Exception:
        return dash


def _fmt_step(step):
    return "inf" if step is None else str(int(step))


def _fmt_ratio(a, b, dash="-/-"):
    try:
        a = int(a)
        b = int(b)
        return f"{a}/{b}" if b != 0 else f"{a}/0"
    except Exception:
        return dash


def _fmt_minmedmax(stats):
    try:
        return "({}/{}/{})".format(int(stats.get("min", 0)), int(stats.get("med", 0)), int(stats.get("max", 0)))
    except Exception:
        return "(-/-/-)"


def _print_truncation(total, shown):
    if total > shown:
        print(f"... ({total - shown} more)")
        

def _fmt_tau(tau):
    return _fmt6(tau)

# ----------------------------------------------------------------------
# Filtration summary table
# ----------------------------------------------------------------------

def print_filtration_header_density_guided(params, n):
    print("\n=== Density guided filtration ===")
    p = params or {}
    try:
        m = int(p.get("m", 0))
    except Exception:
        m = p.get("m", "-")

    try:
        k_density = int(p.get("k_density", 0))
    except Exception:
        k_density = p.get("k_density", "-")

    lambda_scale = _fmt6(p.get("lambda_scale", None))
    schedule = str(p.get("schedule", "-"))

    tau_q_start = _fmt6(p.get("tau_q_start", None))
    tau_q_end = _fmt6(p.get("tau_q_end", None))
    tau_q_exponent = _fmt6(p.get("tau_q_exponent", None))

    print(
        "m = {} | k_density = {} | lambda_scale = {} | schedule = {} | "
        "tau_q_start = {} | tau_q_end = {} | tau_q_exponent = {}".format(
            m, k_density, lambda_scale, schedule, tau_q_start, tau_q_end, tau_q_exponent
        )
    )
    print("")
    print("  i           τ      β₀      β₁   cyc/vert    cyc/edge   anchors  covered  uncovered      |V|       |E|         |T|")
    print("-" * 115)
    

def print_filtration_row_density_guided(step_out):
    i = int(step_out.get("step", step_out.get("i", 0)))
    tau = step_out.get("tau", None)
    beta0 = int(step_out.get("beta0", 0))
    beta1 = int(step_out.get("beta1", 0))
    diag = step_out.get("diagnostics", {}) or {}

    anchors = int(diag.get("anchors", 0))

    covered = int(diag.get("covered", 0))
    uncovered = int(diag.get("uncovered", 0))

    K = step_out.get("K", step_out.get("K", None))
    nV = int(len(K.get("V", []))) if isinstance(K, dict) else beta0
    nE = int(len(K.get("E", []))) if isinstance(K, dict) else 0
    nT = int(len(K.get("T", []))) if isinstance(K, dict) else 0

    cycle_per_vert = (beta1 / int(max(1, nV)))
    cycle_per_edge = (beta1 / int(max(1, nE)))

    tau_str = _fmt_tau(tau)
    print(
        f"{i:3d}  {tau_str:>10}  {beta0:5d}  {beta1:7d}  "
        f"{cycle_per_vert:9.4f}  "
        f"{cycle_per_edge:9.4f}  "
        f"{anchors:8d}  {covered:8d}  {uncovered:8d}  "
        f"{nV:8d}  {nE:8d}  {nT:9d}"
    )



# ----------------------------------------------------------------------
# Mapper printing (tables only)
# ----------------------------------------------------------------------

def _as_int(x, default=0):
    try:
        if x is None:
            return int(default)
        if callable(x):
            return int(default)
        return int(x)
    except Exception:
        return int(default)

    
def _as_float(x, default=None):
    if x is None:
        return default
    try:
        return float(x)
    except Exception:
        return default


def print_mapper_summary(mapper_out, nX=None):
    print("")
    print("=== Mapper baseline (density lens + eps-graph clustering) ===")

    p = {}
    if isinstance(mapper_out, dict):
        p = mapper_out.get("params", {}) or {}

    eps = None
    if isinstance(mapper_out, dict):
        eps = mapper_out.get("eps", None)

    n_str = "?" if nX is None else str(_as_int(nX, default="?"))

    print("n = {} | k_density = {} | n_intervals = {} | overlap = {:.2f} | eps = {:.6f} | min_cluster_size = {}".format(
        n_str,
        _as_int(p.get("k_density", 0), 0),
        _as_int(p.get("n_intervals", 0), 0),
        _as_float(p.get("overlap", 0.0), 0.0),
        _as_float(eps, 0.0),
        _as_int(p.get("min_cluster_size", 0), 0),
    ))

    nodes = mapper_out.get("nodes", []) if isinstance(mapper_out, dict) else []
    edges = mapper_out.get("edges", []) if isinstance(mapper_out, dict) else []

    print("#nodes = {} | #edges = {} | beta0(G) = {} | beta1(G) = {}".format(
        _as_int(len(nodes), 0),
        _as_int(len(edges), 0),
        _as_int(mapper_out.get("beta0", 0) if isinstance(mapper_out, dict) else 0, 0),
        _as_int(mapper_out.get("beta1", 0) if isinstance(mapper_out, dict) else 0, 0),
    ))

def print_mapper_nodes(mapper_out, max_nodes=25):
    nodes = mapper_out.get("nodes", []) if isinstance(mapper_out, dict) else []
    print("")
    print(f"Mapper nodes (first {int(max_nodes)}):")
    print("  id  intv   rep   |members|")
    print("  " + "-" * 26)
    for nd in nodes[:int(max_nodes)]:
        nid = int(nd.get("id", 0))
        itv = int(nd.get("interval_idx", 0))
        rep = int(nd.get("rep", 0))
        msz = int(len(nd.get("members", [])))
        print(f"{nid:4d}  {itv:4d}  {rep:4d}  {msz:9d}")
    if len(nodes) > int(max_nodes):
        print(f"... ({len(nodes) - int(max_nodes)} more)")

# ----------------------------------------------------------------------
# Filtration vs Mapper comparison printing
# ----------------------------------------------------------------------

def print_sweep_compare(best):
    print("")
    print("=== Sweep compare: best step for Filtration vs Mapper ===")
    print(best)

def print_compare_filtration_vs_mapper(cmp, topk=10):
    if cmp is None:
        print("=== Compare Filtration vs Mapper: no data ===")
        return

    step = cmp.get("step", None)
    tau = cmp.get("tau", None)
    mapper_mode = cmp.get("mapper_mode", "components")

    fil_sizes = cmp.get("filtration_sizes", [])
    map_sizes = cmp.get("mapper_sizes", [])

    print("")
    print("=== Compare Filtration(beta0 @ step {}) vs Mapper {} ===".format(step, mapper_mode))
    if tau is not None:
        try:
            print("tau = {:.6f}".format(float(tau)))
        except Exception:
            pass

    print("Filtration regions: {}  |  Mapper {}: {}".format(
        int(cmp.get("n_filtration", 0)), mapper_mode, int(cmp.get("n_mapper", 0))
    ))

    if fil_sizes:
        print("Top region sizes (Filtration):  ", ", ".join(str(int(x)) for x in fil_sizes))
    if map_sizes:
        print("Top region sizes (Mapper):", ", ".join(str(int(x)) for x in map_sizes))

    # -------- table 1 --------
    bestM = cmp.get("best_mapper_for_filtration", [])
    print("\nBest Mapper matches for top Filtration components (Jaccard):")
    print(" fil_rank  |fil|   bestJ   map_rank  |map|")
    print(" " + "-" * 46)
    for r in bestM[:int(topk)]:
        fr = int(r.get("filtration_rank", 0))
        fs = int(r.get("filtration_size", 0))
        bj = float(r.get("bestJ", 0.0))
        mr = int(r.get("mapper_rank", 0))
        ms = int(r.get("mapper_size", 0))
        print(f"{fr:8d}  {fs:5d}  {bj:7.3f}  {mr:8d}  {ms:5d}")

    # -------- table 2 --------
    bestF = cmp.get("best_filtration_for_mapper", [])
    print("\nBest Filtration matches for top Mapper components (Jaccard):")
    print(" map_rank  |map|   bestJ   fil_rank  |fil|")
    print(" " + "-" * 46)
    for r in bestF[:int(topk)]:
        mr = int(r.get("mapper_rank", 0))
        ms = int(r.get("mapper_size", 0))
        bj = float(r.get("bestJ", 0.0))
        fr = int(r.get("filtration_rank", 0))
        fs = int(r.get("filtration_size", 0))
        print(f"{mr:8d}  {ms:5d}  {bj:7.3f}  {fr:8d}  {fs:5d}")



# ----------------------------------------------------------------------
# Barcode tables
# ----------------------------------------------------------------------

def print_beta0_barcodes_summary(bars, top=20):
    def pers(b):
        bt = float(b.get("birth_time", 0.0))
        dt = b.get("death_time", None)
        return 10**18 if dt is None else float(dt) - bt

    bars2 = sorted(bars, key=lambda b: (-pers(b), float(b.get("birth_time", 0.0))))
    bars2 = bars2[:int(top)]

    print("\n=== β0 barcode (top {}) ===".format(len(bars2)))
    print("rank  birth→death   tau_birth→tau_death      C/q      X%     repcls")
    print("-------------------------------------------------------------------")
    for k, b in enumerate(bars2):
        bs = b.get("birth_step", "")
        ds = b.get("death_step", None)
        ds_s = "inf" if ds is None else str(ds)

        bt = b.get("birth_time", None)
        dt = b.get("death_time", None)
        dt_s = "inf" if dt is None else "{:.3f}".format(float(dt))
        bt_s = "-" if bt is None else "{:.3f}".format(float(bt))

        mb = (b.get("meta", {}) or {}).get("birth", {}) or {}
        q = int(mb.get("q", 0))
        C = int(mb.get("comp_num_classes", 0))
        supp = int(mb.get("suppX_size", 0))
        repcls = int(mb.get("rep_class_size", 0))

        Cq = "{}/{}".format(C, q) if q else "-/-"
        Xpct = (100.0 * supp / float(supp if q==0 else max(supp,1)))  

        print("{:>4}  {:>3}→{:>3}     {:>7}→{:>7}       {:>7}   {:>5}   {:>6}".format(
            k, bs, ds_s, bt_s, dt_s, Cq, supp, repcls
        ))

def print_beta1_barcodes_summary(bars, top=20, min_step_life=1):
    def life_steps(b):
        bs = int(b.get("birth_step", 0))
        ds = b.get("death_step", None)
        return 10**18 if ds is None else int(ds) - bs

    barsF = []
    for b in bars:
        bs = int(b.get("birth_step", 0))
        ds = b.get("death_step", None)
        if ds is None or (int(ds) - bs >= int(min_step_life)):
            barsF.append(b)

    bars2 = sorted(barsF, key=lambda b: (-life_steps(b), int(b.get("birth_step", 0))))
    bars2 = bars2[:int(top)]

    print("\n=== β1 barcode (top {}, min_step_life={}) ===".format(len(bars2), int(min_step_life)))
    print("rank  birth→death   tau_birth→tau_death     q   type      i   |C|  |E|  |P|pref")
    print("-" * 86)

    for k, b in enumerate(bars2):
        bs = int(b.get("birth_step", 0))
        ds = b.get("death_step", None)
        ds_s = "inf" if ds is None else str(int(ds))

        bt = b.get("birth_time", None)
        dt = b.get("death_time", None)
        bt_s = "-" if bt is None else "{:.3f}".format(float(bt))
        dt_s = "inf" if dt is None else "{:.3f}".format(float(dt))

        mb = (b.get("meta", {}) or {}).get("birth", {}) or {}
        q = int(mb.get("q", 0))
        btype = str(mb.get("birth_type", "-"))
        h1i = mb.get("h1_basis_index", None)
        h1i_s = "-" if h1i is None else str(int(h1i))

        C = mb.get("cycle_support_classes", [])
        E = mb.get("cycle_edges", [])
        Pp = mb.get("cycle_support_points_prefix", [])

        nC = len(C) if isinstance(C, list) else 0
        nE = len(E) if isinstance(E, list) else 0
        nP = len(Pp) if isinstance(Pp, list) else 0

        print("{rank:>4d}  {bs:>3d}→{ds:>3}     {bt:>7}→{dt:>7}   {q:>4d}  {typ:<8}  {i:>3}  {nC:>4d}  {nE:>4d}  {nP:>7d}".format(
            rank=int(k),
            bs=int(bs),
            ds=str(ds_s),
            bt=str(bt_s),
            dt=str(dt_s),
            q=int(q),
            typ=str(btype)[:8],
            i=str(h1i_s).rjust(3),
            nC=int(nC),
            nE=int(nE),
            nP=int(nP),
        ))


def print_beta1_global_candidates(bars, nX, top=20, min_ratio=0.15):
    def ratio(b):
        mb = (b.get("meta", {}) or {}).get("birth", {}) or {}
        supp = int(mb.get("suppX_size", 0))
        return float(supp) / float(nX) if nX > 0 else 0.0

    cand = [b for b in bars if ratio(b) >= float(min_ratio)]
    cand.sort(key=lambda b: (-ratio(b), int(b.get("birth_step", 0))))

    cand = cand[:int(top)]
    print("\n=== β1 global candidates (top {}, min_ratio={:.2f}) ===".format(len(cand), float(min_ratio)))
    print("rank  birth→death   V  E   X%   diam   tau_birth")
    print("--------------------------------------------------")
    for k, b in enumerate(cand):
        mb = (b.get("meta", {}) or {}).get("birth", {}) or {}
        V = int(mb.get("cycle_num_vertices", 0))
        E = int(mb.get("cycle_num_edges", 0))
        supp = int(mb.get("suppX_size", 0))
        Xpct = 100.0 * supp / float(nX) if nX > 0 else 0.0
        diam = mb.get("suppX_diameter", None)
        diam_s = "-" if diam is None else "{:.2f}".format(float(diam))
        bs = int(b.get("birth_step", 0))
        ds = b.get("death_step", None)
        ds_s = "inf" if ds is None else str(int(ds))
        bt = float(b.get("birth_time", 0.0))

        print("{:>4}  {:>3}→{:>3}   {:>2} {:>3}  {:>4.0f}%  {:>5}   {:>7.3f}".format(
            k, bs, ds_s, V, E, Xpct, diam_s, bt
        ))


def beta1_lifetime_histogram(bars):
    hist = {}
    for b in bars:
        bs = int(b.get("birth_step", 0))
        ds = b.get("death_step", None)
        life = -1 if ds is None else int(ds) - bs
        hist[life] = hist.get(life, 0) + 1
    for k in sorted(hist.keys()):
        label = "inf" if k == -1 else str(k)
        print("  life_steps = {:>3}: {}".format(label, hist[k]))
