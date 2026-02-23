# examples/real/paul15_n1000.sage

from sage.all import set_random_seed, RR
import os
import numpy as np

load("scripts/api.sage")
name = current_example_name()

set_random_seed(0)

# ------------------------------------------------------------
# Load dataset
# ------------------------------------------------------------
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
X_path = os.path.join(ROOT, "data", "paul15", "paul15_n1000.csv")
L_path = os.path.join(ROOT, "data", "paul15", "paul15_n1000_meta.csv")
print("Meta path:", L_path)

X = np.loadtxt(X_path, delimiter=",")
n, d = X.shape

# Optional labels (not needed for filtration)
labels = None
if os.path.exists(L_path):
    try:
        labels = np.loadtxt(L_path, dtype=str, delimiter=",", skiprows=1, usecols=[1])
        labels = labels.flatten()
        print("Loaded labels: {}".format(len(labels)))
    except Exception as e:
        print("Failed to load labels: {}".format(e))

D = distance_matrix_from_array(X)
print("Distance matrix built: {}x{}".format(n, n))


# ------------------------------------------------------------
# Choose topology mode
# ------------------------------------------------------------
topology_mode = "density_guided" 


# ------------------------------------------------------------
# DENSITY-GUIDED
# ------------------------------------------------------------
params_density_guided = dict(
    m = 25,
    k_density = 12,
    lambda_scale = 1.0,
    schedule = "quantile",
    tau_q_start = 0.08,
    tau_q_end = 0.95,
    tau_q_exponent = 2.0,
    tau_min = None,
    tau_max = None,
    print_table = True,
    debug = False,
    debug_schedule = False,
    debug_Uof_vmap = False,
    debug_poset_map = False,
)



# ------------------------------------------------------------
# Run
# ------------------------------------------------------------
if topology_mode == "density_guided":
    outs = run_filtration_density_guided(D, **params_density_guided) 
else:
    raise ValueError("Unknown topology_mode: %s" % topology_mode)


# ------------------------------------------------------------
# Giant vs dust
# ------------------------------------------------------------

def component_sizes_from_K(K):
    # K["V"] should be an iterable of vertex ids (or list of vertices)
    # K["E"] list of edges (u,v)
    V = list(K["V"])
    E = list(K["E"])
    G = Graph()
    G.add_vertices(V)
    G.add_edges(E)
    comps = G.connected_components()
    sizes = sorted([len(C) for C in comps], reverse=True)
    return sizes

def print_giant_dust_timeline(outs, min_step=2):
    print("\n=== Giant/Dust timeline (from K 1-skeleton) ===")
    print(" step  beta0  largest  frac    dust_count  dust_total  dust_max")
    print("--------------------------------------------------------------")
    for out in outs:
        i = out["step"]
        if i < min_step: 
            continue
        sizes = component_sizes_from_K(out["K"])
        if not sizes:
            continue
        largest = sizes[0]
        total = sum(sizes)
        frac = float(largest) / float(total) if total else 0.0
        dust = sizes[1:]
        dust_count = len(dust)
        dust_total = sum(dust)
        dust_max = max(dust) if dust else 0
        print(f"{i:4d}  {out['beta0']:5d}  {largest:7d}  {frac:5.3f}   {dust_count:9d}  {dust_total:9d}  {dust_max:8d}")

def first_step_giant_exceeds(outs, frac_thresh=0.50, min_step=2):
    for out in outs:
        if out["step"] < min_step:
            continue
        sizes = component_sizes_from_K(out["K"])
        if sizes and float(sizes[0]) / float(sum(sizes)) >= frac_thresh:
            return out["step"], sizes[0], sum(sizes)
    return None

# After outs computed:
print_giant_dust_timeline(outs)
print("\nFirst step with giant ≥ 50%:", first_step_giant_exceeds(outs, 0.50))
print("First step with giant ≥ 80%:", first_step_giant_exceeds(outs, 0.80))



# # ------------------------------------------------------------
# # Mapper compute (no print/plot)
# # ------------------------------------------------------------
# mapper_out = run_mapper_density_baseline(
#     D,
#     k_density=params_density_guided["k_density"],
#     n_intervals=14,
#     overlap=0.50,
#     eps_q=0.10,
#     min_cluster_size=3
# )

# # ---------------------------
# # Mapper: print + plot
# # ---------------------------
# print_mapper_summary(mapper_out, nX=len(X))
# print_mapper_nodes(mapper_out, max_nodes=25)

# png = save_mapper_graph_png(mapper_out, name=name)
# print("[mapper] saved graph plot to: {}".format(png))

# # ---------------------------
# # Compare 
# # ---------------------------
# best_comps = sweep_best_step_by_overlap(outs, mapper_out, mapper_mode="components", topk=5)

# if best_comps is None or best_comps.get("step", None) is None:
#     print("=== Sweep compare: no valid step found (filtration too degenerate under current params) ===")
# else:
#     cmp_comps = compare_filtration_beta0_vs_mapper(
#         outs, mapper_out, step=best_comps["step"], mapper_mode="components"
#     )
#     print_sweep_compare(best_comps)
#     print_compare_filtration_vs_mapper(cmp_comps, topk=10)



# ------------------------------------------------------------
# Barcodes (post-processing) + SVG rendering
# ------------------------------------------------------------

bars0 = compute_beta0_barcodes(outs, nX=len(X), start_step=2)
print_beta0_barcodes_summary(bars0, top=20)

svg0 = barcode_outfile(name, 0)
write_enhanced_barcode_svg(
    bars0,
    nX=len(X),
    outfile=svg0,
    title=barcode_title(name, 0),
    max_bars=200,
)
print("[barcode] saved:", svg0)

j0 = barcode_jsonl_path(name, 0)
write_bars_jsonl(bars0, j0)
print("[barcode] wrote:", j0)


# ------------------------------------------------------------
# H1 induced maps (GF(2)) -- standard simplicial maps 
# ------------------------------------------------------------
build_H1_induced_maps_GF2(
    outs,
    start_step=2,
    strict_missing=True,
    debug=False
)


# ------------------------------------------------------------
# β1 barcodes (GF(2))
# ------------------------------------------------------------
bars1 = compute_beta1_barcodes_GF2(outs, start_step=2, min_step_life=1, strict_missing=True, debug=False)
print_beta1_barcodes_summary(bars1, top=20, min_step_life=1)


if len(bars1) > 0:
    svg1 = barcode_outfile(name, 1)
    write_enhanced_barcode_svg(
        bars1,
        nX=len(X),
        outfile=svg1,
        title=barcode_title(name, 1),
        max_bars=200,
    )
    print("[beta1] saved:", svg1)

    j1 = barcode_jsonl_path(name, 1)
    write_bars_jsonl(bars1, j1)
    print("[beta1] wrote:", j1)
else:
    print("[beta1] empty")
