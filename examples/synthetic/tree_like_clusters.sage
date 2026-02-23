# examples/synthetic/tree_like_clusters.sage

import os
import random
from sage.all import RR, matrix, sqrt, set_random_seed, Graphics, point

load("scripts/api.sage")
name = current_example_name()

# ------------------------------------------------------------
# DATA GENERATION PARAMETERS
# ------------------------------------------------------------
SEED = 0
PTS_PER_NODE = 20
BLOB_STD = RR(0.18)
DX = RR(2.2)   
DY = RR(1.8)   
root = (RR(0), RR(0))
c1 = (RR(-DX), -DY)
c2 = (RR(+DX), -DY)
g11 = (RR(-DX - DX/2), -2*DY)
g12 = (RR(-DX + DX/2), -2*DY)
g21 = (RR(+DX - DX/2), -2*DY)
g22 = (RR(+DX + DX/2), -2*DY)
CENTERS = [root, c1, c2, g11, g12, g21, g22]
NUM_NODES = len(CENTERS)

# ------------------------------------------------------------
# Choose topology mode
# ------------------------------------------------------------
topology_mode = "density_guided" 


# ------------------------------------------------------------
# DENSITY-GUIDED
# ------------------------------------------------------------
params_density_guided = dict(
    m = 25,
    k_density = 7,
    lambda_scale = 1.65,
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



def tree_like_blobs(seed=0, centers=None, pts_per_node=20, blob_std=RR(0.18)):
    if centers is None:
        centers = []
    set_random_seed(int(seed))
    random.seed(int(seed))
    pts = []
    for (cx, cy) in centers:
        for _ in range(int(pts_per_node)):
            x = RR(random.gauss(float(cx), float(blob_std)))
            y = RR(random.gauss(float(cy), float(blob_std)))
            pts.append((x, y))
    return pts


# ------------------------------------------------------------
# MAIN EXECUTION
# ------------------------------------------------------------

points = tree_like_blobs(
    seed=SEED,
    centers=CENTERS,
    pts_per_node=PTS_PER_NODE,
    blob_std=BLOB_STD
)

D = distance_matrix_from_points(points)

# Plot data via viz.sage
png = save_scatter_plot(points, name=name)
print("[figure] saved graph plot to: {}".format(png))


# ------------------------------------------------------------
# Run
# ------------------------------------------------------------
if topology_mode == "density_guided":
    outs = run_filtration_density_guided(D, **params_density_guided) 
else:
    raise ValueError("Unknown topology_mode: %s" % topology_mode)


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
# print_mapper_summary(mapper_out, nX=len(points))
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

bars0 = compute_beta0_barcodes(outs, nX=len(points), start_step=2)
print_beta0_barcodes_summary(bars0, top=20)

svg0 = barcode_outfile(name, 0)
write_enhanced_barcode_svg(
    bars0,
    nX=len(points),
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
        nX=len(points),
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
