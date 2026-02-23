# examples/synthetic/six_bumps_on_ring_band.sage

import os
import random
from sage.all import RR, matrix, sqrt, cos, sin, pi, set_random_seed, Graphics, point

load("scripts/api.sage")
name = current_example_name()

# ------------------------------------------------------------
# DATA GENERATION PARAMETERS
# ------------------------------------------------------------
SEED = 0
RING_RADIUS   = RR(3.0)     
RADIAL_STD    = RR(0.18)    
ANGULAR_JITTER_BG = RR(0.0) 
N_BACKGROUND  = 72          
NUM_BUMPS     = 6
PTS_PER_BUMP  = 20         
BUMP_ANG_STD  = RR(0.10)    
BUMP_RAD_STD  = RR(0.12)    


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




def polar_to_cart(r, ang):
    return (r*cos(ang), r*sin(ang))

def ring_band_with_bumps(seed=0):
    set_random_seed(int(seed))
    random.seed(int(seed))

    pts = []

    for t in range(N_BACKGROUND):
        ang = RR(2*pi*t/N_BACKGROUND)
        if ANGULAR_JITTER_BG != 0:
            ang += RR(random.gauss(0.0, float(ANGULAR_JITTER_BG)))
        r = RR(random.gauss(float(RING_RADIUS), float(RADIAL_STD)))
        pts.append(polar_to_cart(r, ang))

    centers = [RR(2*pi*j/NUM_BUMPS) for j in range(NUM_BUMPS)]
    for a0 in centers:
        for _ in range(PTS_PER_BUMP):
            ang = RR(random.gauss(float(a0), float(BUMP_ANG_STD)))
            r   = RR(random.gauss(float(RING_RADIUS), float(BUMP_RAD_STD)))
            pts.append(polar_to_cart(r, ang))

    return pts


# ------------------------------------------------------------
# MAIN EXECUTION
# ------------------------------------------------------------

points = ring_band_with_bumps(seed=SEED)

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
