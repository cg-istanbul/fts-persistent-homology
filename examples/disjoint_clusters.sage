import os
import random
from sage.all import RR, matrix, sqrt, set_random_seed, Graphics, point

load('fts_pipeline.sage')

# ----------------------------
# Parameters (data)
# ----------------------------
SEED = 0
NUM_CLUSTERS = 5
PTS_PER_CLUSTER = 32          
CLUSTER_STD = RR(0.18)
CENTER_SPACING = RR(4.0)

# ----------------------------
# Parameters (filtration)
# ----------------------------
m = 25
k_density = 6
lambda_scale = 2.0

schedule = "quantile"
tau_q_start = 0.10
tau_q_end   = 0.98
tau_min = None
tau_max = None

anchor_cap_topA = None        
anchor_cap_until_step = None  

include_X = False
tiny_k = 1

uncovered_k_every = 3
uncovered_k_max   = 8

do_saturation = True
saturation_passes = 1

coarsen_mode  = "none"     
coarsen_start = None
radius_start  = 0.0
radius_end    = None

merge_maximals = False
merge_every = 1
merge_k = 1
eps_q_start = 0.05
eps_q_end   = 0.90
keep_generators = True

check_continuity = True

# ----------------------------
# Geometry helpers
# ----------------------------
def euclidean_distance(p, q):
    return sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2)

def distance_matrix_from_points(points):
    n = len(points)
    D = matrix(RR, n, n)
    for i in range(n):
        D[i,i] = RR(0)
        for j in range(i+1, n):
            d = euclidean_distance(points[i], points[j])
            D[i,j] = d
            D[j,i] = d
    return D

def disjoint_gaussian_clusters(seed=0,
                               num_clusters=5,
                               pts_per_cluster=32,
                               center_spacing=RR(4.0),
                               std=RR(0.18)):

    set_random_seed(int(seed))
    random.seed(int(seed))

    centers = [(center_spacing * RR(t), RR(0.0)) for t in range(num_clusters)]

    pts = []
    for (cx, cy) in centers:
        for _ in range(pts_per_cluster):
            x = RR(random.gauss(float(cx), float(std)))
            y = RR(random.gauss(float(cy), float(std)))
            pts.append((x, y))
    return pts

def save_scatter_plot(points, out_png):
    G = Graphics()
    G += point(points, size=18)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    G.save(out_png)
    print("[plot] saved to: {}".format(out_png))

# ----------------------------
# Main
# ----------------------------
points = disjoint_gaussian_clusters(
    seed=SEED,
    num_clusters=NUM_CLUSTERS,
    pts_per_cluster=PTS_PER_CLUSTER,
    center_spacing=CENTER_SPACING,
    std=CLUSTER_STD
)

print("Total points n = {}".format(len(points)))

out_png = os.path.join(os.path.dirname(__file__), "figures", "disjoint_clusters.png")
save_scatter_plot(points, out_png)

D = distance_matrix_from_points(points)

run_density_filtration(
    D,
    m=m,
    k_density=k_density,
    lambda_scale=lambda_scale,
    schedule=schedule,
    tau_q_start=tau_q_start,
    tau_q_end=tau_q_end,
    tau_min=tau_min,
    tau_max=tau_max,
    anchor_cap_topA=anchor_cap_topA,
    anchor_cap_until_step=anchor_cap_until_step,
    include_X=include_X,
    tiny_k=tiny_k,
    uncovered_k_every=uncovered_k_every,
    uncovered_k_max=uncovered_k_max,
    do_saturation=do_saturation,
    saturation_passes=saturation_passes,
    coarsen_mode=coarsen_mode,
    coarsen_start=coarsen_start,
    radius_start=radius_start,
    radius_end=radius_end,
    merge_maximals=merge_maximals,
    merge_every=merge_every,
    merge_k=merge_k,
    eps_q_start=eps_q_start,
    eps_q_end=eps_q_end,
    keep_generators=keep_generators,
    check_continuity=check_continuity,
)
